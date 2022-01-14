#include <sstream>

#include "bam_io.h"
#include "error.h"
#include "stringops.h"

std::string BuildCigarString(const std::vector<CigarOp>& cigar_data){
  std::stringstream cigar_string;
  for (auto cigar_iter = cigar_data.begin(); cigar_iter != cigar_data.end(); cigar_iter++)
    cigar_string << cigar_iter->Length << cigar_iter->Type;
  return cigar_string.str();
}

void BamAlignment::ExtractSequenceFields(){
  int32_t length = b_->core.l_qseq;
  bases_         = std::string(length, ' ');
  qualities_     = std::string(length, ' '); 
  if (length == 0)
    return;
  
  // Rebuild the quality string
  // 33 is the reference point for the quality encoding
  uint8_t* quals = bam_get_qual(b_);
  for (int32_t i = 0; i < length; ++i)
    qualities_[i] = (char)(quals[i] + 33);
  
  // Rebuild the sequenced bases string
  uint8_t *bases = bam_get_seq(b_);
  for (int32_t i = 0; i < length; ++i)
    bases_[i] = HTSLIB_INT_TO_BASE[bam_seqi(bases, i)];
  
  // Rebuild the CIGAR operations
  int32_t num_cigar_ops = b_->core.n_cigar;
  uint32_t* cigars      = bam_get_cigar(b_);
  cigar_ops_.clear();
  for (int32_t i = 0; i < num_cigar_ops; ++i)
    cigar_ops_.push_back(CigarOp(bam_cigar_opchr(cigars[i]), bam_cigar_oplen(cigars[i])));
  
  built_ = true;
}


void BamHeader::parse_read_groups(const char *text){
  assert(read_groups_.empty());
  std::stringstream ss; ss << text;
  std::string line;
  std::vector<std::string> tokens;
  while (std::getline(ss, line)){
    if (string_starts_with(line, "@RG")){
      split_by_delim(line, '\t', tokens);
      ReadGroup rg;
      for (int i = 1; i < tokens.size(); i++){
	if (string_starts_with(tokens[i], "ID:"))
	  rg.SetID(tokens[i].substr(3));
	else if (string_starts_with(tokens[i], "SM:"))
	  rg.SetSample(tokens[i].substr(3));
	else if (string_starts_with(tokens[i], "LB:"))
	  rg.SetLibrary(tokens[i].substr(3));
      }
      tokens.clear();
      read_groups_.push_back(rg);
    }
  }
}

void BamCramReader::clear_cram_data_structures(){
  if (iter_ != NULL && in_->is_cram){
    // Change the region for the CRAM file to be a single position and iterate through all of the records
    // We use this to deallocate the memory created during CRAM IO, as deallocation only
    // appears to happen once the iterator reads all of the associated records for the region
    bam1_t *b = bam_init1();
    cram_range cur_range = in_->fp.cram->range;
    cram_range new_range = {cur_range.refid, cur_range.start, cur_range.start+1};
    in_->fp.cram->range  = new_range;
    while (sam_itr_next(in_, iter_, b) >= 0){}
    bam_destroy1(b);

    hts_itr_destroy(iter_);
    iter_ = NULL;
  }
}

BamCramReader::BamCramReader(const std::string& path, std::string fasta_path)
  : path_(path), chrom_(""){

  // Open the file itself
  if (!file_exists(path))
    printErrorAndDie("File " + path + " doest not exist");
  in_ = sam_open(path.c_str(), "r");
  if (in_ == NULL)
    printErrorAndDie("Failed to open file " + path);

  if (in_->is_cram){
    if (fasta_path.empty())
      printErrorAndDie("Must specify a FASTA reference file path for CRAM file " + path);

    // Open the FASTA reference file for the CRAM
    char* fasta = new char[fasta_path.size()+1];
    for (size_t i = 0; i < fasta_path.size(); ++i)
      fasta[i] = fasta_path[i];
    fasta[fasta_path.size()] = '\0';

    if (cram_load_reference(in_->fp.cram, fasta) < 0)
      printErrorAndDie("Failed to open FASTA reference file for CRAM file");
    delete [] fasta;
  }

  // Read the header
  if ((hdr_ = sam_hdr_read(in_)) == 0)
    printErrorAndDie("Failed to read the header for file " + path);
  header_        = new BamHeader(hdr_);
  shared_header_ = false;

  // Open the index
  idx_ = sam_index_load(in_, path.c_str());
  if (idx_ == NULL)
    printErrorAndDie("Failed to load the index for file " + path);

  iter_            = NULL;
  start_           = -1;
  end_             = -1;
  min_offset_      = 0;
  reuse_first_aln_ = false;
  cram_done_       = false;
}

BamCramReader::~BamCramReader(){
  if (!shared_header_){
    bam_hdr_destroy(hdr_);
    delete header_;
  }

  clear_cram_data_structures();
  hts_idx_destroy(idx_);
  sam_close(in_);

  if (iter_ != NULL)
    hts_itr_destroy(iter_);
}

bool BamCramReader::SetChromosome(const std::string& chrom){
  iter_            = sam_itr_querys(idx_, hdr_, chrom.c_str());
  chrom_           = chrom;
  min_offset_      = 0;
  reuse_first_aln_ = false;

  if (iter_ != NULL){
    start_     = 0;
    end_       = INT32_MAX;
    cram_done_ = false;
    return true;
  }
  else {
    start_     = -1;
    end_       = -1;
    cram_done_ = true;
    return false;
  }
}

bool BamCramReader::SetRegion(const std::string& chrom, int32_t start, int32_t end){
  if (in_->is_cram && iter_ != NULL && chrom.compare(chrom_) == 0 && start >= start_){
    // Determine if we can reuse the CRAM iterator from the previous region
    // and if so, modify the iterator accordingly
    bool can_reuse = false;
    if (min_offset_ != 0 && in_->fp.cram != NULL && in_->fp.cram->ctr != NULL){
      cram_container *ctr = in_->fp.cram->ctr;
      if (ctr->slice != NULL){
	cram_slice *slice = ctr->slice;
	if (slice->ref_start+1 < first_aln_.Position() && slice->ref_end > start){
	  can_reuse     = true;
	  // When the min_offset_ is 0, we set it to 1 as 0 is the unset value. We therefore undo that change here
	  // If the offset really was 1, we'll simply consider 1 more alignment that won't overlap the region
	  slice->curr_rec = (min_offset_ == 1 ? 0 : min_offset_);

	  // Modify the iterator's start coordinate to reflect the current region
	  cram_range cur_range = in_->fp.cram->range;
	  cram_range new_range = {cur_range.refid, start+1, cur_range.end};
	  in_->fp.cram->range  = new_range;
	}
      }
    }

    if (can_reuse){
      chrom_           = chrom;
      start_           = start;
      end_             = end;
      cram_done_       = false;
      min_offset_      = 0;
      reuse_first_aln_ = false;
      return true;
    }

    // If we can't reuse the CRAM iterator, we'll set the region as in all other cases by proceeding to the rest of the function
  }

  if (iter_ != NULL){
    clear_cram_data_structures();
    hts_itr_destroy(iter_);
    iter_ = NULL;
  }

  // For CRAM files, as we try to reuse the iterator for future regions, we cannot set an end coordinate
  // Otherwise, the CRAM container will be freed after reading the last alignment in the current region
  std::stringstream region;
  if (in_->is_cram)
    region << chrom << ":" << start+1;
  else
    region << chrom << ":" << start+1 << "-" << end;

  std::string region_str = region.str();
  iter_ = sam_itr_querys(idx_, hdr_, region_str.c_str());
  if (iter_ != NULL){
    chrom_     = chrom;
    start_     = start;
    end_       = end;
    cram_done_ = false;

    bool reuse_offset = (!in_->is_cram && min_offset_ != 0 && chrom.compare(chrom_) == 0 && start >= start_);
    if (reuse_offset)
      if (iter_->n_off == 1 && min_offset_ >= iter_->off[0].u && min_offset_ <= iter_->off[0].v)
	iter_->off[0].u = min_offset_;

    if (reuse_offset && first_aln_.GetEndPosition() > start && first_aln_.Position() < end){
      // NOTE: min_offset_ remains unchanged, as the first valid alignment is the current first alignment
      reuse_first_aln_ = true;
    }
    else {
      min_offset_      = 0;
      reuse_first_aln_ = false;
    }
    return true;
  }
  else {
    chrom_           = "";
    start_           = -1;
    end_             = -1;
    min_offset_      = 0;
    reuse_first_aln_ = false;
    cram_done_       = true;
    return false;
  }
}

bool BamCramReader::GetNextAlignment(BamAlignment& aln){
  if (iter_ == NULL) return false;
  if (cram_done_)    return false;

  if (reuse_first_aln_){
    reuse_first_aln_ = false;
    aln = first_aln_;
    return true;
  }

  int ret = sam_itr_next(in_, iter_, aln.b_);
  if ((ret < 0) || aln.b_->core.pos > end_+1){
    if (ret < -1)
      printErrorAndDie("Invalid record encountered in " + path_ + ". Please ensure the BAM/CRAM is not truncated and is properly formatted");

    if (in_->is_cram){
      // Don't destroy the iterator, as we may want to reuse it later
      cram_done_ = true;
      return false;
    }
    else {
      hts_itr_destroy(iter_);
      iter_ = NULL;
      return false;
    }
  }

  // Set up alignment instance variables
  aln.built_    = false;
  aln.file_     = path_;
  aln.ref_      = header_->ref_name(aln.b_->core.tid);
  aln.mate_ref_ = header_->ref_name(aln.b_->core.mtid);
  aln.length_   = aln.b_->core.l_qseq;
  aln.pos_      = aln.b_->core.pos;
  aln.end_pos_  = bam_endpos(aln.b_);

  if (min_offset_ == 0){
    if (in_->is_cram){
      assert(in_->fp.cram != NULL && in_->fp.cram->ctr != NULL);
      cram_container *ctr = in_->fp.cram->ctr;
      assert(ctr->slice != NULL and ctr->slice->curr_rec > 0);
      first_aln_  = aln;
      min_offset_ = ctr->slice->curr_rec - 1;
      if (min_offset_ == 0)
	min_offset_++;
    }
    else {
      first_aln_  = aln;
      min_offset_ = iter_->curr_off;
    }
  }

  return true;
}



bool BamCramMultiReader::SetRegion(const std::string& chrom, int32_t start, int32_t end){
  aln_heap_.clear();
  chrom_ = chrom;
  start_ = start;
  end_   = end;
  if (merge_type_ == ORDER_ALNS_BY_POSITION){
    for (int32_t reader_index = 0; reader_index < bam_readers_.size(); reader_index++){
      if (!bam_readers_[reader_index]->SetRegion(chrom, start, end))
	return false;
      if (bam_readers_[reader_index]->GetNextAlignment(cached_alns_[reader_index]))
	aln_heap_.push_back(std::pair<int32_t, int32_t>(-cached_alns_[reader_index].Position(), reader_index));
    }
    reader_unset_ = std::vector<bool>(bam_readers_.size(), false);
  }
  else if (merge_type_ == ORDER_ALNS_BY_FILE){
    for (int32_t reader_index = 0; reader_index < bam_readers_.size(); reader_index++){
      // We avoid doing any region setting here and instead will set it when the first call to GetNextAlignment() requires the reader's data
      // For CRAMs, setting the region for all files will load all of their containers into memory
      aln_heap_.push_back(std::pair<int32_t, int32_t>(-reader_index, reader_index));
    }
    reader_unset_ = std::vector<bool>(bam_readers_.size(), true);
  }
  else
    printErrorAndDie("Invalid merge order in SetRegion()");
  std::make_heap(aln_heap_.begin(), aln_heap_.end());
  return true;
}

bool BamCramMultiReader::GetNextAlignment(BamAlignment& aln){
  if (aln_heap_.empty())
    return false;
  std::pop_heap(aln_heap_.begin(), aln_heap_.end());
  int32_t reader_index = aln_heap_.back().second;
  aln_heap_.pop_back();

  // Sometimes we don't set a reader's region until invoking GetNextAlignment() to reduce memory usage
  if (reader_unset_[reader_index]){
    reader_unset_[reader_index] = false;
    if (!bam_readers_[reader_index]->SetRegion(chrom_, start_, end_))
      return false;
    if (!bam_readers_[reader_index]->GetNextAlignment(cached_alns_[reader_index]))
      return GetNextAlignment(aln);
  }

  // Assign optimal alignment to provided reference
  aln = cached_alns_[reader_index];

  // Add reader's next alignment to the cache
  if (bam_readers_[reader_index]->GetNextAlignment(cached_alns_[reader_index])){
    if (merge_type_ == ORDER_ALNS_BY_POSITION)
      aln_heap_.push_back(std::pair<int32_t, int32_t>(-cached_alns_[reader_index].Position(), reader_index));
    else if (merge_type_ == ORDER_ALNS_BY_FILE)
      aln_heap_.push_back(std::pair<int32_t, int32_t>(-reader_index, reader_index));
    else
      printErrorAndDie("Invalid merge order in GetNextAlignment()");

    std::push_heap(aln_heap_.begin(), aln_heap_.end());
  }
  return true;
}


void compare_bam_headers(const BamHeader* hdr_a, const BamHeader* hdr_b, const std::string& file_a, const std::string& file_b){
  std::stringstream error_msg;
  if (hdr_a->num_seqs() != hdr_b->num_seqs()){
    error_msg << "BAM header mismatch issue. BAM headers for files " << file_a << " and " << file_b << " must have the same number of reference sequences";
    printErrorAndDie(error_msg.str());
  }

  for (int32_t i = 0; i < hdr_a->num_seqs(); ++i){
    if (hdr_a->ref_name(i).compare(hdr_b->ref_name(i)) != 0){
      error_msg << "BAM header mismatch issue. Order of reference sequences in BAM headers for files " << file_a << " and " << file_b << " must match";
      printErrorAndDie(error_msg.str());
    }
    if (hdr_a->ref_length(i) != hdr_b->ref_length(i)){
      error_msg << "BAM header mismatch issue. Length of reference sequences in BAM headers for files " << file_a << " and " << file_b << " must match";
      printErrorAndDie(error_msg.str());
    }
  }
}


void BamAlignment::TrimAlignment(int32_t min_read_start, int32_t max_read_stop, char min_base_qual){
  if (!built_)
    ExtractSequenceFields();
  assert(bases_.size() == qualities_.size());

  int ltrim = 0;
  int32_t start_pos = pos_;
  while ((start_pos < min_read_start) && cigar_ops_.size() > 0){
    // Check if we should stop trimming b/c the quality score is above the threshold
    bool qual_above_thresh = false;
    switch (cigar_ops_.front().Type){
    case 'M': case '=': case 'X': case 'I': case 'S':
      qual_above_thresh = (qualities_[ltrim] > min_base_qual);
      break;
    case 'D': case 'H':
      break;
    default:
      printErrorAndDie("Invalid CIGAR option encountered in trimAlignment");
      break;
    }
    if (qual_above_thresh)
      break;

    switch(cigar_ops_.front().Type){
    case 'M': case '=': case 'X':
      ltrim++;
      start_pos++;
      break;
    case 'D':
      start_pos++;
      break;
    case 'I': case 'S':
      ltrim++;
      break;
    case 'H':
      break;
    default:
      printErrorAndDie("Invalid CIGAR option encountered in TrimAlignment");
      break;
    }
    if (cigar_ops_.front().Length == 1)
      cigar_ops_.erase(cigar_ops_.begin(), cigar_ops_.begin()+1);
    else
      cigar_ops_.front().Length--;
  }

  int rtrim = 0, qual_string_len = qualities_.size()-1;
  int32_t end_pos = end_pos_;
  while ((end_pos > max_read_stop) && cigar_ops_.size() > 0){
    // Check if we should stop trimming b/c the quality score is above the threshold
    bool qual_above_thresh = false;
    switch(cigar_ops_.back().Type){
    case 'M': case '=': case 'X': case 'I': case 'S':
      qual_above_thresh = (qualities_[qual_string_len-rtrim] > min_base_qual);
      break;
    case 'D': case 'H':
      break;
    default:
      printErrorAndDie("Invalid CIGAR option encountered in TrimAlignment");
      break;
    }
    if (qual_above_thresh)
      break;

    switch(cigar_ops_.back().Type){
    case 'M': case '=': case 'X':
      rtrim++;
      end_pos--;
      break;
    case 'D':
      end_pos--;
      break;
    case 'I': case 'S':
      rtrim++;
      break;
    case 'H':
      break;
    default:
      printErrorAndDie("Invalid CIGAR option encountered in trimAlignment");
      break;
    }
    if (cigar_ops_.back().Length == 1)
      cigar_ops_.pop_back();
    else
      cigar_ops_.back().Length--;
  }

  assert(ltrim+rtrim <= bases_.size());
  bases_     = bases_.substr(ltrim, bases_.size()-ltrim-rtrim);
  qualities_ = qualities_.substr(ltrim, qualities_.size()-ltrim-rtrim);
  length_   -= (ltrim + rtrim);
  pos_       = start_pos;
  end_pos_   = end_pos;
}

void BamAlignment::TrimLowQualityEnds(char min_base_qual){
  return TrimAlignment(end_pos_+1, pos_-1, min_base_qual);
}
