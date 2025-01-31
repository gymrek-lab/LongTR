#include <iomanip>
#include <iostream>
#include <time.h>

//#include "sys/sysinfo.h"
//#include "sys/types.h"

#include "extract_indels.h"
#include "genotyper_bam_processor.h"

int parseLine(char* line){
  int i = strlen(line);
  while (*line < '0' || *line > '9') line++;
  line[i-3] = '\0';
  i = atoi(line);
  return i;
}

int getUsedPhysicalMemoryKB(){
  FILE* file = fopen("/proc/self/status", "r");
  int result = -1;
  char line[128];

  while (fgets(line, 128, file) != NULL){
    if (strncmp(line, "VmRSS:", 6) == 0){
      result = parseLine(line);
      break;
    }
  }
  fclose(file);
  return result;
}

/*
  Left align BamAlignments in the provided vector and store those that successfully realign in the provided vector.
  Also extracts other information for successfully realigned reads into provided vectors.
 */
void GenotyperBamProcessor::left_align_reads(const RegionGroup& region_group, const std::string& chrom_seq, std::vector<BamAlnList>& alignments,
					     const std::vector< std::vector<double> >& log_p1, const std::vector< std::vector<double> >& log_p2,
					     std::vector< std::vector<double> >& filt_log_p1,  std::vector< std::vector<double> >& filt_log_p2, std::vector<int>& n_p1s, std::vector<int>& n_p2s,
					     std::vector<Alignment>& left_alns){
  locus_left_aln_time_ = clock();
  selective_logger() << "Trimming reads" << std::endl;
  std::map<std::string, int> seq_to_alns;
  int32_t align_fail_count = 0, total_reads = 0;
  left_alns.clear(); filt_log_p1.clear(); filt_log_p2.clear();

  std::vector<bool> passes_region_filters; passes_region_filters.reserve(region_group.num_regions());
  for (unsigned int i = 0; i < alignments.size(); ++i){
    filt_log_p1.push_back(std::vector<double>());
    filt_log_p2.push_back(std::vector<double>());
    int n_p1_sample = 0;
    int n_p2_sample = 0;
    for (unsigned int j = 0; j < alignments[i].size(); ++j, ++total_reads){
      // Discard read if it doesn't entirely overlap with the repeat
      if (alignments[i][j].pos_  > region_group.start() || alignments[i][j].end_pos_ < region_group.stop()) {
        align_fail_count++;
        continue;
      }
      // Trim alignment if it extends very far upstream or downstream of the STR. For tractability, we limit it to 200bp
      alignments[i][j].TrimAlignment((region_group.start() > FLANK_SIZE ? region_group.start()-FLANK_SIZE : 1), region_group.stop()+FLANK_SIZE);
      if (alignments[i][j].Length() == 0){ // if string is deleted, add it as a deleted alignment
        Alignment new_aln(region_group.start(), region_group.stop(), alignments[i][j].IsReverseStrand(), true, alignments[i][j].Name(), "", "", "");
        left_alns.push_back(new_aln);
        filt_log_p1[i].push_back(log_p1[i][j]);
        filt_log_p2[i].push_back(log_p2[i][j]);
        std::vector<bool> region_passes;
        region_passes.push_back(true); // use this deleted alignment for haplotype generation
        left_alns.back().set_hap_gen_info(region_passes);
        continue;
      }
    Alignment new_aln(alignments[i][j].Position(), alignments[i][j].GetEndPosition() - 1,
                    alignments[i][j].IsReverseStrand(), alignments[i][j].GetDeleted(), alignments[i][j].Name(),
                    alignments[i][j].Qualities(), uppercase(alignments[i][j].QueryBases()), "");

    std::string read_sequence = uppercase(alignments[i][j].QueryBases());
    int32_t seq_index = 0, ref_index = alignments[i][j].Position();
    bool soft_clipped = 0;
    std::stringstream aln_ss;
    for (auto cigar_iter = alignments[i][j].CigarData().begin(); cigar_iter != alignments[i][j].CigarData().end(); cigar_iter++){
        int32_t cigar_index    = 0;
        char prev_cigar_type   = '=';
        int32_t prev_cigar_num = 0;

        switch (cigar_iter->Type){
        case 'H':
          break;
        case 'S':
          new_aln.add_cigar_element(CigarElement(cigar_iter->Type, cigar_iter->Length));
          seq_index += cigar_iter->Length;
          soft_clipped = 1;
          break;
        case 'I':
          new_aln.add_cigar_element(CigarElement(cigar_iter->Type, cigar_iter->Length));
          aln_ss << read_sequence.substr(seq_index, cigar_iter->Length);
          seq_index += cigar_iter->Length;
          break;
        case 'D':
          new_aln.add_cigar_element(CigarElement(cigar_iter->Type, cigar_iter->Length));
          aln_ss << std::string(cigar_iter->Length, '-');
          ref_index += cigar_iter->Length;
          break;
        case 'M': case '=': case 'X':
          while (cigar_index < cigar_iter->Length){
        if (read_sequence[seq_index] == static_cast<char>(toupper(chrom_seq[ref_index]))){
          if (prev_cigar_type == '=')
            prev_cigar_num++;
          else {
            if (prev_cigar_num != 0)
              new_aln.add_cigar_element(CigarElement(prev_cigar_type, prev_cigar_num));
            prev_cigar_type = '=';
            prev_cigar_num  = 1;
          }
        }
        else {
          if (prev_cigar_type == 'X')
            prev_cigar_num++;
          else {
            if (prev_cigar_num != 0)
              new_aln.add_cigar_element(CigarElement(prev_cigar_type, prev_cigar_num));
            prev_cigar_type = 'X';
            prev_cigar_num  = 1;
          }
        }
        aln_ss << read_sequence[seq_index];
        cigar_index++; ref_index++; seq_index++;
          }
          if (prev_cigar_num != 0)
        new_aln.add_cigar_element(CigarElement(prev_cigar_type, prev_cigar_num));
          break;
        default:
          printErrorAndDie("Invalid CIGAR option encountered in convertAlignment");
          break;
        }
   }
   new_aln.set_alignment(aln_ss.str());
   if (soft_clipped){
    align_fail_count++;
    continue;
    }

   left_alns.push_back(new_aln);
   seq_to_alns[alignments[i][j].QueryBases()] = left_alns.size()-1;

   int64_t haplotype;
    if (alignments[i][j].HasTag(HAPLOTYPE_TAG.c_str())){
        alignments[i][j].GetIntTag(HAPLOTYPE_TAG.c_str(), haplotype);
        if (haplotype == 1) n_p1_sample++;
        if (haplotype == 2) n_p2_sample++;
    }

    left_alns.back().check_CIGAR_string(); // Ensure alignment is properly formatted
    filt_log_p1[i].push_back(log_p1[i][j]);
    filt_log_p2[i].push_back(log_p2[i][j]);

    passes_region_filters.clear();
    BamProcessor::passes_filters(alignments[i][j], passes_region_filters);
    left_alns.back().set_hap_gen_info(passes_region_filters);
    }
    n_p1s.push_back(n_p1_sample);
    n_p2s.push_back(n_p2_sample);
  }

  locus_left_aln_time_  = (clock() - locus_left_aln_time_)/CLOCKS_PER_SEC;
  total_left_aln_time_ += locus_left_aln_time_;
  if (align_fail_count != 0)
    selective_logger() << "Failed to trim align " << align_fail_count << " out of " << total_reads << " reads" << std::endl;
}

StutterModel* GenotyperBamProcessor::learn_stutter_model(std::vector<BamAlnList>& alignments,
							 const std::vector< std::vector<double> >& log_p1s,
							 const std::vector< std::vector<double> >& log_p2s,
							 bool haploid, const std::vector<std::string>& rg_names, const Region& region){
  std::vector< std::vector<int> > str_bp_lengths(alignments.size());
  std::vector< std::vector<double> > str_log_p1s(alignments.size()), str_log_p2s(alignments.size());
  int inf_reads = 0;
  const int MAX_INF_READS = 10000;

  // Extract bp differences and phasing probabilities for each read if we need to train a stutter model
  for (unsigned int i = 0; i < alignments.size(); ++i){
    for (unsigned int j = 0; j < alignments[i].size(); ++j){
      int bp_diff;
      bool got_size = ExtractCigar(alignments[i][j].CigarData(), alignments[i][j].Position(), region.start()-region.period(), region.stop()+region.period(), bp_diff);
      if (got_size){
	if (bp_diff < -(int)(region.stop()-region.start()+1))
	  continue;
	inf_reads++;
	str_bp_lengths[i].push_back(bp_diff);
	if (log_p1s.size() == 0){
	  str_log_p1s[i].push_back(0); str_log_p2s[i].push_back(0); // Assign equal phasing LLs as no SNP info is available
	}
	else {
	  str_log_p1s[i].push_back(log_p1s[i][j]); str_log_p2s[i].push_back(log_p2s[i][j]);
	}
      }
    }
    if (inf_reads > MAX_INF_READS)
      break;
  }

  if (inf_reads < MIN_TOTAL_READS){
    full_logger() << "Skipping locus with too few informative reads for stutter training: TOTAL=" << inf_reads << ", MIN=" << MIN_TOTAL_READS << std::endl;
    too_few_reads_++;
    return NULL;
  }

  selective_logger() << "Building EM stutter model" << std::endl;
  EMStutterGenotyper length_genotyper(haploid, region.motif(), str_bp_lengths, str_log_p1s, str_log_p2s, rg_names, 0);
  selective_logger() << "Training EM stutter model" << std::endl;
  bool trained = length_genotyper.train(MAX_EM_ITER, ABS_LL_CONVERGE, FRAC_LL_CONVERGE, false, selective_logger());
  if (trained){
    if (output_stutter_models_)
      length_genotyper.get_stutter_model()->write_model(region.chrom(), region.start(), region.stop(), stutter_model_out_);
    num_em_converge_++;
    StutterModel* stutter_model = length_genotyper.get_stutter_model()->copy();
    selective_logger() << "Learned stutter model " << *stutter_model;
    return stutter_model;
  }
  else {
    num_em_fail_++;
    full_logger() << "Stutter model training failed for locus " << region.chrom() << ":" << region.start() << "-" << region.stop()
		  << " with " << inf_reads << " informative reads" << std::endl;
    return NULL;
  }
}

void GenotyperBamProcessor::analyze_reads_and_phasing(std::vector<BamAlnList>& alignments,
						      std::vector< std::vector<double> >& log_p1s,
						      std::vector< std::vector<double> >& log_p2s,
						      const std::vector<std::string>& rg_names, const RegionGroup& region_group, const std::string& chrom_seq){
  int32_t total_reads = 0;
  for (unsigned int i = 0; i < alignments.size(); i++)
    total_reads += alignments[i].size();
  if (total_reads < MIN_TOTAL_READS){
    full_logger() << "Skipping locus with too few reads: TOTAL=" << total_reads << ", MIN=" << MIN_TOTAL_READS << std::endl;
    too_few_reads_++;
    return;
  }
  // Can't simply check the total number of reads because the bam processor may have stopped reading at the threshold and then removed PCR duplicates
  // Instead, we check this flag which it sets when too many reads are encountered during filtering
  if (TOO_MANY_READS){
    full_logger() << "Skipping locus with too many reads: TOTAL=" << total_reads << ", MAX=" << MAX_TOTAL_READS << std::endl;
    too_many_reads_++;
    return;
  }

  assert(alignments.size() == log_p1s.size() && alignments.size() == log_p2s.size() && alignments.size() == rg_names.size());
  bool haploid = (haploid_chroms_.find(region_group.chrom()) != haploid_chroms_.end());
  const std::vector<Region>& regions = region_group.regions();

  // Clip the reads using de Bruijn graph assembly principles
  //assembly_based_read_clipping(alignments, region_group, chrom_seq);

  // Learn the stutter model for each region
  std::vector<StutterModel*> stutter_models;
  locus_stutter_time_  = clock();
  bool stutter_success = true;
  for (auto region_iter = regions.begin(); region_iter != regions.end(); region_iter++){
    StutterModel* stutter_model = NULL;
    if (def_stutter_model_ != NULL){
      stutter_model = def_stutter_model_->copy();
      stutter_model->set_period(region_iter->period());
    }
    else if (read_stutter_models_){
      // Attempt to extact model from dictionary
      auto model_iter = stutter_models_.find(*region_iter);
      if (model_iter != stutter_models_.end())
	stutter_model = model_iter->second->copy();
      else {
	full_logger() << "WARNING: No stutter model found for " << region_iter->chrom() << ":" << region_iter->start() << "-" << region_iter->stop() << std::endl;
	num_missing_models_++;
      }
    }
    else {
      // Learn stutter model using length-based EM algorithm
      stutter_model = learn_stutter_model(alignments, log_p1s, log_p2s, haploid, rg_names, *region_iter);
    }
    stutter_models.push_back(stutter_model);
    stutter_success &= (stutter_model != NULL);
  }
  locus_stutter_time_  = (clock() - locus_stutter_time_)/CLOCKS_PER_SEC;
  total_stutter_time_ += locus_stutter_time_;

  // Genotype the regions, if requested
  locus_genotype_time_ = clock();
  SeqStutterGenotyper* seq_genotyper = NULL;
  if (vcf_writer_.is_open() && stutter_success) {
    std::vector<Alignment> left_alignments;
    std::vector< std::vector<double> > filt_log_p1s, filt_log_p2s;
    std::vector<int> n_p1s, n_p2s;
    left_align_reads(region_group, chrom_seq, alignments, log_p1s, log_p2s, filt_log_p1s,
		     filt_log_p2s, n_p1s, n_p2s, left_alignments);
    bool run_assembly = (REQUIRE_SPANNING == 0);
    seq_genotyper = new SeqStutterGenotyper(region_group, haploid, 1, left_alignments, filt_log_p1s, filt_log_p2s, n_p1s, n_p2s, rg_names, chrom_seq,
					    stutter_models, ref_vcf_, selective_logger(), skip_assembly_, INDEL_FLANK_LEN, SWITCH_OLD_ALIGN_LEN, alignment_parameters_);
    if (seq_genotyper->genotype(MAX_TOTAL_HAPLOTYPES, MAX_FLANK_HAPLOTYPES, MIN_FLANK_FREQ, selective_logger())) {
      bool pass = true;
      // If appropriate, recalculate the stutter model using the haplotype ML alignments,
      // realign the reads and regenotype the samples
//      if (recalc_stutter_model_)
//	pass = seq_genotyper->recompute_stutter_models(selective_logger(), MAX_TOTAL_HAPLOTYPES, MAX_FLANK_HAPLOTYPES, MIN_FLANK_FREQ, MAX_EM_ITER, ABS_LL_CONVERGE, FRAC_LL_CONVERGE);

      if (pass){
	num_genotype_success_++;
	seq_genotyper->write_vcf_record(samples_to_genotype_, chrom_seq, output_viz_, (VIZ_LEFT_ALNS == 1), viz_out_, &vcf_writer_, selective_logger());
      }
      else
	num_genotype_fail_++;
    }
    else
      num_genotype_fail_++;
  }
  locus_genotype_time_  = (clock() - locus_genotype_time_)/CLOCKS_PER_SEC;
  total_genotype_time_ += locus_genotype_time_;

  selective_logger() << "Locus timing:"                                          << "\n"
		     << " BAM seek time       = " << locus_bam_seek_time()       << " seconds\n"
		     << " Read filtering      = " << locus_read_filter_time()    << " seconds\n"
		     << " SNP info extraction = " << locus_snp_phase_info_time() << " seconds\n"
		     << " Stutter estimation  = " << locus_stutter_time()        << " seconds\n";
  if (stutter_success && vcf_writer_.is_open()){
    selective_logger() << " Genotyping          = " << locus_genotype_time()       << " seconds\n";
    if (vcf_writer_.is_open()){
      assert(seq_genotyper != NULL);
      selective_logger() << "\t" << "Trim alignment        = "  << locus_left_aln_time_             << " seconds\n"
			 << "\t" << " Haplotype generation  = "  << seq_genotyper->hap_build_time()  << " seconds\n"
			 << "\t" << " Haplotype alignment   = "  << seq_genotyper->hap_aln_time()    << " seconds\n";
//			 << "\t" << " Flank assembly        = "  << seq_genotyper->assembly_time()   << " seconds\n"
//			 << "\t" << " Posterior computation = "  << seq_genotyper->posterior_time()  << " seconds\n"
//			 << "\t" << " Alignment traceback   = "  << seq_genotyper->aln_trace_time()  << " seconds\n";

      process_timer_.add_time("Trimming alignment",        locus_left_aln_time_);
      process_timer_.add_time("Haplotype generation",  seq_genotyper->hap_build_time());
      process_timer_.add_time("Haplotype alignment",   seq_genotyper->hap_aln_time());
      //process_timer_.add_time("Flank assembly",        seq_genotyper->assembly_time());
      //process_timer_.add_time("Posterior computation", seq_genotyper->posterior_time());
      //process_timer_.add_time("Alignment traceback",   seq_genotyper->aln_trace_time());
    }
  }

  /*
  logger() << "Total memory in use = " << getUsedPhysicalMemoryKB() << " KB"
	   << std::endl;
  */

  full_logger() << "\n";

  delete seq_genotyper;
  for (int i = 0; i < stutter_models.size(); i++)
    delete stutter_models[i];
}
