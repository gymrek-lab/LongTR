#include <algorithm>
#include <climits>
#include <map>
#include <set>
#include <sstream>
#include <iostream>
#include "AlignmentModel.h"
#include "AlignmentTraceback.h"
#include "HapAligner.h"
#include "HapBlock.h"
#include "../mathops.h"
#include "RepeatBlock.h"
#include "StutterAlignerClass.h"
#include "HaplotypeGenerator.h"


// Minimum distance of a seed base from an indel, mismatch or a repetitive region
const int32_t MIN_SEED_DIST = 5;

// Large negative value to prevent impossible or undesirable configurations 
const double IMPOSSIBLE = -1000000000;

// Only consider a base as a SNP if the minimum log-probability that its sequence is correct
// is above this threshold
const double MIN_SNP_LOG_PROB_CORRECT = -0.0043648054;


void HapAligner::align_seq_to_hap(Haplotype* haplotype, bool reuse_alns,
				  const char* seq_0, int seq_len, double& left_prob, const double* base_log_wrong, const double* base_log_correct){
  // NOTE: Input matrix structure: Row = Haplotype position, Column = Read index



  double* L_log_probs = new double[seq_len];

  // Initialize first row of matrix (each base position matched with leftmost haplotype base)
  left_prob = 0.0;
  char first_hap_base = haplotype->get_first_char();
  std::string read_seq = seq_0;
  if (haplotype->get_seq().size() <= 60){ //TODO, it usually happens in case of big deletions
    left_prob = IMPOSSIBLE;
    return;
  }
  int REF_FLANK_LEN = 35; //from HaplotypeGenerator.h
  std::string haplotype_seq = haplotype->get_seq().substr(REF_FLANK_LEN - 5, haplotype->get_seq().size() - (REF_FLANK_LEN - 5)*2);
  int n = haplotype_seq.size();
  int m = read_seq.size();
  double match_matrix_[n][m];
  double deletion_matrix_[n][m];
  double insert_matrix_[n][m];

  int homopolymer_len = 1;

  float MISMATCH = -9.0;
  float MATCH = -0.000100005;

  //double coefficient = 0.5;
  double coefficient = 1.0;
  for (int i = 0; i < m; ++i){
    match_matrix_[0][i]    = (read_seq[i] == haplotype_seq[0] ? MATCH : MISMATCH) + left_prob;
    insert_matrix_[0][i]   = IMPOSSIBLE;
    deletion_matrix_[0][i] = coefficient*LOG_MATCH_TO_DEL[homopolymer_len] + left_prob;
    left_prob += coefficient*LOG_DEL_TO_DEL;
  }



  left_prob = 0.0;
  for (int i = 0; i < n; ++i){
    match_matrix_[i][0]    = (read_seq[0] == haplotype_seq[i] ? MATCH : MISMATCH) + left_prob;
    insert_matrix_[i][0]   = MATCH + coefficient*LOG_MATCH_TO_INS[homopolymer_len] + left_prob;
    deletion_matrix_[i][0] = IMPOSSIBLE;
    left_prob += coefficient*LOG_INS_TO_INS;
  }

  for (int j = 1; j < n; j++){
    for (int i = 1; i < m; i++){

	  double match_emit     = (haplotype_seq[j] == read_seq[i] ? MATCH : MISMATCH);
	  match_matrix_[j][i]    = coefficient*match_emit + std::max(match_matrix_[j-1][i-1] + coefficient*LOG_MATCH_TO_MATCH[homopolymer_len],
	                                       std::max(deletion_matrix_[j-1][i-1] + coefficient*LOG_MATCH_TO_DEL[homopolymer_len],
	                                                  insert_matrix_[j-1][i-1] + coefficient*LOG_MATCH_TO_INS[homopolymer_len]));
	  insert_matrix_[j][i]   = MATCH + std::max(match_matrix_[j][i-1] + coefficient*LOG_INS_TO_MATCH,
									insert_matrix_[j][i-1] + coefficient*LOG_INS_TO_INS);
	  deletion_matrix_[j][i] = std::max(match_matrix_[j-1][i]  + coefficient*LOG_DEL_TO_MATCH,
						   deletion_matrix_[j-1][i] + coefficient*LOG_DEL_TO_DEL);
    }

 }
 left_prob = std::max(deletion_matrix_[n-1][m-1], std::max(insert_matrix_[n-1][m-1], match_matrix_[n-1][m-1]));
 //std::cout << "alignment of sequence with size " << read_seq << " to haplotype with size " << haplotype_seq << " with l_prob " << left_prob << std::endl;
}


void HapAligner::trim_alignment(const Alignment& aln, std::string& trimmed_seq){
    int32_t pos          = aln.get_start();
    int32_t padding = 5;
    int32_t min_read_start = repeat_starts_[0] - padding;
    int32_t max_read_stop = repeat_ends_[0] + padding;
    int32_t start_pos = aln.get_start() + 1;
    int32_t alignment_size = aln.get_sequence().size();
    int32_t end_pos = aln.get_stop() + 1;
    int32_t rtrim = 0;
    int32_t ltrim = 0;

    std::vector<CigarElement> cigar_ops_ = aln.get_cigar_list();

    // left region
    while ((start_pos <= min_read_start) && cigar_ops_.size() > 0){
        switch(cigar_ops_.front().get_type()){
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
              printErrorAndDie("Invalid CIGAR option encountered in trim_alignment");
              break;
        }
        if (cigar_ops_.front().get_num() == 1)
          cigar_ops_.erase(cigar_ops_.begin(), cigar_ops_.begin()+1);
        else
          cigar_ops_.front().set_num(cigar_ops_.front().get_num()-1);
        }

    // left flank
    int mid_pointer = start_pos;
    while ((mid_pointer > min_read_start) && (mid_pointer <= min_read_start + padding) && cigar_ops_.size() > 0){
       switch(cigar_ops_.front().get_type()){
            case 'M': case '=': case 'X':
              mid_pointer++;
              break;
            case 'D':
              ltrim--;
              mid_pointer++;
              break;
            case 'I': case 'S': //If an insertion happened in the flanking region, we will include it
              //ltrim++;
              break;
            case 'H':
              break;
            default:
              printErrorAndDie("Invalid CIGAR option encountered in trim_alignment");
              break;
        }
        if (cigar_ops_.front().get_num() == 1)
          cigar_ops_.erase(cigar_ops_.begin(), cigar_ops_.begin()+1);
        else
          cigar_ops_.front().set_num(cigar_ops_.front().get_num()-1);
    }

    // right region
    while ((end_pos > max_read_stop) && cigar_ops_.size() > 0 && cigar_ops_.size() > 0){
        switch(cigar_ops_.back().get_type()){
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
              printErrorAndDie("Invalid CIGAR option encountered in trim_alignment");
              break;
        }
        if (cigar_ops_.back().get_num() == 1)
          cigar_ops_.pop_back();
        else
          cigar_ops_.back().set_num(cigar_ops_.back().get_num()-1);
    }

    // right flank
    mid_pointer = end_pos;
    while ((mid_pointer > max_read_stop - padding) && (mid_pointer <= max_read_stop) && cigar_ops_.size() > 0){
       switch(cigar_ops_.back().get_type()){
            case 'M': case '=': case 'X':
              mid_pointer--;
              break;
            case 'D':
              rtrim--;
              mid_pointer--;
              break;
            case 'I': case 'S':
              break;
            case 'H':
              break;
            default:
              printErrorAndDie("Invalid CIGAR option encountered in trim_alignment");
              break;
        }
        if (cigar_ops_.back().get_num() == 1)
          cigar_ops_.pop_back();
        else
          cigar_ops_.back().set_num(cigar_ops_.back().get_num()-1);
    }


  if (ltrim < 0) ltrim = 0;
  if (rtrim < 0) rtrim = 0;
  assert(ltrim+rtrim <= alignment_size);
  trimmed_seq     = aln.get_sequence().substr(ltrim, alignment_size-ltrim-rtrim);
  }


void HapAligner::process_reads(const std::vector<Alignment>& alignments, int init_read_index, const BaseQuality* base_quality, const std::vector<bool>& realign_read,
			       double* aln_probs, int* seed_positions){
  assert(alignments.size() == realign_read.size());
  AlignmentTrace trace(fw_haplotype_->num_blocks());
  double* prob_ptr = aln_probs + (init_read_index*fw_haplotype_->num_combs());
  for (unsigned int i = 0; i < alignments.size(); i++){
    if (!realign_read[i]){
      prob_ptr += fw_haplotype_->num_combs();
      continue;
    }

    int seed_base = alignments[i].get_sequence().size() - 1;
    seed_positions[init_read_index+i] = seed_base;
    process_read(alignments[i], seed_base, base_quality, false, prob_ptr, trace);
    prob_ptr += fw_haplotype_->num_combs();

  }
}

const double TRACE_LL_TOL = 0.001;
inline int triple_min_index(double v1, double v2, double v3){
  if (v1 > v2+TRACE_LL_TOL)
    return (v1 > v3+TRACE_LL_TOL ? 0 : 2);
  else
    return (v2 > v3+TRACE_LL_TOL ? 1 : 2);
}

inline int rev_triple_min_index(double v1, double v2, double v3){
  if (v3 > v2+TRACE_LL_TOL)
    return (v3 > v1+TRACE_LL_TOL ? 2 : 0);
  else
    return (v2 > v1+TRACE_LL_TOL ? 1 : 0);
}

inline int     pair_min_index(double v1, double v2){ return (v1 > v2+TRACE_LL_TOL ? 0 : 1); }
inline int rev_pair_min_index(double v1, double v2){ return (v2 > v1+TRACE_LL_TOL ? 1 : 0); }

std::string HapAligner::retrace(Haplotype* haplotype, const char* read_seq, const double* base_log_correct,
				int seq_len, int block_index, int base_index, int matrix_index,
				double* match_matrix, double* insert_matrix, double* deletion_matrix, int* best_artifact_size, int* best_artifact_pos,
				AlignmentTrace& trace){
  const int MATCH = 0, DEL = 1, INS = 2, NONE = -1; // Types of matrices
  int seq_index   = seq_len-1;
  int matrix_type = MATCH;
  std::stringstream aln_ss;

  int (*pair_index_fn)(double, double);
  int (*triple_index_fn)(double, double, double);

  if (!haplotype->reversed()){
    pair_index_fn    = &pair_min_index;
    triple_index_fn  = &triple_min_index;
  }
  else {
    pair_index_fn   = &rev_pair_min_index;
    triple_index_fn = &rev_triple_min_index;
  }

  while (block_index >= 0){
    bool stutter_block = haplotype->get_block(block_index)->get_repeat_info() != NULL;
    if (stutter_block){
      int* artifact_size_ptr = best_artifact_size + seq_len*block_index;
      int* artifact_pos_ptr  = best_artifact_pos  + seq_len*block_index;
      std::stringstream str_ss;
      const std::string& block_seq = haplotype->get_seq(block_index);
      int block_len    = block_seq.size();
      int stutter_size = artifact_size_ptr[seq_index];
      assert(matrix_type == MATCH && base_index+1 == block_len);

      int i = 0;
      for (; i < std::min(seq_index+1, artifact_pos_ptr[seq_index]); i++){
	aln_ss      << "M";
	str_ss      << read_seq[seq_index-i];
      }
      if (artifact_size_ptr[seq_index] < 0)
	aln_ss << std::string(-artifact_size_ptr[seq_index], 'D');
      else
	for (; i < std::min(seq_index+1, artifact_pos_ptr[seq_index] + artifact_size_ptr[seq_index]); i++){
	  aln_ss      << "I";
	  str_ss      << read_seq[seq_index-i];
	}
      for (; i < std::min(block_len + artifact_size_ptr[seq_index], seq_index+1); i++){
	aln_ss      << "M";
	str_ss      << read_seq[seq_index-i];
      }
      std::string str_seq = str_ss.str();

      // Add STR data to trace instance
      if (haplotype->reversed()){
	// Alignment for sequence to right of seed. Block indexes are reversed, but alignment is correct
	trace.add_str_data(haplotype->num_blocks()-1-block_index, stutter_size, str_seq);
      }
      else {
	// Alignment for sequence to left of seed. Block indexes are correct, but alignment is reversed
	std::reverse(str_seq.begin(), str_seq.end());
	trace.add_str_data(block_index, stutter_size, str_seq);
      }

      if (block_len + artifact_size_ptr[seq_index] >= seq_index+1)
	return aln_ss.str(); // Sequence doesn't span stutter block
      else {
	matrix_index -= (block_len + artifact_size_ptr[seq_index] + seq_len*block_len);
	matrix_type   = MATCH;
	seq_index    -= (block_len + artifact_size_ptr[seq_index]);
      }
    }
    else {
      int prev_matrix_type    = NONE;
      std::string block_seq   = haplotype->get_seq(block_index);
      int32_t pos             = haplotype->get_block(block_index)->start() + (haplotype->reversed() ? -base_index : base_index);
      const int32_t increment = (haplotype->reversed() ? 1 : -1);
      int32_t indel_seq_index = -1, indel_position = -1;
      std::stringstream flank_ss;

      // Retrace flanks while tracking any indels that occur
      // Indels are ultimately reported as (position, size) tuples, where position is the left-most
      // start coordinate and sizes are +/- for insertions and deletions, respectively
      while (base_index >= 0 && seq_index >= 0){
	// Update the homopolymer tract length
	int homopolymer_len = std::min(MAX_HOMOP_LEN, std::max(haplotype->homopolymer_length(block_index, base_index),
							       haplotype->homopolymer_length(block_index, std::max(0, base_index-1))));

	if (matrix_type != prev_matrix_type){
	  // Record any processed indels
	  if (prev_matrix_type == DEL){
	    if (haplotype->reversed())
	      trace.add_flank_indel(std::pair<int32_t,int32_t>(indel_position, indel_position - pos));
	    else
	      trace.add_flank_indel(std::pair<int32_t,int32_t>(pos+1, pos - indel_position));
	  }
	  else if (prev_matrix_type == INS)
	    trace.add_flank_indel(std::pair<int32_t,int32_t>(indel_position + (haplotype->reversed() ? 0 : 1), indel_seq_index - seq_index));

	  // Initialize fields for indels that just began
	  if (matrix_type == DEL || matrix_type == INS){
	    indel_seq_index = seq_index;
	    indel_position  = pos;
	  }
	  prev_matrix_type = matrix_type;
        }

	// Extract alignment character for current base and update indices
	switch (matrix_type){
	case MATCH:
	  if (block_seq[base_index] != read_seq[seq_index] &&  base_log_correct[seq_index] > MIN_SNP_LOG_PROB_CORRECT)
	    trace.add_flank_snp(pos, read_seq[seq_index]);
	  flank_ss << read_seq[seq_index];
	  aln_ss << "M";
	  seq_index--;
	  base_index--;
	  pos += increment;
	  break;
	case DEL:
	  trace.inc_flank_del();
	  aln_ss << "D";
	  base_index--;
	  pos += increment;
	  break;
	case INS:
	  trace.inc_flank_ins();
	  flank_ss << read_seq[seq_index];
	  aln_ss << "I";
	  seq_index--;
	  break;
	default:
	  printErrorAndDie("Invalid matrix type when retracing alignments");
	  break;
	}

	if (seq_index == -1 || (base_index == -1 && block_index == 0)){
	  while(seq_index != -1){
	    aln_ss << "S";
	    seq_index--;
	  }

	  std::string flank_seq = flank_ss.str();
	  if (haplotype->reversed())
	    trace.add_flank_data(haplotype->num_blocks()-1-block_index, flank_seq);
	  else {
	    std::reverse(flank_seq.begin(), flank_seq.end());
	    trace.add_flank_data(block_index, flank_seq);
	  }
	  return aln_ss.str();
	}

	int best_opt;
	switch (matrix_type){
	case MATCH:
	  best_opt = triple_index_fn(insert_matrix[matrix_index-1]           + LOG_MATCH_TO_INS[homopolymer_len],
				     deletion_matrix[matrix_index-seq_len-1] + LOG_MATCH_TO_DEL[homopolymer_len],
				     match_matrix[matrix_index-seq_len-1]    + LOG_MATCH_TO_MATCH[homopolymer_len]);
	  if (best_opt == 0){
	    matrix_type   = INS;
	    matrix_index -= 1;
	  }
	  else if (best_opt == 1){
	    matrix_type   = DEL;
	    matrix_index -= (seq_len + 1);
	  }
	  else {
	    matrix_type   = MATCH;
	    matrix_index -= (seq_len + 1);
	  }
	  break;
	case DEL:
	  best_opt = pair_index_fn(deletion_matrix[matrix_index-seq_len] + LOG_DEL_TO_DEL,
				   match_matrix[matrix_index-seq_len]    + LOG_DEL_TO_MATCH);
	  if (best_opt == 0){
	    matrix_type   = DEL;
	    matrix_index -= seq_len;
	  }
	  else {
	    matrix_type   = MATCH;
	    matrix_index -= seq_len;
	  }
	  break;
	case INS:
	  best_opt = pair_index_fn(insert_matrix[matrix_index-1]        + LOG_INS_TO_INS,
				   match_matrix[matrix_index-seq_len-1] + LOG_INS_TO_MATCH);
	  if (best_opt == 0){
	    matrix_type   = INS;
	    matrix_index -= 1;
	  }
	  else {
	    matrix_type   = MATCH;
	    matrix_index -= (seq_len + 1);
	  }
	  break;
	default:
	  printErrorAndDie("Invalid matrix type when retracing alignments");
	  break;
	}
      }

      std::string flank_seq = flank_ss.str();
      if (haplotype->reversed())
	trace.add_flank_data(haplotype->num_blocks()-1-block_index, flank_seq);
      else {
	std::reverse(flank_seq.begin(), flank_seq.end());
	trace.add_flank_data(block_index, flank_seq);
      }
    }
    base_index = haplotype->get_seq(--block_index).size()-1;
  }
  return aln_ss.str();
}

void HapAligner::process_read(const Alignment& aln, int seed_base, const BaseQuality* base_quality, bool retrace_aln,
			      double* prob_ptr, AlignmentTrace& trace){
  //assert(seed_base != -1);
  assert(aln.get_sequence().size() == aln.get_base_qualities().size());
  // Extract probabilites related to base quality scores
  double* base_log_wrong   = new double[aln.get_sequence().size()]; // log10(Prob(error))
  double* base_log_correct = new double[aln.get_sequence().size()]; // log10(Prob(correct))
  const std::string& qual_string = aln.get_base_qualities();
  for (unsigned int j = 0; j < qual_string.size(); j++){
    base_log_wrong[j]   = base_quality->log_prob_error(qual_string[j]);
    base_log_correct[j] = base_quality->log_prob_correct(qual_string[j]);
//    base_log_wrong[j] = -10.5392;
//    base_log_correct[j] = -7.9436e-05;
}
  // Here we trim the read sequence because the flanking regions are too long for alignment to haplotypes
  std::string base_seq_str;
  trim_alignment(aln, base_seq_str);
  if (base_seq_str.size() == 0){
        base_seq_str += fw_haplotype_->get_first_block()->get_seq(0).substr(fw_haplotype_->get_first_block()->get_seq(0).size() - 5, 5);
        base_seq_str += fw_haplotype_->get_last_block()->get_seq(0).substr(0,5);
    }
  int base_seq_len = (int)base_seq_str.size();
  const char* base_seq = base_seq_str.c_str();
  seed_base = base_seq_len - 1;

  // Allocate scoring matrices based on the maximum haplotype size
  int max_hap_size          = fw_haplotype_->max_size();
  int num_hap_blocks        = fw_haplotype_->num_blocks();
  double max_LL             = -100000000;

  // Reverse bases and quality scores for the right flank
  std::string rev_rseq = base_seq_str.substr(seed_base+1);
  std::reverse(rev_rseq.begin(), rev_rseq.end());
  std::reverse(base_log_wrong+seed_base+1,   base_log_wrong+base_seq_len);
  std::reverse(base_log_correct+seed_base+1, base_log_correct+base_seq_len);

  // True iff we should reuse alignment information from the previous haplotype to accelerate computations
  bool reuse_alns = false;

  do {
    if (!realign_to_hap_[fw_haplotype_->cur_index()]){
      prob_ptr++;
      reuse_alns = false;
      continue;
    }
    // Perform alignment to current haplotype
    double l_prob, r_prob;
    int max_index;

    align_seq_to_hap(fw_haplotype_, reuse_alns, base_seq, seed_base, l_prob, base_log_wrong, base_log_correct);
    double LL = l_prob;
    *prob_ptr = LL;
    prob_ptr++;
    reuse_alns = true;

    if (LL > max_LL){
      max_LL = LL;
    }
  } while (fw_haplotype_->next() && rev_haplotype_->next());
  fw_haplotype_->reset();
  rev_haplotype_->reset();
}

AlignmentTrace* HapAligner::trace_optimal_aln(const Alignment& orig_aln, int seed_base, int best_haplotype, const BaseQuality* base_quality){
  fw_haplotype_->go_to(best_haplotype);
  fw_haplotype_->fix();
  rev_haplotype_->go_to(best_haplotype);
  fw_haplotype_->fix();
  double prob;
  AlignmentTrace* trace = new AlignmentTrace(fw_haplotype_->num_blocks());
  process_read(orig_aln, seed_base, base_quality, true, &prob, *trace);
  fw_haplotype_->unfix();
  rev_haplotype_->unfix();
  return trace;
}
