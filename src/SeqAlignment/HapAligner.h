#ifndef HAP_ALIGNER_H_
#define HAP_ALIGNER_H_

#include <assert.h>
#include <string>
#include <vector>

#include "AlignmentTraceback.h"
#include "../base_quality.h"
#include "Haplotype.h"

class AlignmentModel {
public:

    unsigned int MAX_HOMOP_LEN;
    float LOG_INS_TO_INS;
    float LOG_INS_TO_MATCH;
    float LOG_DEL_TO_DEL;
    float LOG_DEL_TO_MATCH;
    float LOG_MATCH_TO_MATCH;
    float LOG_MATCH_TO_INS;
    float LOG_MATCH_TO_DEL;

    // Constructor
    AlignmentModel(unsigned int max_homop_len, float log_ins_to_ins, float log_ins_to_match,
                   float log_del_to_del, float log_del_to_match, float log_match_to_match,
                   float log_match_to_ins, float log_match_to_del){
    MAX_HOMOP_LEN = max_homop_len;
    LOG_INS_TO_INS = log_ins_to_ins;
    LOG_INS_TO_MATCH = log_ins_to_match;
    LOG_DEL_TO_DEL = log_del_to_del;
    LOG_DEL_TO_MATCH = log_del_to_match;
    LOG_MATCH_TO_MATCH = log_match_to_match;
    LOG_MATCH_TO_INS = log_match_to_ins;
    LOG_MATCH_TO_DEL = log_match_to_del;
    }
};

class HapAligner {
 private:
  Haplotype* fw_haplotype_;
  Haplotype* rev_haplotype_;
  std::vector<bool> realign_to_hap_;
  std::vector<HapBlock*> rev_blocks_;
  std::vector<int32_t> repeat_starts_;
  std::vector<int32_t> repeat_ends_;
  int INDEL_FLANK_LEN;
  int SWITCH_OLD_ALIGN_LEN;
  AlignmentModel* AlnModel;


  void needleman_wunsch(const std::string& cent_seq, const std::string& read_seq, int& score) const;

  void align_seq_to_hap(Haplotype* haplotype, bool reuse_alns,
			const char* seq_0, int seq_len, double& left_prob);

   void align_seq_to_hap_short(Haplotype* haplotype, bool reuse_alns,
			const char* seq_0, int seq_len,
			const double* base_log_wrong, const double* base_log_correct,
			double* match_matrix, double* insert_matrix, double* deletion_matrix,
			int* best_artifact_size, int* best_artifact_pos, double& left_prob);

  /**
   * Compute the log-probability of the alignment given the alignment matrices for the left and right segments.
   * Stores the index of the haplotype position with which the seed base is aligned in the maximum likelihood alignment
   **/
  double compute_aln_logprob(int base_seq_len, int seed_base,
			     char seed_char, double log_seed_wrong, double log_seed_correct,
			     double* l_match_matrix, double* l_insert_matrix, double* l_deletion_matrix, double l_prob,
			     double* r_match_matrix, double* r_insert_matrix, double* r_deletion_matrix, double r_prob,
			     int& max_index);

  std::string retrace(Haplotype* haplotype, const char* read_seq, const double* base_log_correct,
		      int seq_len, int block_index, int base_index, int matrix_index, double* l_match_matrix,
		      double* l_insert_matrix, double* l_deletion_matrix, int* best_artifact_size, int* best_artifact_pos,
		      AlignmentTrace& trace);

  void calc_best_seed_position(int32_t region_start, int32_t region_end,
			       int32_t& best_dist, int32_t& best_pos);


  /**
  * Trim alignment before aligning to haplotype
  **/

  void trim_alignment(const Alignment& aln, std::string& trimmed_seq);


  // Private unimplemented copy constructor and assignment operator to prevent operations
  HapAligner(const HapAligner& other);
  HapAligner& operator=(const HapAligner& other);

 public:
  HapAligner(Haplotype* haplotype, std::vector<bool>& realign_to_haplotype, int INDEL_FLANK_LEN_, int SWITCH_OLD_ALIGN_LEN_,
                    std::vector<float>& alignment_model_params_){
    assert(realign_to_haplotype.size() == haplotype->num_combs());
    fw_haplotype_   = haplotype;
    rev_haplotype_  = haplotype->reverse(rev_blocks_);
    realign_to_hap_ = realign_to_haplotype;
    INDEL_FLANK_LEN = INDEL_FLANK_LEN_;
    SWITCH_OLD_ALIGN_LEN = SWITCH_OLD_ALIGN_LEN_;

    for (int i = 0; i < fw_haplotype_->num_blocks(); i++){
      HapBlock* block = fw_haplotype_->get_block(i);
      if (block->get_repeat_info() != NULL){
        repeat_starts_.push_back(block->start());
        repeat_ends_.push_back(block->end());
      }
    }

    if (!alignment_model_params_.empty()){
        AlnModel = new AlignmentModel(10, alignment_model_params_[0], alignment_model_params_[1],
                                          alignment_model_params_[2], alignment_model_params_[3],
                                          alignment_model_params_[4], alignment_model_params_[5],
                                          alignment_model_params_[6]);
    }
    else {
        AlnModel = new AlignmentModel(10, -1.0, -0.458675, -1.0, -0.458675, -0.00005800168, -10.448214728, -10.448214728); //default values from Dindel model, works best for Illumina reads and PacBio HiFi
    }
  }

  ~HapAligner(){
    for (unsigned int i = 0; i < rev_blocks_.size(); i++)
      delete rev_blocks_[i];
    rev_blocks_.clear();
    delete rev_haplotype_;
  }

  /** 
   * Returns the 0-based index into the sequence string that should be used as the seed for alignment or -1 if no valid seed exists
   **/
  int calc_seed_base(const Alignment& alignment);

  void process_read(const Alignment& aln, int seed_base, const BaseQuality* base_quality, bool retrace_aln,
		    double* prob_ptr, AlignmentTrace& traced_aln, int short_);

  void process_reads(const std::vector<Alignment>& alignments, int init_read_index, const BaseQuality* base_quality, const std::vector<bool>& realign_read,
		     double* aln_probs, int* seed_positions);

  /*
    Retraces the Alignment's optimal alignment to the provided haplotype.
    Returns the result as a new Alignment relative to the reference haplotype
   */
  AlignmentTrace* trace_optimal_aln(const Alignment& orig_aln, int seed_base, int best_haplotype, const BaseQuality* base_quality);
};

#endif
