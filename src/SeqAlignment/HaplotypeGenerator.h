#ifndef HAPLOTYPE_GENERATOR_H_
#define HAPLOTYPE_GENERATOR_H_

#include <chrono>
#include <cstdint>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "AlignmentData.h"
#include "../error.h"
#include "../region.h"
#include "../stutter_model.h"
#include "Haplotype.h"
#include "HapBlock.h"

// Default seed for the POA cluster-subsampling RNG (HaplotypeGenerator::poa()). A fixed default makes
// runs bit-reproducible; override with --seed to probe how sensitive calls are to the subsample.
constexpr uint32_t DEFAULT_POA_SEED = 42;
// Cluster size at or above which poa() subsamples instead of using every sequence. Trades POA cost
// against consensus quality; also lets tests force the subsampling path on small fixtures.
constexpr int DEFAULT_POA_CLUSTER_LIMIT = 30;

class HaplotypeGenerator {
 private:
  // Criteria used to determine whether a candidate sequence should be identified as an allele
  double MIN_FRAC_READS;
  double MIN_FRAC_SAMPLES;
  double MIN_FRAC_STRONG_SAMPLE;
  double MIN_READS_STRONG_SAMPLE;
  double MIN_STRONG_SAMPLES;

  // When extracting alleles in regions, we pad by these amounts to improve the capture of proximal indels
  int LEFT_PAD, RIGHT_PAD;

  int32_t MIN_BLOCK_SPACING; // Minimum distance (bp) between variant haplotype blocks
  int32_t REF_FLANK_LEN;     // Maximum length of reference sequences flanking the variant haplotype blocks

  bool finished_; // True iff the underlying haplotype blocks are ready for downstream use
  std::string failure_msg_;
  int32_t min_aln_start_, max_aln_stop_;
  std::vector<HapBlock*> hap_blocks_;
  bool blocks_released_ = false; // True once get_haplotype_blocks() hands ownership to the caller (genotyper)

  // Wall-clock watchdog: poll during POA/clustering; aborts once past the per-locus deadline.
  // wd_aborted_ is mutable so wd_expired() can flip it from the const POA methods (gen_candidate_seqs, poa, ...).
  bool wd_enabled_ = false;
  mutable bool wd_aborted_ = false;
  std::chrono::steady_clock::time_point wd_deadline_;
  bool wd_expired() const {
    if (wd_enabled_ && std::chrono::steady_clock::now() > wd_deadline_){ wd_aborted_ = true; return true; }
    return false;
  }

  // Seed for the POA cluster subsampling in poa(), which previously used std::random_device and so
  // reseeded from system entropy on every call -- making genotypes differ between identical runs at
  // any locus with a cluster of >= cluster_size_limit sequences.
  //
  // poa() constructs a generator from this seed on EVERY call rather than drawing from one long-lived
  // stream. That matters: a per-locus stream would make each cluster's subsample depend on how many
  // clusters were drawn before it, so any change in cluster ordering (e.g. the read arrival order
  // differing between N single-sample BAMs and one multi-sample BAM) would change the calls. Reseeding
  // per call makes a cluster's subsample a function of the seed and that cluster alone.
  uint32_t poa_seed_ = DEFAULT_POA_SEED;
  int poa_cluster_limit_ = DEFAULT_POA_CLUSTER_LIMIT;

  void trim(int ideal_min_length,
	    int32_t& region_start, int32_t& region_end, std::vector<std::pair<std::string, bool>>& sequences) const;

  bool extract_sequence(const Alignment& aln, int32_t start, int32_t end, std::string& seq) const;

  void poa(const std::vector<std::string>& seqs, std::string& consensus) const;

  void needleman_wunsch(const std::string& cent_seq, const std::string& read_seq, int& score, int T) const;

  bool greedy_clustering(const std::vector<std::string>& seqs, std::map<std::string, std::vector<std::string>>& clusters, int t) const;

  bool merge_clusters(const std::vector<std::string>& new_centeroids, std::map<std::string, std::vector<std::string>>& clusters, int t) const;

  void gen_candidate_seqs(const std::string& ref_seq, int ideal_min_length,
			  const std::vector< std::vector<Alignment> >& alignments, const std::vector<std::string>& vcf_alleles,
			  int32_t& region_start, int32_t& region_end, std::vector<std::pair<std::string, bool>>& sequences,
			  bool log_alt, std::map<std::string, std::string>* admission_out) const;

  void get_aln_bounds(const std::vector< std::vector<Alignment> >& alignments,
		      int32_t& min_aln_start, int32_t& max_aln_stop) const;

  // Private unimplemented copy constructor and assignment operator to prevent operations
  HaplotypeGenerator(const HaplotypeGenerator& other);
  HaplotypeGenerator& operator=(const HaplotypeGenerator& other);

 public:
  HaplotypeGenerator(int32_t min_aln_start, int32_t max_aln_stop, int INDEL_FLANK_LEN_){
    finished_                = false;
    MIN_FRAC_READS           = 0.05;
    MIN_FRAC_SAMPLES         = 0.05;
    MIN_FRAC_STRONG_SAMPLE   = 0.2;
    MIN_READS_STRONG_SAMPLE  = 2;
    MIN_STRONG_SAMPLES       = 1;
    LEFT_PAD                 = INDEL_FLANK_LEN_;
    RIGHT_PAD                = INDEL_FLANK_LEN_;
    MIN_BLOCK_SPACING        = 10;
    REF_FLANK_LEN            = 35;
    min_aln_start_           = min_aln_start;
    max_aln_stop_            = max_aln_stop;
  }

  ~HaplotypeGenerator(){
    // Free any haplotype blocks we still own. On the success path they are handed to the genotyper via
    // get_haplotype_blocks() (which sets blocks_released_) and freed by its destructor; on any failure
    // path ownership was never transferred, so we free them here to avoid leaking on build failure.
    if (!blocks_released_)
      for (unsigned int i = 0; i < hap_blocks_.size(); i++)
        delete hap_blocks_[i];
  }

  bool add_vcf_haplotype_block(int32_t pos, const std::string& chrom_seq,
			       const std::vector<std::string>& vcf_alleles, const StutterModel* stutter_model);

  bool add_haplotype_block(const Region& region, const std::string& chrom_seq, const std::vector< std::vector<Alignment> >& alignments,
			   const std::vector<std::string>& vcf_alleles, const StutterModel* stutter_model,
			   bool log_alt, std::map<std::string, std::string>* admission_out);

  bool fuse_haplotype_blocks(const std::string& chrom_seq);

  // Seed for the POA cluster-subsampling RNG (see poa()).
  void set_seed(uint32_t seed){ poa_seed_ = seed; }

  // Cluster size at/above which poa() subsamples. Lower is cheaper but coarser; also lets small
  // test fixtures reach the subsampling path that otherwise only triggers at realistic depth.
  void set_poa_cluster_limit(int limit){ poa_cluster_limit_ = limit; }

  // Wall-clock watchdog hooks. set_deadline() shares the per-locus deadline; timed_out() reports abort.
  void set_deadline(std::chrono::steady_clock::time_point deadline, bool enabled){ wd_deadline_ = deadline; wd_enabled_ = enabled; }
  bool timed_out() const { return wd_aborted_; }

  const std::string& failure_msg(){ return failure_msg_; }

  const std::vector<HapBlock*> get_haplotype_blocks() {
    if (!finished_)
      printErrorAndDie("Haplotype blocks are not ready for downstream use");
    blocks_released_ = true; // ownership transferred to the caller; our destructor must not free these
    return hap_blocks_;
  }
};

#endif
