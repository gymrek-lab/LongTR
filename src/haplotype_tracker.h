#ifndef HAPLOTYPE_TRACKER_H_
#define HAPLOTYPE_TRACKER_H_

#include <climits>
#include <deque>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "vcf_reader.h"
#include "pedigree.h"

class DiploidEditDistance {
 private:
  int distances_[4];

 public:
  DiploidEditDistance(int d11, int d12, int d21, int d22){
    distances_[0] = d11;
    distances_[1] = d12;
    distances_[2] = d21;
    distances_[3] = d22;
  }

  int distance(int index_a, int index_b) const {
    if (index_a < 0 || index_a > 1 || index_b < 0 || index_b > 1)
      printErrorAndDie("Index for distance() function in DiplodEditDistance class must be 0 or 1");
    return distances_[index_a*2 + index_b];
  }

  void min_distance(int& dist, int& index) const {
    dist  = distances_[0];
    index = 0;
    for (int i = 1; i < 4; i++)
      if (distances_[i] < dist){
	dist  = distances_[i];
	index = i;
      }
  }

  void second_min_distance(int& second_dist, int& second_index) const {
    int dist    = INT_MAX, index = -1;
    second_dist = INT_MAX;

    for (int i = 0; i < 4; i++){
      if (distances_[i] < dist){
	second_dist  = dist;
	second_index = index;
	dist  = distances_[i];
	index = i;
      }
      else if (distances_[i] < second_dist){
	second_dist  = distances_[i];
	second_index = i;
      }
    }
  }
};

class DiploidHaplotype {
 private:
  std::deque<int64_t> snps_1_, snps_2_;
  int64_t erase_mask_, set_mask_;
  int num_removed_;
  const static int64_t MAX_DIGIT = (1ULL << 63);

  int num_set_bits(int64_t value) const {
    int count = 0;
    while (value != 0){
      value &= (value-1);
      count++;
    }
    return count;
  }

  int edit_distance(const std::deque<int64_t>& snp_hap_1, const std::deque<int64_t>& snp_hap_2) const {
    assert(snp_hap_1.size() == snp_hap_2.size());
    int distance = 0;
    auto iter_a = snp_hap_1.begin();
    auto iter_b = snp_hap_2.begin();
    while (iter_a != snp_hap_1.end()){
      distance += num_set_bits((*iter_a) ^ (*iter_b));
      iter_a++; iter_b++;
    }
    return distance;
  }

  void extract_set_bits(int64_t value, int offset, std::set<int>& mismatch_indices) const;

 public:  
  DiploidHaplotype(){
    reset();
  }

  DiploidEditDistance edit_distances(const DiploidHaplotype& other_hap) const {
    int d11 = edit_distance(snps_1_, other_hap.snps_1_);
    int d12 = edit_distance(snps_1_, other_hap.snps_2_);
    int d21 = edit_distance(snps_2_, other_hap.snps_1_);
    int d22 = edit_distance(snps_2_, other_hap.snps_2_);
    return DiploidEditDistance(d11, d12, d21, d22);
  }

  void add_mismatched_sites(int hap_index, DiploidHaplotype& other_hap, int other_index,
			    std::set<int>& mismatch_indices) const;

  void add_snp(int gt_a, int gt_b);

  void remove_next_snp();

  void reset(){
    snps_1_.clear();
    snps_2_.clear();
    snps_1_.push_back(0);
    snps_2_.push_back(0);
    erase_mask_  = -2;
    set_mask_    = 1;
    num_removed_ = 0;
  }
};

class HaplotypeTracker {
 private:
  std::string chrom_;
  std::vector<NuclearFamily> families_;
  std::vector<std::string> samples_;
  std::vector<int> vcf_indices_;
  std::map<std::string, int> sample_indices_;
  std::vector<DiploidHaplotype> snp_haplotypes_;
  VCF::VCFReader snp_vcf_;
  int32_t window_size_;
  int32_t num_snps_;
  std::deque<int32_t> positions_;
  int32_t prev_window_start_, prev_window_end_;

  // Private unimplemented copy constructor and assignment operator to prevent operations
  HaplotypeTracker(const HaplotypeTracker& other);
  HaplotypeTracker& operator=(const HaplotypeTracker& other);

  int32_t next_snp_position() const {
    if (num_snps_ == 0)
      return -1;
    return positions_.front();
  }

  int32_t last_snp_position() const {
    if (num_snps_ == 0)
      return -1;
    return positions_.back();
  }

  void remove_next_snp(){
    if (num_snps_ == 0)
      return;
    num_snps_--;
    for (unsigned int i = 0; i < snp_haplotypes_.size(); i++)
      snp_haplotypes_[i].remove_next_snp();
    positions_.pop_front();
  }

  void reset(){
    num_snps_  = 0;
    positions_ =  std::deque<int32_t>();
    prev_window_start_ = -1;
    prev_window_end_   = -1;
    for (unsigned int i = 0; i < snp_haplotypes_.size(); i++)
      snp_haplotypes_[i].reset();
  }

  void add_snp(const VCF::Variant& variant);

 public:
 HaplotypeTracker(const std::vector<NuclearFamily>& families, const std::string& snp_vcf_file, int32_t window_size)
   : families_(families), snp_vcf_(snp_vcf_file){
    window_size_ = window_size;
    for (auto family_iter = families_.begin(); family_iter != families_.end(); family_iter++){
      family_iter->load_vcf_indices(snp_vcf_);
      samples_.insert(samples_.end(),  family_iter->get_samples().begin(),  family_iter->get_samples().end());
    }
    for (unsigned int i = 0; i < samples_.size(); i++){
      if (!snp_vcf_.has_sample(samples_[i]))
	printErrorAndDie("No sample data available in VCF");
      vcf_indices_.push_back(snp_vcf_.get_sample_index(samples_[i]));
      sample_indices_[samples_[i]] = i;
    }

    snp_haplotypes_    = std::vector<DiploidHaplotype>(samples_.size(), DiploidHaplotype());
    num_snps_          = 0;
    prev_window_start_ = -1;
    prev_window_end_   = -1;
  }

  const std::vector<NuclearFamily>& families() const {
    return families_;
  }

  int32_t num_stored_snps() const { return num_snps_; }

  DiploidEditDistance edit_distances(const std::string& sample_1, const std::string& sample_2) const {
    int index_1 = sample_indices_.find(sample_1)->second;
    int index_2 = sample_indices_.find(sample_2)->second;
    return snp_haplotypes_[index_1].edit_distances(snp_haplotypes_[index_2]);
  }

  void advance(const std::string& chrom, int32_t pos, const std::set<std::string>& sites_to_skip);

  bool infer_haplotype_inheritance(const NuclearFamily& family, int max_best_score, int min_second_best_score,
				   std::vector<int>& maternal_indices, std::vector<int>& paternal_indices, std::set<int32_t>& bad_sites);
};

#endif
