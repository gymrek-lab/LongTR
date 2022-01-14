#ifndef ALIGNMENT_DATA_H_
#define ALIGNMENT_DATA_H_

#include <assert.h>
#include <sstream>
#include <string>
#include <vector>

#include "../base_quality.h"
#include "../error.h"

class CigarElement {
 private:
  char type_;
  int  num_;

 public:
  CigarElement(char type, int num){
    type_ = type;
    num_  = num;
  }

  inline void set_type(char type){ type_ = type;}
  inline void set_num(int num)   {  num_ = num; }
  inline char get_type()   const { return type_; }
  inline int  get_num()    const { return num_;  }
};

class Alignment {
 private:
  int32_t start_;
  int32_t stop_;
  std::vector<CigarElement> cigar_list_;
  std::string name_;
  std::string base_qualities_;
  std::string sequence_;
  std::string alignment_;
  std::vector<bool> use_for_haps_;
  bool rev_strand_;

 public:
   Alignment(int32_t start, int32_t stop, bool rev_strand,
	    const std::string& name,
	    const std::string& base_qualities,
	    const std::string& sequence,
	    const std::string& alignment)
    : name_(name), base_qualities_(base_qualities), sequence_(sequence), alignment_(alignment){
    start_      = start;
    stop_       = stop;
    rev_strand_ = rev_strand;
  }

  explicit Alignment(const std::string& name)
    : name_(name){
    start_      = 0;
    stop_       = -1;
    rev_strand_ = false;
  }

  inline const std::string& get_name()   const { return name_;   }
  inline int32_t get_start()             const { return start_;  }
  inline int32_t get_stop()              const { return stop_;   }
  inline void set_start(int32_t start)         { start_ = start; }
  inline void set_stop(int32_t stop)           { stop_  = stop;  }

  bool operator<(const Alignment &aln) const {
    if (start_ != aln.get_start())
      return start_ < aln.get_start();
    if (stop_ != aln.get_stop())
      return stop_ < aln.get_stop();
    return false;
  }

  void check_CIGAR_string(){
    unsigned int num = 0;
    for (auto iter = cigar_list_.begin(); iter != cigar_list_.end(); ++iter)
      if (iter->get_type() != 'D' && iter->get_type() != 'H')
	num += iter->get_num();
    if (num != sequence_.size()){
      std::cerr << "CIGAR check failed for read " << name_ << ": "
		<< num << " " << sequence_.size() << std::endl
		<< sequence_  << std::endl
		<< alignment_ << std::endl
		<< getCigarString() << std::endl;
      assert(false);
    }
  }

  int num_indels() const {
    int num = 0;
    for (auto iter = cigar_list_.begin(); iter != cigar_list_.end(); ++iter)
      if (iter->get_type() == 'I' || iter->get_type() == 'D')
	num++;
    return num;
  }
  
  int num_mismatches() const {
    int num = 0;
    for (auto iter = cigar_list_.begin(); iter != cigar_list_.end(); ++iter)
      if (iter->get_type() == 'X')
	num++;
    return num;
  }

  int num_matched_bases() const {
    int num = 0;
    for (auto iter = cigar_list_.begin(); iter != cigar_list_.end(); ++iter)
      if (iter->get_type() == 'M' || iter->get_type() == '=')
	num += iter->get_num();
    return num;
  }

  inline void set_base_qualities(const std::string& base_qualities)       { base_qualities_.assign(base_qualities); }
  inline void set_sequence(const std::string& sequence)                   { sequence_.assign(sequence);             }
  inline void set_alignment(const std::string& alignment)                 { alignment_.assign(alignment);           }
  inline void set_hap_gen_info(const std::vector<bool>& use_for_haps)     { use_for_haps_ = use_for_haps;           }
  inline void add_cigar_element(CigarElement e)                           { cigar_list_.push_back(e);               }
  inline void set_cigar_list(const std::vector<CigarElement>& cigar_list) {
    cigar_list_.clear();
    for (unsigned int i = 0; i < cigar_list.size(); i++)
      cigar_list_.push_back(cigar_list[i]);
  }

  inline const std::string& get_base_qualities()           const { return base_qualities_; }
  inline const std::string& get_sequence()                 const { return sequence_;       }
  inline const std::string& get_alignment()                const { return alignment_;      }
  inline const std::vector<CigarElement>& get_cigar_list() const { return cigar_list_;     }
  bool use_for_hap_generation(int region_index) const { return use_for_haps_[region_index]; }
  bool is_from_reverse_strand() const { return rev_strand_; }

  std::string getCigarString() const {
    std::stringstream cigar_str;
    for (auto iter = cigar_list_.begin(); iter != cigar_list_.end(); iter++)
      cigar_str << iter->get_num() << iter->get_type();
    return cigar_str.str();
  }
};

#endif
