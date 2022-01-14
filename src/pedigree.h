#ifndef PEDIGREE_H_
#define PEDIGREE_H_

#include <algorithm>
#include <assert.h>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "error.h"
#include "vcf_reader.h"

class NuclearFamily {
 private:
  std::string family_id_;
  std::string mother_, father_;
  std::vector<std::string> children_;
  std::vector<std::string> samples_;
  std::vector<int> vcf_indices_;

 public:
  NuclearFamily(const std::string& family_id, const std::string& mother, const std::string& father, const std::vector<std::string>& children)
    : family_id_(family_id), mother_(mother), father_(father), children_(children){
    samples_.push_back(mother_);
    samples_.push_back(father_);
    samples_.insert(samples_.end(), children_.begin(), children_.end());
  }

  void load_vcf_indices(const VCF::VCFReader& vcf_reader){
    vcf_indices_.clear();
    for (int i = 0; i < samples_.size(); i++){
      if (!vcf_reader.has_sample(samples_[i]))
	printErrorAndDie("No sample data available in VCF");
      vcf_indices_.push_back(vcf_reader.get_sample_index(samples_[i]));
    }
  }

  void clear_vcf_indices(){
    vcf_indices_.clear();
  }

  const std::string& get_family_id() const { return family_id_; }
  const std::string& get_mother()    const { return mother_; }
  const std::string& get_father()    const { return father_; }
  const std::vector<std::string>& get_children() const { return children_; }
  const std::vector<std::string>& get_samples()  const { return samples_; }
  int size()         const { return 2 + children_.size(); }
  int num_children() const { return children_.size();     }

  bool is_missing_sample(const std::set<std::string>& samples) const {
    for (auto sample_iter = samples_.begin(); sample_iter != samples_.end(); sample_iter++)
      if (samples.find(*sample_iter) == samples.end())
        return true;
    return false;
  }

  bool is_missing_genotype(const VCF::Variant& variant) const {
    if (vcf_indices_.empty())
      printErrorAndDie("No VCF indices were preloaded in the NuclearFamily");
    for (auto index_iter = vcf_indices_.begin(); index_iter != vcf_indices_.end(); index_iter++)
      if (variant.sample_call_missing(*index_iter))
	return true;
    return false;
  }

  bool is_mendelian(const VCF::Variant& variant) const {
    if (vcf_indices_.empty())
      printErrorAndDie("No VCF indices were preloaded in the NuclearFamily");

    int m_1, m_2, f_1, f_2, c_1, c_2;
    variant.get_genotype(vcf_indices_[0], m_1, m_2);
    variant.get_genotype(vcf_indices_[1], f_1, f_2);

    for (int i = 2; i < vcf_indices_.size(); i++){
      variant.get_genotype(vcf_indices_[i], c_1, c_2);
      if ((c_1 != m_1 && c_1 != m_2) || (c_2 != f_1 && c_2 != f_2))
	if ((c_1 != f_1 && c_1 != f_2) || (c_2 != m_1 && c_2 != m_2))
	  return false;
    }
    return true;
  }
};

class PedigreeNode {
 private:
  std::string name_;
  PedigreeNode* mother_;
  PedigreeNode* father_;
  std::vector<PedigreeNode*> children_;
  std::string family_id_;

 public:
  PedigreeNode(const std::string& name, const std::string& family_id)
    : name_(name), family_id_(family_id){
    mother_ = NULL;
    father_ = NULL;
  }

  ~PedigreeNode(){
    children_.clear();
  }

  bool has_mother()                 const { return mother_ != NULL;  }
  bool has_father()                 const { return father_ != NULL;  }
  PedigreeNode* get_mother()        const { return mother_;          }
  PedigreeNode* get_father()        const { return father_;          }
  const std::string& get_name()     const { return name_;            }
  const std::string& get_family()   const { return family_id_;       }
  std::vector<PedigreeNode*>& get_children() { return children_;     }

  void set_mother(PedigreeNode* mother) { mother_ = mother;           }
  void set_father(PedigreeNode* father) { father_ = father;           }
  void add_child (PedigreeNode* child)  { children_.push_back(child); }
  void del_child (PedigreeNode* child)  {
    auto iter = std::find(children_.begin(), children_.end(), child);;
    if (iter == children_.end())
      printErrorAndDie("Can't delete child from node as it is not among the existing children");
    children_.erase(iter);
  }
  
  void print(std::ostream& out) const {
    out << "NAME:"     << name_
	<< "\tFATHER:" << (father_ == NULL ? "NONE" : father_->get_name())
	<< "\tMOTHER:" << (mother_ == NULL ? "NONE" : mother_->get_name()) << std::endl; 
  }
};

class PedigreeGraph {
 private:
  // Nodes that don't have any ancestors 
  std::vector<PedigreeNode*> no_ancestors_;
  
  // Nodes that don't have any descendants
  std::vector<PedigreeNode*> no_descendants_;

  // Nodes in graph sorted in topological order
  std::vector<PedigreeNode*> nodes_;

  bool topological_sort(std::vector<PedigreeNode*>& nodes);
  bool build(const std::string& input_file);
  void init_no_ancestors();
  void init_no_descendants();
  bool build_subgraph(std::vector<PedigreeNode*>& sorted_nodes);

  // Private unimplemented copy constructor and assignment operator to prevent operations
  PedigreeGraph(const PedigreeGraph& other);
  PedigreeGraph& operator=(const PedigreeGraph& other);
    
 public:
  PedigreeGraph(){}

  explicit PedigreeGraph(const std::string& input_file){
    bool success = build(input_file);
    if (!success)
      printErrorAndDie("Supplied pedigree file " + input_file + " contains cycles");
    init_no_ancestors();
    init_no_descendants();
  }

  explicit PedigreeGraph(std::vector<PedigreeNode*>& subgraph_nodes){
    if (!build_subgraph(subgraph_nodes))
      printErrorAndDie("Subgraph in pedigree contains a cycle");
    init_no_ancestors();
    init_no_descendants();
  }

  ~PedigreeGraph(){
    for (int i = 0; i < nodes_.size(); i++)
      delete nodes_[i];
    nodes_.clear();
    no_ancestors_.clear();
    no_descendants_.clear();
  }

  int size() const { return nodes_.size(); }

  void print(std::ostream& out) const {
    out << "Pedigree graph contains " << nodes_.size() << " nodes" << std::endl;
  }

  void prune(const std::set<std::string>& sample_set);

  void split_into_connected_components(std::vector<PedigreeGraph*>& components);

  bool is_nuclear_family() const;

  NuclearFamily convert_to_nuclear_family() const;
};

void extract_pedigree_nuclear_families(const std::string& pedigree_fam_file, const std::set<std::string>& samples_with_data,
                                       std::vector<NuclearFamily>& nuclear_families, std::ostream& logger);

#endif
