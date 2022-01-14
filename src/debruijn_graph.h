#ifndef DEBRUIJN_GRAPH_
#define DEBRUIJN_GRAPH_

#include <assert.h>
#include <algorithm>
#include <string>
#include <vector>

#include "directed_graph.h"

class DebruijnPath;

class DebruijnGraph : public DirectedGraph {
 protected:
  int k_;
  std::string ref_seq_;
  std::string source_kmer_;
  std::string sink_kmer_;
  int32_t num_strings_;
  std::vector<bool> ref_edge_; // True iff the edge at the corresponding index is from the reference sequence

  void get_alt_kmer_nodes(std::string& kmer, bool source, bool sink, std::vector<Node*>& nodes);

  void prune_edges(std::vector<bool>& remove_edges);

 public:
 DebruijnGraph(int k, const std::string& ref_seq) : ref_seq_(ref_seq){
    assert(ref_seq.size() > k);
    k_           = k;
    source_kmer_ = ref_seq.substr(0, k_);
    sink_kmer_   = ref_seq.substr(ref_seq.size()-k, k_);
    num_strings_ = 0;

    // Add the reference path with a weight of 2
    add_string(ref_seq, 2);
    ref_edge_.clear();
    ref_edge_.resize(edges_.size(), true);
  }

  void add_string(const std::string& seq, int weight=1);

  void enumerate_paths(int min_weight, int max_paths, std::vector<std::pair<std::string, int> >& paths);

  static bool calc_kmer_length(const std::string& ref_seq, int min_kmer, int max_kmer, int& kmer);

  bool is_source_ok();

  bool is_sink_ok();

  void prune_edges(double min_edge_freq, int min_weight);
};

class DebruijnPath {
 private:
  DebruijnPath* parent_;
  int node_id_;
  int min_weight_, max_weight_;
  int depth_;

 public:
  explicit DebruijnPath(int node_id){
    parent_     = NULL;
    min_weight_ = 1000000;
    max_weight_ = 0;
    node_id_    = node_id;
    depth_      = 0;
  }

  int get_min_weight()       const { return min_weight_; }
  int get_max_weight()       const { return max_weight_; }
  int get_node_id()          const { return node_id_;    }
  int get_depth()            const { return depth_;      }
  DebruijnPath* get_parent() const { return parent_;     }

  DebruijnPath* add_edge(Edge* edge){
    DebruijnPath* new_path = new DebruijnPath(edge->get_destination());
    new_path->parent_      = this;
    new_path->min_weight_  = std::min(min_weight_, edge->get_weight());
    new_path->max_weight_  = std::max(max_weight_, edge->get_weight());
    new_path->depth_       = depth_ + 1;
    return new_path;
  }
  
  std::string get_sequence(DebruijnGraph* graph) const;
};

bool path_comparator(DebruijnPath* p1, DebruijnPath* p2);

#endif
