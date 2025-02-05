#ifndef REGION_H_
#define REGION_H_

#include <assert.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <cstdint>

#include "error.h"


class Region {
private:
    std::string chrom_, name_, motifs_;
    int32_t start_, stop_;
    int period_;

    // Helper function to compute the period

    std::vector<std::string> splitMotifs(const std::string& motifs, char delimiter = ',') const {
        std::vector<std::string> result;
        std::stringstream ss(motifs);
        std::string item;

        while (std::getline(ss, item, delimiter)) {
            result.push_back(item);
        }

        return result;
    }
    int computePeriod(const std::string& motifs) {
        std::vector<std::string> motif_list = splitMotifs(motifs);
        std::set<int> period_list;
        for (const auto& motif : motif_list) {
            period_list.insert(motif.size());
        }
        return (period_list.size() == 1) ? *period_list.begin() : -1;
    }




public:
    // Constructors
    Region(const std::string& chrom, int32_t start, int32_t stop, const std::string& motifs)
        : chrom_(chrom), start_(start), stop_(stop), motifs_(motifs), period_(computePeriod(motifs)) {
        assert(stop > start);
    }

    Region(const std::string& chrom, int32_t start, int32_t stop, const std::string& motifs, const std::string& name)
        : chrom_(chrom), name_(name), start_(start), stop_(stop), motifs_(motifs), period_(computePeriod(motifs)) {
        assert(stop > start);
    }

    // Getters
    const std::string& chrom() const { return chrom_; }
    const std::string& name()  const { return name_; }
    int32_t start()           const { return start_; }
    int32_t stop()            const { return stop_; }
    int period()              const { return period_; }
    const std::string& motif() const { return motifs_; }
    std::string period_str() const {
        std::vector<std::string> motif_list = splitMotifs(motif());
        std::ostringstream oss;
        for (int i = 0; i < motif_list.size(); ++i) {
            if (i > 0) oss << ",";
            oss << motif_list[i].size();
        }
        return oss.str();
    }


    // Clone function
    Region* copy() const { return new Region(chrom_, start_, stop_, motifs_, name_); }

    // Setters
    void set_start(int32_t start) { start_ = start; }
    void set_stop(int32_t stop) { stop_ = stop; }

    // String representation
    std::string str() const {
        std::stringstream ss;
        ss << chrom_ << ":" << start_ << "-" << stop_;
        return ss.str();
    }

    // Comparison operator for sorting
    bool operator<(const Region& r) const {
        if (chrom_ != r.chrom()) return chrom_ < r.chrom();
        if (start_ != r.start()) return start_ < r.start();
        return stop_ < r.stop();
    }
};

void readRegions(const std::string& input_file, uint32_t max_regions, const std::string& chrom_limit, std::vector<Region>& regions, std::ostream& logger);

void orderRegions(std::vector<Region>& regions);


bool isValidMotif(const std::string& motif);

class RegionGroup {
  std::vector<Region> regions_;
  std::string chrom_;
  int32_t start_;
  int32_t stop_;

 public:
  explicit RegionGroup(const Region& region){
    regions_.push_back(region);
    chrom_ = region.chrom();
    start_ = region.start();
    stop_  = region.stop();
  }

  const std::vector<Region>& regions() const {
    return regions_;
  }

  const std::string& chrom() const { return chrom_;          }
  int32_t  start()           const { return start_;          }
  int32_t  stop()            const { return stop_;           }
  int      num_regions()     const { return regions_.size(); }

  RegionGroup* copy() const {
    RegionGroup* clone = new RegionGroup(regions_[0]);
    for (int i = 1; i < regions_.size(); i++)
      clone->add_region(regions_[i]);
    return clone;
  }

  void add_region(const Region& region){
    if (region.chrom().compare(chrom_) != 0)
      printErrorAndDie("RegionGroup can only consist of regions on a single chromosome");
    start_ = std::min(start_, region.start());
    stop_  = std::max(stop_,  region.stop());
    regions_.push_back(region);
    std::sort(regions_.begin(), regions_.end());
  }
};

#endif
