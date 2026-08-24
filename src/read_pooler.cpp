#include "read_pooler.h"
#include <string>

int32_t ReadPooler::add_alignment(Alignment& aln){
  if (pooled_)
    printErrorAndDie("Cannot call add_alignment function once pool() function has been invoked");
  
  // Pool key. Two reads may share a pool only if they pose the SAME alignment problem to HapAligner,
  // which reads the sequence, the coordinates AND the CIGAR off this Alignment: process_read()
  // (HapAligner.cpp:840) calls trim_alignment(), which walks the CIGAR to decide how many bases to clip
  // from each end before the DP, and calc_seed_base() walks it to choose the seed base.
  //
  // Keying on the sequence alone pooled reads that trim to the same window sequence (left_align_reads()
  // clips every read to locus +/- FLANK_SIZE) but carry different indel placements from their reference
  // alignments. The pool kept whichever read arrived first and copied its likelihoods to the rest, so a
  // legitimate alignment was silently discarded and the result depended on read -- hence sample --
  // ordering. Measured: 0.10% of pools, and at chr21:14161569 a 2273 log-unit shift in one sample's
  // alignment LL, a flipped MAP genotype and a different ALT set.
  //
  // Including the alignment in the key restores the premise pooling relies on: members are now
  // genuinely identical inputs, so sharing one DP result is exact and the output matches what you would
  // get with pooling disabled.
  //
  // Cost, measured on 10 ONT samples: +0.02% pools on an unbiased 492-locus sample (91,096 -> 91,114;
  // 485 of 492 loci unchanged, worst locus +2.8%). Pooling still removes 45.5% of alignments, unchanged
  // to one decimal. The tail is real but narrow: at loci selected *because* they had cross-sample
  // alignment conflicts the growth reaches +2.3%, and the worst single locus observed was +9.9%.
  //
  // '|' separates the fields so the concatenation is unambiguous; it cannot occur in a sequence (ACGTN)
  // or in a CIGAR string (digits and MIDNSHP=X).
  std::string pool_key = aln.get_sequence() + '|'
                       + std::to_string(aln.get_start()) + '|'
                       + std::to_string(aln.get_stop())  + '|'
                       + aln.getCigarString();

  auto pool_iter = seq_to_pool_.find(pool_key);
  if (pool_iter == seq_to_pool_.end()){
    seq_to_pool_[pool_key] = pool_index_;
    pooled_alns_.push_back(Alignment(aln.get_start(), aln.get_stop(), false, aln.get_deleted(), "READPOOL", "", aln.get_sequence(), aln.get_alignment()));
    pooled_alns_.back().set_cigar_list(aln.get_cigar_list());
    qualities_by_pool_.push_back(std::vector<const std::string*>());
    qualities_by_pool_.back().push_back(new std::string(aln.get_base_qualities()));
    return pool_index_++;
  }
  else{
    qualities_by_pool_[pool_iter->second].push_back(new std::string(aln.get_base_qualities()));
    return pool_iter->second;
  }  
}
