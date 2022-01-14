#include <assert.h>

#include <set>

#include "denovos/denovo_scanner.h"
#include "error.h"
#include "snp_tree.h"

bool in_any_region(const VCF::Variant& variant, const std::vector<Region>& skip_regions, int32_t skip_padding){
  for (auto region_iter = skip_regions.begin(); region_iter != skip_regions.end(); region_iter++)
    if (variant.get_position() >= region_iter->start() - skip_padding)
      if (variant.get_position() <= region_iter->stop() + skip_padding)
	return true;
  return false;
}

void filter_snps(std::vector<SNP>& snps, const std::set<int32_t>& bad_sites){
  int insert_index = 0;
  for (int i = 0; i < snps.size(); i++)
    if (bad_sites.find(snps[i].pos()+1) == bad_sites.end()) // +1 required b/c bad sites are 1-based, while SNPs are 0-based
      snps[insert_index++] = snps[i];
  snps.resize(insert_index);
}

bool create_snp_trees(const std::string& chrom, uint32_t start, uint32_t end, const std::vector<Region>& skip_regions, int32_t skip_padding, VCF::VCFReader* snp_vcf, HaplotypeTracker* tracker,
                      std::map<std::string, unsigned int>& sample_indices, std::vector<SNPTree*>& snp_trees, std::ostream& logger){
  logger << "Building SNP tree for region " << chrom << ":" << start << "-" << end << std::endl;
  assert(sample_indices.size() == 0 && snp_trees.size() == 0);

  if (!snp_vcf->set_region(chrom, start, end))
    return false;

  // Index samples
  unsigned int sample_count = 0;
  const std::vector<std::string>& vcf_samples = snp_vcf->get_samples();
  for (auto sample_iter = vcf_samples.begin(); sample_iter != vcf_samples.end(); sample_iter++)
    sample_indices[*sample_iter] = sample_count++;

  std::vector< std::set<int32_t> > bad_sites_by_family(tracker != NULL ? tracker->families().size() : 0);

  // Iterate through all VCF entries
  std::vector< std::vector<SNP> > snps_by_sample(vcf_samples.size());
  VCF::Variant variant;
  uint32_t locus_count = 0;
  while (snp_vcf->get_next_variant(variant)){
    if (!variant.is_biallelic_snp() || in_any_region(variant, skip_regions, skip_padding))
      continue;

    // When performing pedigree-based filtering, we need to identify sites with any Mendelian
    // inconsistencies or missing genotypes as these won't be detected by the haplotype tracker
    if (tracker != NULL){
      const std::vector<NuclearFamily>& families = tracker->families();
      int family_index = 0;
      for (auto family_iter = families.begin(); family_iter != families.end(); ++family_iter, ++family_index)
	if (family_iter->is_missing_genotype(variant) || !family_iter->is_mendelian(variant))
	  bad_sites_by_family[family_index].insert(variant.get_position());
    }

    ++locus_count;
    int gt_a, gt_b;
    for (int i = 0; i < vcf_samples.size(); i++){
      if (variant.sample_call_missing(i) || !variant.sample_call_phased(i))
	continue;
      variant.get_genotype(i, gt_a, gt_b);
      if (gt_a != gt_b){
	char a1 = variant.get_allele(gt_a)[0];
	char a2 = variant.get_allele(gt_b)[0];

	// IMPORTANT NOTE: VCFs are 1-based, but BAMs are 0-based. Decrease VCF coordinate by 1 for consistency
	snps_by_sample[i].push_back(SNP(variant.get_position()-1, a1, a2));
      }
    }
  }
  logger << "Region contained a total of " << locus_count << " valid SNPs" << std::endl;

  // Filter out SNPs on a per-sample basis using any available pedigree information
  if (tracker != NULL){
    int32_t filt_count = 0, unfilt_count = 0;
    const std::vector<NuclearFamily>& families = tracker->families();
    int family_index = 0;
    for (auto family_iter = families.begin(); family_iter != families.end(); ++family_iter, ++family_index){
      std::vector<int> maternal_indices, paternal_indices;
      bool good_haplotypes = tracker->infer_haplotype_inheritance(*family_iter, DenovoScanner::MAX_BEST_SCORE, DenovoScanner::MIN_SECOND_BEST_SCORE,
								  maternal_indices, paternal_indices, bad_sites_by_family[family_index]);

      // If the family haplotypes aren't good enough, clear all of the sample's SNPs. Otherwise, remove only the bad sites from each sample's list
      for (auto sample_iter = family_iter->get_samples().begin(); sample_iter != family_iter->get_samples().end(); sample_iter++){
	auto sample_index = sample_indices.find(*sample_iter);
	if (sample_index != sample_indices.end()){
	  filt_count += snps_by_sample[sample_index->second].size();
	  if (!good_haplotypes)
	    snps_by_sample[sample_index->second].clear();
	  else
	    filter_snps(snps_by_sample[sample_index->second], bad_sites_by_family[family_index]);
	  filt_count   -= snps_by_sample[sample_index->second].size();
	  unfilt_count += snps_by_sample[sample_index->second].size();
	}
      }
    }
    logger << "Removed " << filt_count << " out of " << filt_count+unfilt_count << " individual heterozygous SNP calls due to pedigree uncertainties or inconsistencies" << std::endl;
  }

  // Create SNP trees
  for (unsigned int i = 0; i < snps_by_sample.size(); i++){
    //logger << "Building SNP tree for " << variant_file.sampleNames[i] << " containing " << snps_by_sample[i].size() << " heterozygous SNPs" << std::endl;
    snp_trees.push_back(new SNPTree(snps_by_sample[i]));
  }

  // Discard SNPs
  snps_by_sample.clear();
  
  return true;
}

void destroy_snp_trees(std::vector<SNPTree*>& snp_trees){
  for (unsigned int i = 0; i < snp_trees.size(); i++)
    delete snp_trees[i];
  snp_trees.clear();
}
