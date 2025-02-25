#include <climits>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#include "bam_io.h"
#include "error.h"
#include "genotyper_bam_processor.h"
#include "pedigree.h"
#include "stringops.h"
#include "vcf_reader.h"
#include "version.h"

bool file_exists(const std::string& path){
  return (access(path.c_str(), F_OK) != -1);
}

bool is_file(const std::string& name){
  struct stat st_buf;
  stat(name.c_str(), &st_buf);
  return (S_ISREG (st_buf.st_mode));
}

void print_usage(int def_mdist, int def_min_reads, int def_max_reads, int def_max_str_len, int def_max_haplotypes, int def_max_flanks, double def_min_flank_freq, int def_indel_flank_len, int def_switch_old_align_len, double def_min_mapq, double def_mean_qual){
  std::cerr << "Usage: LongTR --bams <list_of_bams> --fasta <genome.fa> --regions <region_file.bed> --tr-vcf <tr_gts.vcf.gz> [OPTIONS]" << "\n" << "\n"
    
	    << "Required parameters:" << "\n"
	    << "\t" << "--bams          <list_of_bams>        "  << "\t" << "Comma separated list of BAM/CRAM files. Either --bams or --bam-files must be specified"   << "\n"
	    << "\t" << "--fasta         <genome.fa>           "  << "\t" << "FASTA file containing all of the relevant sequences for your organism   "                 << "\n"
	    << "\t" << "                                      "  << "\t" << "When analyzing CRAMs, this FASTA file must match the file used for compression"         << "\n"
	    << "\t" << "--regions       <region_file.bed>     "  << "\t" << "BED file containing coordinates for each TR region"                                      << "\n"
	    << "\t" << "--tr-vcf        <tr_gts.vcf.gz>       "  << "\t" << "Bgzipped VCF file to which TR genotypes will be written"                                 << "\n" << "\n"

	    << "Optional input parameters:" << "\n"
	    << "\t" << "--bam-files  <bam_files.txt>          "  << "\t" << "File containing BAM/CRAM files to analyze, one per line "                              << "\n"
	    << "\t" << "--ref-vcf    <tr_ref_panel.vcf.gz>    "  << "\t" << "Bgzipped input VCF file of a reference panel of TR genotypes. VCF alleles will be"    << "\n"
	    << "\t" << "                                      "  << "\t" << " used as candidate variants instead of finding candidates in the BAMs/CRAMs (Default)" << "\n"
	    << "\t" << "--snp-vcf    <phased_snps.vcf.gz>     "  << "\t" << "Bgzipped input VCF file containing phased SNP genotypes for the samples"               << "\n"
	    << "\t" << "                                      "  << "\t" << " to be genotyped. These SNPs will be used to physically phase TRs "                   << "\n"

	    << "\t" << "--min-mean-qual       <threshold>     "  << "\t" << "Minimum average quality threshold for sequencing read data (Default = " << def_mean_qual << ")" << "\n"
	    << "\t" << "--min-mapq            <threshold>     "  << "\t" << "Minimum MAPQ per read (Default = " << def_min_mapq << ")" << "\n"
	    << "\t" << "--stutter-align-len   <threshold>     "  << "\t" << "Use stutter alignment for repeats with length less than threshold (Default = " << def_switch_old_align_len << ")" << "\n"
	    << "\t" << "--phased-bam	                      "  << "\t" << "Use phasing information from haplotagged sequencing data." << "\n"
	    << "\t" << "--indel-flank-len     <max_bp>        "  << "\t" << "Include InDels in max_bp base pair around repeat as part of the repeath (Default = " << def_indel_flank_len << ")" << "\n"
	    << "\t" << "--alignment-params    <list_of_params>"  << "\t" << "Comma-separated parameters for alignment model. Refer to GitHub for more details." << "\n"
	    << "\t" << "--stutter-in <stutter_models.txt>     "  << "\t" << "Use stutter models in the file to genotype STRs (Default = Learn via EM algorithm)"    << "\n" << "\n"

	    << "Optional output parameters:" << "\n"
	    << "\t" << "--log           <log.txt>             "  << "\t" << "Output the log information to the provided file (Default = Standard error)"         << "\n"
	    << "\t" << "--viz-out       <aln_viz.gz>          "  << "\t" << "Output a file of each locus' alignments for visualization with VizAln or VizAlnPdf" << "\n"
	    << "\t" << "--stutter-out   <stutter_models.txt>  "  << "\t" << "Output stutter models learned by the EM algorithm to the provided file"             << "\n" << "\n"
        << "\t" << "--pass-bam      <used_reads.bam>      "  << "\t" << "Output a BAM file containing the reads used to genotype each region"                 << "\n"
        << "\t" << "--filt-bam      <filt_reads.bam>      "  << "\t" << "Output a BAM file containing the reads filtered in each region. Each BAM entry"      << "\n"
        << "\t" << "                                      "  << "\t" << " has an FT tag specifying the reason for filtering"                                  << "\n" << "\n"

	    << "Optional VCF formatting parameters:" << "\n"
	    << "\t" << "--max-flank-indel <max_flank_frac>    "  << "\t" << "Don't output genotypes for a sample if the fraction of reads containing an indel"    << "\n"
	    << "\t" << "                                      "  << "\t" << " in the sequence flanking the STR is greater than MAX_FLANK_FRAC (Default = 0.15)"   << "\n"
        << "\t" << "--hide-allreads                       "  << "\t" << "Don't output the ALLREADS  FORMAT field to the VCF. By default, it will be output"   << "\n"
        << "\t" << "--hide-mallreads                      "  << "\t" << "Don't output the MALLREADS FORMAT field to the VCF. By default, it will be output"   << "\n"
	    << "\t" << "--output-gls                          "  << "\t" << "Write genotype likelihoods to the VCF (Default = False)"                             << "\n"
	    << "\t" << "--output-pls                          "  << "\t" << "Write phred-scaled genotype likelihoods to the VCF (Default = False)"                << "\n"
	    << "\t" << "--output-phased-gls                   "  << "\t" << "Write phased genotype likelihoods to the VCF (Default = False)"                      << "\n"
	    << "\t" << "--output-filters                      "  << "\t" << "Write why individual calls were filtered to the VCF (Default = False)"               << "\n" << "\n"

	    << "Optional BAM/CRAM tweaking parameters:" << "\n"
	    << "\t" << "--bam-samps     <list_of_samples>     "  << "\t" << "Comma separated list of read groups in same order as BAM/CRAM files. "               << "\n"
	    << "\t" << "                                      "  << "\t" << "  Assign each read the read group corresponding to its file. By default, "           << "\n"
	    << "\t" << "                                      "  << "\t" << "  each read must have an RG tag and the sample is determined from the SM field"      << "\n"
	    << "\t" << "--bam-libs      <list_of_libraries>   "  << "\t" << "Comma separated list of libraries in same order as BAM/CRAM files. "                 << "\n"
	    << "\t" << "                                      "  << "\t" << "  Assign each read the library corresponding to its file. By default, "              << "\n"
	    << "\t" << "                                      "  << "\t" << "  each read must have an RG tag and the library is determined from the LB field"     << "\n"
	    << "\t" << "--lib-from-samp                       "  << "\t" << "Assign each read the library corresponding to its sample name. By default,  "       << "\n"
	    << "\t" << "                                      "  << "\t" << "  each read must have an RG tag and the library is determined from the LB field"     << "\n" << "\n"

	    << "Optional haplotype filtering parameters:" << "\n"
	    << "\t" << "--max-haps <max_haplotypes>           "  << "\t" << "Maximum allowable candidate haplotypes for an TR (Default = " << def_max_haplotypes << ")" << "\n"
	    << "\t" << "                                      "  << "\t" << " Loci with more candidate haplotypes will not be genotyped" << "\n"
	    << "\t" << "--max-hap-flanks <max_flanks>         "  << "\t" << "Maximum allowable non-reference flanking sequences for an TR (Default = " << def_max_flanks << ")" << "\n"
	    << "\t" << "                                      "  << "\t" << " Loci with more candidate flanks will not be genotyped"                              << "\n"
	    << "\t" << "--min-flank-freq <min_freq>           "  << "\t" << "Filter a flank if its fraction of supporting samples < MIN_FREQ (Default = " << def_min_flank_freq  << ")" << "\n" << "\n"

	    << "Other optional parameters:" << "\n"
	    << "\t" << "--help                                "  << "\t" << "Print this help message and exit"                                                     << "\n"
	    << "\t" << "--version                             "  << "\t" << "Print LongTR version and exit"                                                        << "\n"
	    << "\t" << "--quiet                               "  << "\t" << "Only output terse logging messages (Default = output all messages)"                   << "\n"
	    << "\t" << "--silent                              "  << "\t" << "Don't output any logging messages  (Default = output all messages)"                   << "\n"
	    << "\t" << "--def-stutter-model                   "  << "\t" << "For each locus, use a stutter model with PGEOM=0.9 and UP=DOWN=0.05 for in-frame"     << "\n"
	    << "\t" << "                                      "  << "\t" << " artifacts and PGEOM=0.9 and UP=DOWN=0.01 for out-of-frame artifacts"                 << "\n"
	    << "\t" << "--chrom              <chrom>          "  << "\t" << "Only consider TRs on this chromosome"                                                << "\n"
	    << "\t" << "--haploid-chrs       <list_of_chroms> "  << "\t" << "Comma separated list of chromosomes to treat as haploid (Default = all diploid)"      << "\n"
	    << "\t" << "--hap-chr-file       <hap_chroms.txt> "  << "\t" << "File containing chromosomes to treat as haploid, one per line"                        << "\n"
	    << "\t" << "--min-reads          <num_reads>      "  << "\t" << "Minimum total reads required to genotype a locus (Default = " << def_min_reads << ")" << "\n"
	    << "\t" << "--max-reads          <num_reads>      "  << "\t" << "Skip a locus if it has more than NUM_READS reads (Default = " << def_max_reads << ")" << "\n"
	    << "\t" << "--max-tr-len         <max_bp>         "  << "\t" << "Only genotype TRs in the provided BED file with length < MAX_BP (Default = " << def_max_str_len << ")" << "\n" << std::endl;
    //<< "\t" << "--skip-genotyping                     "  << "\t" << "Don't perform any STR genotyping and merely compute the stutter model for each STR"  << "\n"
    //<< "\t" << "--dont-use-all-reads                  "  << "\t" << "Only utilize the reads HipSTR thinks will be informative for genotyping"   << "\n"
    //<< "\t" << "                                      "  << "\t" << " Enabling this option usually slightly decreases accuracy but shortens runtimes (~2x)"      << "\n"
    //<< "\t" << "--read-qual-trim     <min_qual>       "  << "\t" << "Trim both ends of a read until a base has quality score > MIN_QUAL (Default = 5)"    << "\n"
	//    << "\t" << "--fam <fam_file>                      "  << "\t" << "FAM file containing pedigree information for samples of interest. Use the pedigree"  << "\n"
     //       << "\t" << "                                      "  << "\t" << "  information to filter SNPs prior to phasing TRs (Default = use all SNPs)"         << "\n"

	    std::cout << "\n" << "\n"
	    << "\t A more detailed documentation of LongTR is available at https://github.com/gymrek-lab/LongTR"  << "\n"
	    << "\t Please file an issue on GitHub if you find a bug/issue or have a feature request." << "\n" << std::endl;
}

void parse_command_line_args(int argc, char** argv,
			     std::string& bamfile_string,     std::string& bamlist_string,    std::string& rg_sample_string,  std::string& rg_lib_string,
			     std::string& haploid_chr_string, std::string& hap_chr_file,      std::string& fasta_file,        std::string& region_file,   std::string& snp_vcf_file,
			     std::string& chrom,              std::string& bam_pass_out_file, std::string& bam_filt_out_file, std::string& ref_vcf_file,
			     std::string& str_vcf_out_file,   std::string& fam_file,          std::string& log_file,
			     int& bam_lib_from_samp, int& skip_genotyping, GenotyperBamProcessor& bam_processor, std::string& alignment_params_string){
  int def_mdist             = bam_processor.MAX_MATE_DIST;
  int def_min_reads         = bam_processor.MIN_TOTAL_READS;
  int def_max_reads         = bam_processor.MAX_TOTAL_READS;
  int def_max_str_len       = bam_processor.MAX_STR_LENGTH;
  int def_max_flanks        = bam_processor.MAX_FLANK_HAPLOTYPES;
  int def_max_haplotypes    = bam_processor.MAX_TOTAL_HAPLOTYPES;
  double def_min_flank_freq = bam_processor.MIN_FLANK_FREQ;
  int def_indel_flank_len = bam_processor.INDEL_FLANK_LEN;
  int def_switch_old_align_len = bam_processor.SWITCH_OLD_ALIGN_LEN;
  double def_min_mapq = bam_processor.MIN_MAPQ;
  double def_mean_qual = bam_processor.MIN_SUM_QUAL_LOG_PROB;

  if (argc == 1 || (argc == 2 && std::string("-h").compare(std::string(argv[1])) == 0)){
    print_usage(def_mdist, def_min_reads, def_max_reads, def_max_str_len, def_max_haplotypes, def_max_flanks, def_min_flank_freq, def_indel_flank_len, def_switch_old_align_len, def_min_mapq, def_mean_qual);
    exit(0);
  }

  int print_help = 0, print_version = 0, quiet_log = 0, silent_log = 0, def_stutter_model = 1, phased_bam = 0, skip_assembly = 1, indel_flank_len = 5, switch_old_align_len = 0;

  static struct option long_options[] = {
    {"bams",            required_argument, 0, 'b'},
    {"bam-files",       required_argument, 0, 'B'},
    {"chrom",           required_argument, 0, 'c'},
    {"max-mate-dist",   required_argument, 0, 'd'},
    {"fam",             required_argument, 0, 'D'},
    {"fasta",           required_argument, 0, 'f'},
    {"bam-samps",       required_argument, 0, 'g'},
    {"max-hap-flanks",  required_argument, 0, 'G'},
    {"max-haps",        required_argument, 0, 'J'},
    {"bam-libs",        required_argument, 0, 'q'},
    {"min-reads",       required_argument, 0, 'i'},
    {"min-flank-freq",  required_argument, 0, 'I'},
    {"read-qual-trim",  required_argument, 0, 'j'},
    {"log",             required_argument, 0, 'l'},
    {"max-reads",       required_argument, 0, 'n'},
    {"max-flank-indel", required_argument, 0, 'F'},
    {"tr-vcf",         required_argument, 0, 'o'},
    {"ref-vcf",         required_argument, 0, 'p'},
    {"regions",         required_argument, 0, 'r'},
    {"snp-vcf",         required_argument, 0, 'v'},
    {"stutter-in",      required_argument, 0, 'm'},
    {"stutter-out",     required_argument, 0, 's'},
    {"sample-list",     required_argument, 0, 'S'},
    {"haploid-chrs",    required_argument, 0, 't'},
    {"hap-chr-file",    required_argument, 0, 'u'},
    {"pass-bam",        required_argument, 0, 'w'},
    {"max-tr-len",     required_argument, 0, 'x'},
    {"filt-bam",        required_argument, 0, 'y'},
    {"viz-out",         required_argument, 0, 'z'},
    {"min-mean-qual",	required_argument, 0, 'W'},
    {"min-mapq",	    required_argument, 0, 'A'},
    {"phased-bam",           no_argument, &phased_bam, 1},
    {"h",                  no_argument, &print_help, 1},
    {"help",               no_argument, &print_help, 1},
    {"lib-from-samp",      no_argument, &bam_lib_from_samp, 1},
    {"hide-allreads",      no_argument, &(Genotyper::OUTPUT_ALLREADS),      0},
    {"hide-mallreads",     no_argument, &(Genotyper::OUTPUT_MALLREADS),     0},
    {"output-gls",         no_argument, &(Genotyper::OUTPUT_GLS),           1},
    {"output-pls",         no_argument, &(Genotyper::OUTPUT_PLS),           1},
    {"output-phased-gls",  no_argument, &(Genotyper::OUTPUT_PHASED_GLS),    1},
    {"output-filters",     no_argument, &(Genotyper::OUTPUT_FILTERS),       1},
    {"no-rmdup",           no_argument, &(bam_processor.REMOVE_PCR_DUPS),      0},
    {"use-unpaired",       no_argument, &(bam_processor.REQUIRE_PAIRED_READS), 0},
    {"dont-use-all-reads", no_argument, &(bam_processor.REQUIRE_SPANNING),     1},
    {"viz-left-alns",      no_argument, &(bam_processor.VIZ_LEFT_ALNS),        1},
    {"def-stutter-model",  no_argument, &def_stutter_model, 1},
    {"version",            no_argument, &print_version, 1},
    {"quiet",              no_argument, &quiet_log, 1},
    {"silent",             no_argument, &silent_log, 1},
    {"skip-genotyping",    no_argument, &skip_genotyping, 1},
    {"skip-assembly",	   no_argument, &skip_assembly, 0},
    {"stutter-align-len",    required_argument, 0, 'O'},
    {"indel-flank-len",    required_argument, 0, 'L'},
    {"alignment-params", required_argument, 0, 'P'},
    {0, 0, 0, 0}
  };

  std::string filename;
  while (true){
    int option_index = 0;
    int c = getopt_long(argc, argv, "b:B:c:d:D:e:f:F:g:G:i:I:j:k:l:m:n:o:p:q:r:s:S:t:u:v:w:x:y:z:W:O:L:A:P", long_options, &option_index);
    if (c == -1)
      break;

    if (optarg != NULL){
      std::string val(optarg);
      if (string_starts_with(val, "--"))
	printErrorAndDie("Argument to option --" + std::string(long_options[option_index].name) + " cannot begin with \"--\"\n\tBad argument: " + val);
    }

    switch(c){
    case 0:
      break;
    case 'b':
      bamlist_string = std::string(optarg);
      break;
    case 'B':
      bamfile_string = std::string(optarg);
      break;
    case 'c':
      chrom = std::string(optarg);
      break;
    case 'd':
      bam_processor.MAX_MATE_DIST = atoi(optarg);
      break;
    case 'D':
      fam_file = std::string(optarg);
      break;
    case 'f':
      fasta_file = std::string(optarg);
      break;
    case 'g':
      rg_sample_string = std::string(optarg);
      break;
    case 'G':
      bam_processor.MAX_FLANK_HAPLOTYPES = atoi(optarg);
      if (bam_processor.MAX_FLANK_HAPLOTYPES < 1)
	printErrorAndDie("--max-hap-flanks must be greater than 0");
      break;
    case 'i':
      bam_processor.MIN_TOTAL_READS = atoi(optarg);
      if (bam_processor.MIN_TOTAL_READS < 0)
	printErrorAndDie("--min-total-reads must be >= 0");
      break;
    case 'I':
      bam_processor.MIN_FLANK_FREQ = atof(optarg);
      if (bam_processor.MIN_FLANK_FREQ < 0 || bam_processor.MIN_FLANK_FREQ > 1)
	printErrorAndDie("--min-flank-freq must be between 0 and 1");
      break;
    case 'j':
      if (std::string(optarg).size() != 1)
	printErrorAndDie("--read-qual-trim requires a single character argument");
      bam_processor.BASE_QUAL_TRIM = std::string(optarg)[0];
      break;
    case 'J':
      bam_processor.MAX_TOTAL_HAPLOTYPES = atoi(optarg);
      if (bam_processor.MAX_TOTAL_HAPLOTYPES <= 1)
	printErrorAndDie("--max-haps must be greater than 1");
      break;
    case 'l':
      log_file = std::string(optarg);
      break;
    case 'm':
      filename = std::string(optarg);
      bam_processor.set_input_stutter(filename);
      break;
    case 'n':
      bam_processor.MAX_TOTAL_READS = atoi(optarg);
      break;
    case 'o':
      str_vcf_out_file = std::string(optarg);
      break;
    case 'p':
      ref_vcf_file = std::string(optarg);
      break;
    case 'q':
      rg_lib_string = std::string(optarg);
      break;
    case 'r':
      region_file = std::string(optarg);
      break;
    case 's':
      filename = std::string(optarg);
      bam_processor.set_output_stutter(filename);
      break;
    case 'S':
      bam_processor.set_sample_set(std::string(optarg));
      break;
    case 't':
      haploid_chr_string = std::string(optarg);
      break;
    case 'u':
      hap_chr_file = std::string(optarg);
      break;
    case 'v':
      snp_vcf_file = std::string(optarg);
      break;
    case 'w':
      bam_pass_out_file = std::string(optarg);
      break;
    case 'x':
      bam_processor.MAX_STR_LENGTH = atoi(optarg);
      break;
    case 'y':
      bam_filt_out_file = std::string(optarg);
      break;
    case 'z':
      filename = std::string(optarg);
      if (!string_ends_with(filename, ".gz"))
	printErrorAndDie("Path for alignment visualization file must end in .gz as it will be bgzipped");
      bam_processor.set_output_viz(filename);
      break;
    case 'F':
      Genotyper::MAX_FLANK_INDEL_FRAC = atof(optarg);
      break;
    case 'W':
	bam_processor.MIN_SUM_QUAL_LOG_PROB = atof(optarg);
	break;
	case 'A':
	bam_processor.MIN_MAPQ = atof(optarg);
	break;
	case 'L':
	bam_processor.INDEL_FLANK_LEN = atof(optarg);
	break;
	case 'O':
	bam_processor.SWITCH_OLD_ALIGN_LEN = atof(optarg);
	break;
    case 'P':
    alignment_params_string = std::string(optarg);
    break;
    case '?':
      printErrorAndDie("Unrecognized command line option");
      break;
    default:
      abort();
      break;
    }
  }
  if (optind < argc) {
    std::stringstream msg;
    msg << "Did not recognize the following command line arguments:" << "\n";
    while (optind < argc)
      msg << "\t" << argv[optind++] << "\n";
    msg << "Please check your command line syntax or type ./LongTR--help for additional information" << "\n";
    printErrorAndDie(msg.str());
  }

  if (print_version == 1){
    std::cerr << "LongTR version " << VERSION << std::endl;
    exit(0);
  }
  if (print_help){
    print_usage(def_mdist, def_min_reads, def_max_reads, def_max_str_len, def_max_haplotypes, def_max_flanks, def_min_flank_freq, def_indel_flank_len, def_switch_old_align_len, def_min_mapq, def_mean_qual);
    exit(0);
  }
  if (quiet_log)
    bam_processor.suppress_most_logging();
  if (silent_log)
    bam_processor.suppress_all_logging();
  if (def_stutter_model == 1)
    bam_processor.set_default_stutter_model(0.95, 0.05, 0.05, 0.95, 0.01, 0.01);
  if (phased_bam){
    bam_processor.use_phased_bam_tags();
    bam_processor.full_logger() << "Using phased BAM tags to genotype and phase TRs (WARNING: Any arguments provided to --snp-vcf will be ignored)" << std::endl;
  }
  if (skip_assembly == 1){
    bam_processor.skip_assembly();
  }
}	

int main(int argc, char** argv){
  double total_time = clock();
  precompute_integer_logs(); // Calculate and cache log of integers from 1 -> 999

  std::stringstream full_command_ss;
  full_command_ss << "LongTR-" << VERSION;
  for (int i = 1; i < argc; i++)
    full_command_ss << " " << argv[i];
  std::string full_command = full_command_ss.str();

  GenotyperBamProcessor bam_processor(true, false);

  int bam_lib_from_samp = 0, skip_genotyping = 0;
  std::string bamfile_string="", bamlist_string="", rg_sample_string="", rg_lib_string="", hap_chr_string="", hap_chr_file="", alignment_params_string="";
  std::string region_file="", fasta_file="", chrom="", snp_vcf_file="";
  std::string bam_pass_out_file="", bam_filt_out_file="", str_vcf_out_file="", fam_file = "", log_file = "", ref_vcf_file="";

  parse_command_line_args(argc, argv,
			  bamfile_string, bamlist_string, rg_sample_string, rg_lib_string, hap_chr_string, hap_chr_file, fasta_file, region_file, snp_vcf_file,
			  chrom, bam_pass_out_file, bam_filt_out_file, ref_vcf_file, str_vcf_out_file, fam_file, log_file, bam_lib_from_samp, skip_genotyping, bam_processor, alignment_params_string);

  if (!log_file.empty())
    bam_processor.set_log(log_file);

  if (bamfile_string.empty() && bamlist_string.empty())
    printErrorAndDie("You must specify either the --bams or --bam-files option");
  else if ((!bamfile_string.empty()) && (!bamlist_string.empty()))
    printErrorAndDie("You can only specify one of the --bams or --bam-files options");
  else if (region_file.empty()){
    std::stringstream err;
    err << "--regions option required" << "\n";
    printErrorAndDie(err.str());
  }
  else if (fasta_file.empty())
    printErrorAndDie("--fasta option required");
  else if (!skip_genotyping && str_vcf_out_file.empty())
    printErrorAndDie("--tr-vcf option required");

  // Check that the FASTA file exists
  if (!file_exists(fasta_file) || !is_file(fasta_file))
    printErrorAndDie("FASTA file " + fasta_file + " does not exist. Please ensure that the path provided to --fasta is a valid FASTA file");

  std::vector<std::string> bam_files;
  if (!bamlist_string.empty()) split_by_delim(bamlist_string, ',', bam_files);


  // Checking and parsing alignment parameters
  std::vector<std::string> alignment_params;
  if (!alignment_params_string.empty()){
    split_by_delim(alignment_params_string, ',', alignment_params);
    if (alignment_params.size() != 7) printErrorAndDie("Number of alignment parameters is not correct");
    std::vector<float> alignment_params_float;

    for (auto p : alignment_params) {
        try {
            float f = std::stof(p);
            if (f >= 0)
                printErrorAndDie("LOG value " + p + " can not be positive.");
            alignment_params_float.push_back(f); // Convert and add to the list
        } catch (const std::invalid_argument& e) {
            printErrorAndDie("Invalid argument: " + p + " could not be converted to float.");
        } catch (const std::out_of_range& e) {
            printErrorAndDie("Out of range: " + p + " is out of the range for float.");
        }
    }
    bam_processor.set_alignment_params(alignment_params_float);
  }

  // Read file containing BAM/CRAM paths line-by-line
  if (!bamfile_string.empty()){
    if (!file_exists(bamfile_string))
      printErrorAndDie("File containing BAM/CRAM file names does not exist: " + bamfile_string);
    std::ifstream input(bamfile_string.c_str());
    if (!input.is_open())
      printErrorAndDie("Failed to open file containing BAM/CRAM file names: " + bamfile_string);
    std::string line;
    while (std::getline(input, line))
      if (!line.empty())
	bam_files.push_back(line);
    input.close();
  }
  bam_processor.full_logger() << "Detected " << bam_files.size() << " BAM/CRAM files" << std::endl;

  // Open all BAM/CRAM files
  std::string cram_fasta_path = fasta_file;
  int merge_type = BamCramMultiReader::ORDER_ALNS_BY_FILE;
  BamCramMultiReader reader(bam_files, cram_fasta_path, merge_type);

  // Construct filename->read group map (if one has been specified) and determine the list
  // of samples of interest based on either the specified names or the RG tags in the BAM/CRAM headers
  std::set<std::string> rg_samples, rg_libs;
  std::map<std::string, std::string> rg_ids_to_sample, rg_ids_to_library;
  if (!rg_sample_string.empty()){
    if ((bam_lib_from_samp == 0) && rg_lib_string.empty())
      printErrorAndDie("--bam-libs option required when --bam-samps option specified");

    std::vector<std::string> read_groups, libraries;
    split_by_delim(rg_sample_string, ',', read_groups);
    split_by_delim(rg_lib_string, ',', libraries);
    if (bam_files.size() != read_groups.size())
      printErrorAndDie("Number of BAM/CRAM files in --bams and samples in --bam-samps must match");
    if ((bam_lib_from_samp == 0) && (bam_files.size() != libraries.size()))
      printErrorAndDie("Number of BAM/CRAM files in --bams and libraries in --bam-libs must match");

    for (unsigned int i = 0; i < bam_files.size(); i++){
      rg_ids_to_sample[bam_files[i]]  = read_groups[i];
      rg_ids_to_library[bam_files[i]] = (bam_lib_from_samp == 0 ? libraries[i]: read_groups[i]);
      rg_samples.insert(read_groups[i]);
    }
    bam_processor.use_custom_read_groups();
    bam_processor.full_logger() << "User-specified read groups for " << rg_samples.size() << " unique samples" << std::endl;
  }
  else {
    for (unsigned int i = 0; i < bam_files.size(); i++){
      const std::vector<ReadGroup>& read_groups = reader.bam_header()->read_groups(i);
      if (read_groups.empty())
	printErrorAndDie("Provided BAM/CRAM files don't contain read groups in the header and the --bam-samps flag was not specified");

      for (auto rg_iter = read_groups.begin(); rg_iter != read_groups.end(); rg_iter++){
	if (!rg_iter->HasID())     printErrorAndDie("RG in BAM/CRAM header is lacking the ID tag");
	if (!rg_iter->HasSample()) printErrorAndDie("RG in BAM/CRAM header is lacking the SM tag");
	if ((bam_lib_from_samp == 0) && !rg_iter->HasLibrary())
	  printErrorAndDie("RG in BAM/CRAM header is lacking the LB tag");
	std::string rg_library = (bam_lib_from_samp == 0 ? rg_iter->GetLibrary() : rg_iter->GetSample());

	// Ensure that there aren't identical read group ids that map to different samples or libraries
	if (rg_ids_to_sample.find(rg_iter->GetID()) != rg_ids_to_sample.end())
	  if (rg_ids_to_sample[rg_iter->GetID()].compare(rg_iter->GetSample()) != 0)
	    printErrorAndDie("Read group id " + rg_iter->GetID() + " maps to more than one sample");
	if (rg_ids_to_library.find(rg_iter->GetID()) != rg_ids_to_library.end())
	  if (rg_ids_to_library[rg_iter->GetID()].compare(rg_library) != 0)
	    printErrorAndDie("Read group id " + rg_iter->GetID() + " maps to more than one library");

	rg_ids_to_sample[bam_files[i] + rg_iter->GetID()]  = rg_iter->GetSample();
	rg_ids_to_library[bam_files[i] + rg_iter->GetID()] = rg_library;
	rg_samples.insert(rg_iter->GetSample());
	rg_libs.insert(rg_library);
      }
    }

    bam_processor.full_logger() << "BAMs/CRAMs contain unique read group IDs for "
				<< rg_libs.size()    << " unique libraries and "
				<< rg_samples.size() << " unique samples" << std::endl;
  }

  BamWriter* bam_pass_writer = NULL;
  if (!bam_pass_out_file.empty())
    bam_pass_writer = new BamWriter(bam_pass_out_file, reader.bam_header());

  BamWriter* bam_filt_writer = NULL;
  if (!bam_filt_out_file.empty())
    bam_filt_writer = new BamWriter(bam_filt_out_file, reader.bam_header());

  if (!ref_vcf_file.empty()){
    if (!string_ends_with(ref_vcf_file, ".gz"))
      printErrorAndDie("Ref VCF file must be bgzipped (and end in .gz)");

    // Check that the VCF exists
    if (!file_exists(ref_vcf_file)) 
      printErrorAndDie("Ref VCF file " + ref_vcf_file + " does not exist. Please ensure that the path provided to --ref-vcf is valid");

    // Check that tabix index exists
    if (!file_exists(ref_vcf_file + ".tbi"))
	printErrorAndDie("No .tbi index found for the ref VCF file. Please index using tabix and rerun LongTR");

    bam_processor.set_ref_vcf(ref_vcf_file);
  }

  if (!snp_vcf_file.empty()){
    if (!string_ends_with(snp_vcf_file, ".gz"))
      printErrorAndDie("SNP VCF file must be bgzipped (and end in .gz)");
    
    // Check that the VCF exists
    if (!file_exists(snp_vcf_file))
      printErrorAndDie("SNP VCF file " + snp_vcf_file + " does not exist. Please ensure that the path provided to --snp-vcf is valid");

    // Check that tabix index exists
    if (!file_exists(snp_vcf_file + ".tbi"))
	printErrorAndDie("No .tbi index found for the SNP VCF file. Please index using tabix and rerun LongTR");

    bam_processor.set_input_snp_vcf(snp_vcf_file);
  }

  if (!skip_genotyping){
    if (!string_ends_with(str_vcf_out_file, ".gz"))
      printErrorAndDie("Path for TR VCF output file must end in .gz as it will be bgzipped");
    bam_processor.set_output_str_vcf(str_vcf_out_file, fasta_file, full_command, rg_samples);
  }

  if (!hap_chr_string.empty()){
    std::vector<std::string> haploid_chroms;
    split_by_delim(hap_chr_string, ',', haploid_chroms);
    for (auto chrom_iter = haploid_chroms.begin(); chrom_iter != haploid_chroms.end(); chrom_iter++)
      bam_processor.add_haploid_chrom(*chrom_iter);
  }
  if (!hap_chr_file.empty()){
    if (!file_exists(hap_chr_file))
      printErrorAndDie("File containing haploid chromosome names does not exist: " + hap_chr_file);
    std::ifstream input(hap_chr_file.c_str());
    if (!input.is_open())
      printErrorAndDie("Failed to open file containing haploid chromosome names: " + hap_chr_file);
    std::string line;
    while (std::getline(input, line))
      if (!line.empty())
	bam_processor.add_haploid_chrom(line);
    input.close();
  }

  // Extract any relevant pedigree information to be used to filter SNPs before physically phasing STRs
  if (!fam_file.empty()){
    if (snp_vcf_file.empty())
      printErrorAndDie("--fam option only applies if --snp-vcf option has been specified as well");

    // Determine what samples are in the SNP VCF
    VCF::VCFReader snp_vcf(snp_vcf_file);
    std::set<std::string> samples_with_data(snp_vcf.get_samples().begin(), snp_vcf.get_samples().end());

    std::vector<NuclearFamily> families;
    extract_pedigree_nuclear_families(fam_file, samples_with_data, families, bam_processor.full_logger());
    if (families.size() != 0)
      bam_processor.use_pedigree_to_filter_snps(families, snp_vcf_file);
  }

  // Run analysis
  bam_processor.process_regions(reader, region_file, fasta_file, rg_ids_to_sample, rg_ids_to_library, full_command, bam_pass_writer, bam_filt_writer, 10000000, chrom);
  bam_processor.finish();

  if (bam_pass_writer != NULL) delete bam_pass_writer;
  if (bam_filt_writer != NULL) delete bam_filt_writer;


  total_time = (clock() - total_time)/CLOCKS_PER_SEC;
  bam_processor.full_logger() << "LongTR execution finished: Total runtime = " << total_time << " sec" << "\n"
			      << "-----------------\n\n" << std::endl;
  return 0;  
}
