import argparse
import collections
import sys

try:
     import vcf
except ImportError:
     exit("This script requires the PyVCF python package. Please install using pip or conda")

def filter_call(sample, filters):
     if sample['DP'] < filters.DEPTH:
          return "Depth"
     elif sample['Q'] < filters.QUAL:
          return "Quality"
     else:
          d_1, d_2 = map(float, sample['PDP'].split('|'))
          if d_1 == 0 or d_2 == 0:
               return "Allele depth"
          if min(d_1, d_2) < filters.ALLELE_DEPTH:
               return "Allele depth"
          elif min(1.0*d_1/d_2, 1.0*d_2/d_1) < filters.ALLELE_RATIO:
               return "Allele ratio"
          elif filters.FLANK_INDEL_FRAC < 1 and 1.0*sample['DFLANKINDEL']/sample['DP'] > filters.FLANK_INDEL_FRAC:
               return "Flank indels"
          elif filters.STUTTER_FRAC < 1 and 1.0*sample['DSTUTTER']/sample['DP'] > filters.STUTTER_FRAC:
               return "Stutter fraction"
          elif filters.ALLELE_BIAS > -100 and sample['AB'] < filters.ALLELE_BIAS:
               return "Allele bias"
          elif filters.STRAND_BIAS > -100 and sample['FS'] < filters.STRAND_BIAS:
               return "Strand bias"

          if filters.SPAN_DEPTH > 0:
               if sample["MALLREADS"] is None:
                    return "Spanning depth"
               gb_a, gb_b = map(int, sample["GB"].split("|"))
               span_depth_dict = dict(map(lambda x: map(int, x.split("|")), sample["MALLREADS"].split(";")))
               dp_a = 0 if gb_a not in span_depth_dict else span_depth_dict[gb_a]
               dp_b = 0 if gb_b not in span_depth_dict else span_depth_dict[gb_b]
               if min(dp_a, dp_b) < filters.SPAN_DEPTH:
                    return "Spanning depth"
          return None

def main():
     parser = argparse.ArgumentParser()
     parser.add_argument("--vcf", help="Input VCF to filter", type=str, required=True, dest="VCF")

     parser.add_argument("--min-call-depth", type=int, required=False, default=0, dest="DEPTH",
                         help="Omit a sample's call if DP < DEPTH")

     parser.add_argument("--min-call-qual", type=float, required=False, default=0.0, dest="QUAL",
                         help="Omit a sample's call if Q < QUAL")

     parser.add_argument("--min-call-allele-depth", type=float, required=False, default=0.0, dest="ALLELE_DEPTH",
                         help="Omit a sample's call if the minimum allele depth (smallest value in PDP FORMAT field) < ALLELE_DEPTH")

     parser.add_argument("--min-call-depth-ratio", type=float, required=False, default=0.0, dest="ALLELE_RATIO",
                         help="Omit a sample's call if the minimum ratio of allele depths (values in PDP FORMAT field) < ALLELE_RATIO")                         

     parser.add_argument("--max-call-flank-indel", type=float, required=False, default=1.0, dest="FLANK_INDEL_FRAC",
                         help="Omit a sample's call if the fraction of reads with flank indels (FORMAT fields DFLANKINDEL/DP) > FLANK_INDEL_FRAC")
                         
     parser.add_argument("--max-call-stutter", type=float, required=False, default=1.0, dest="STUTTER_FRAC",
                         help="Omit a sample's call if the fraction of reads with a stutter artifact (FORMAT fields DSTUTTER/DP) > STUTTER_FRAC")
                         
     parser.add_argument("--min-call-allele-bias", type=float, required=False, default=-100.0, dest="ALLELE_BIAS",
                         help="Omit a sample's call if AB < ALLELE_BIAS")

     parser.add_argument("--min-call-strand-bias", type=float, required=False, default=-100.0, dest="STRAND_BIAS",
                         help="Omit a sample's call if FS < STRAND_BIAS")

     parser.add_argument("--min-call-spanning-depth", type=int, required=False, default=0.0, dest="SPAN_DEPTH",
                         help="Omit a sample's call if the minimum number of spanning reads supporting an allele < SPAN_DEPTH")

     parser.add_argument("--min-loc-depth", type=int, required=False, default=0, dest="MIN_LOC_DEPTH",
                         help="Omit locus if the total depth (DP INFO field) < MIN_LOC_DEPTH")

     parser.add_argument("--max-loc-depth", type=int, required=False, default=1000000000, dest="MAX_LOC_DEPTH",
                         help="Omit locus if the total depth (DP INFO field) > MAX_LOC_DEPTH")

     parser.add_argument("--max-loc-flank-indel", type=float, required=False, default=1.0, dest="LOC_FLANK_INDEL_FRAC",
                         help="Omit locus if the total fraction of reads with an indel in the flanks (INFO fields DFLANKINDEL/DP) > LOC_FLANK_INDEL_FRAC")

     parser.add_argument("--max-loc-stutter", type=float, required=False, default=1.0, dest="LOC_STUTTER",
                         help="Omit locus if the total fraction of reads with a stutter artifact (INFO fields DSTUTTER/DP) > LOC_STUTTER")

     parser.add_argument("--min-loc-calls",type=int,   required=False, default=0, dest="MIN_CALLS",
                         help="Omit locus if the number of valid samples after filtering < MIN_CALLS")
                         
     args = parser.parse_args()
     if args.VCF == "-":
          vcf_reader = vcf.Reader(sys.stdin)
     else:
          vcf_reader = vcf.Reader(filename=args.VCF)
     vcf_writer    = vcf.Writer(sys.stdout, vcf_reader)
     total_counts  = collections.defaultdict(int)
     filter_counts = {}
     for sample in vcf_reader.samples:
          filter_counts[sample] = collections.defaultdict(int)

     for record in vcf_reader:
          if record.INFO['DP'] < args.MIN_LOC_DEPTH:
               continue
          if record.INFO['DP'] > args.MAX_LOC_DEPTH:
               continue
          if args.LOC_FLANK_INDEL_FRAC < 1 and 1.0*record.INFO["DFLANKINDEL"]/record.INFO["DP"] > args.LOC_FLANK_INDEL_FRAC:
               continue
          if args.LOC_STUTTER < 1 and 1.0*record.INFO["DSTUTTER"]/record.INFO["DP"] > args.LOC_STUTTER:
               continue

          samp_fmt = vcf.model.make_calldata_tuple(record.FORMAT.split(':'))
          for fmt in samp_fmt._fields:
               entry_type = vcf_reader.formats[fmt].type
               entry_num  = vcf_reader.formats[fmt].num
               samp_fmt._types.append(entry_type)
               samp_fmt._nums.append(entry_num)

          # Determine if we can remove alleles with 0 counts based on whether any genotype-likelhood associated
          # format fields are present
          fmt_tokens         = record.FORMAT.split(":")
          can_remove_alleles = not ("GL" in fmt_tokens or "PL" in fmt_tokens or "PHASEDGL" in fmt_tokens)
               
          nfields       = len(samp_fmt._fields)
          num_filt      = 0
          allele_counts = len(record.alleles)*[0]

          for sample in record:
               if sample['GT'] is None or sample['GT'] == "./." or sample['GT'] == ".":
                    continue

               filter_reason = filter_call(sample, args)
               if filter_reason is None:
                    gt_a, gt_b = map(int, sample['GT'].split('|'))
                    allele_counts[gt_a] += 1
                    allele_counts[gt_b] += 1
               else:
                    filter_counts[sample.sample][filter_reason] += 1
                    total_counts[filter_reason] += 1

          # Build allele index mapping
          allele_indices   = {0:0}
          filt_num_alleles = 1
          for i in range(1, len(allele_counts)):
              if allele_counts[i] != 0 or not can_remove_alleles:
                  allele_indices[i] = filt_num_alleles
                  filt_num_alleles += 1

          new_samples       = []
          num_filt          = 0
          num_kept          = 0
          total_dp          = 0
          total_dstutter    = 0
          total_dflankindel = 0
          for sample in record:
               if sample['GT'] is None or sample['GT'] == "./." or sample['GT'] == ".":
                    new_samples.append(sample)
                    continue

               filter_reason = filter_call(sample, args)
               if filter_reason is not None:
                    num_filt += 1
                    sampdat   = []
                    for i in range(len(samp_fmt._fields)):
                         key = samp_fmt._fields[i]
                         if key == "GT":
                              sampdat.append("./.")
                         else:
                              if key == "FILTER":
                                   sampdat.append(filter_reason.replace(" ", "_").upper())
                              else:
                                   sampdat.append(None)
               else:
                    num_kept    += 1
                    gt_a, gt_b   = map(int, sample['GT'].split('|'))
                    new_gt       = "%d|%d"%(allele_indices[gt_a], allele_indices[gt_b])
                    sampdat      = []
                    for i in range(len(samp_fmt._fields)):
                        key = samp_fmt._fields[i]
                        if key == "GT":
                            sampdat.append(new_gt)
                        else:
                            sampdat.append(sample[key])

                    total_dp          += sample['DP']
                    total_dstutter    += sample['DSTUTTER']
                    total_dflankindel += sample['DFLANKINDEL']

               call = vcf.model._Call(record, sample.sample, samp_fmt(*sampdat))
               new_samples.append(call)
          record.samples = new_samples

          # Don't output the filtered record if we don't meet the minimum number of calls
          if num_kept < args.MIN_CALLS:
               continue
          
          # Fix set of alleles
          new_alleles = [record.alleles[0]]
          for i in range(1, len(record.alleles)):
              if allele_counts[i] != 0 or not can_remove_alleles:
                  new_alleles.append(record.alleles[i])

          record.alleles = new_alleles
          record.ALT     = new_alleles[1:]
          if len(record.ALT) == 0:
              record.ALT = ["."]
          if 'NFILT' in record.INFO:
               record.INFO['NFILT'] += num_filt

          record.INFO['DP']          = total_dp
          record.INFO['DSTUTTER']    = total_dstutter
          record.INFO['DFLANKINDEL'] = total_dflankindel

          # Recompute and set INFO fields
          if 'BPDIFFS' in record.INFO:
              if len(new_alleles) == 1:
                  record.INFO.pop("BPDIFFS", None)
              else:
                  record.INFO['BPDIFFS'] = list(map(lambda x: len(x) - len(new_alleles[0]), new_alleles[1:]))
          record.INFO['REFAC'] = allele_counts[0]
          if 'AC' in record.INFO:
              if len(new_alleles) == 1:
                  record.INFO.pop("AC", None)
              else:
                   if not can_remove_alleles:
                        record.INFO['AC'] = allele_counts[1:]
                   else:
                        record.INFO['AC'] = list(filter(lambda x: x != 0, allele_counts[1:]))
          if 'AN' in record.INFO:
               record.INFO['AN'] = sum(allele_counts)
                  
          vcf_writer.write_record(record)

if __name__ == "__main__":
    main()
