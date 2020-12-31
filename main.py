#! /usr/bin/env python

import subprocess, argparse
import func

def main(args):
        
    with open(args.position_bed_file, 'r') as hin, open(args.output_file + ".tmp.polished.fq", 'w') as hout:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            tchr, tstart, tend = F[0], int(F[1]), int(F[2])

            tread_num = func.get_read_segments(args.bam_file, args.output_file + ".tmp.seg.fa", 
                tchr, tstart, tend, args.read_num)
            
            # mafft based
            with open(args.output_file + ".tmp.seg.mafft.fa", 'w') as hout2:
                subprocess.check_call(["mafft", args.output_file + ".tmp.seg.fa"], 
                stdout = hout2, stderr = subprocess.DEVNULL)

            tconsensus = func.get_consensus_from_mafft_result(args.output_file + ".tmp.seg.mafft.fa")
            print("@%s_%s_%s_%d_mafft\n%s\n+\n%s" % 
                (tchr, tstart, tend, tread_num, tconsensus, "I" * len(tconsensus)), 
                file = hout)

            # racon based
            with open(args.output_file + ".tmp.seg.first.fa", 'w') as hout2:
                subprocess.check_call(["head", "-n", "2", args.output_file + ".tmp.seg.fa"], 
                stdout = hout2)

            """
            with open(args.output_file + ".tmp.seg.minimap2.paf", 'w') as hout2:
                subprocess.check_call(["minimap2", "-x", "map-ont", args.output_file + ".tmp.seg.first.fa",
                    args.output_file + ".tmp.seg.fa"], stdout = hout2, stderr = subprocess.DEVNULL)
            """
            func.generate_paf_file(args.output_file + ".tmp.seg.fa", 
                args.output_file + ".tmp.seg.first.fa",
                args.output_file + ".tmp.seg.parasail.paf")

            with open(args.output_file + ".tmp.seg.racon.fa", 'w') as hout2:
                subprocess.check_call(["racon", "-u", args.output_file + ".tmp.seg.fa",
                    args.output_file + ".tmp.seg.parasail.paf", args.output_file + ".tmp.seg.first.fa"], 
                    stdout = hout2, stderr = subprocess.DEVNULL)

            with open(args.output_file + ".tmp.seg.racon.fa", 'r') as hin2:
                header = hin2.readline()
                tconsensus = hin2.readline().rstrip('\n')
                print("@%s_%s_%s_%d_racon\n%s\n+\n%s" % 
                    (tchr, tstart, tend, tread_num, tconsensus, "I" * len(tconsensus)), file = hout)


    with open(args.output_file + ".tmp.polished.sam", 'w') as hout:
        subprocess.check_call(["bwa", "mem", args.reference_fasta, args.output_file + ".tmp.polished.fq"],
            stdout = hout, stderr = subprocess.DEVNULL) 

    func.summarize_bwa_alignment(args.output_file + ".tmp.polished.sam", args.output_file)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog = "local_polish_test")

    parser.add_argument("bam_file", default = None, type = str,
                        help = "Path to input BAM file")

    parser.add_argument("position_bed_file", default = None, type = str,
                        help = "Path to BED file to target region for polishing")

    parser.add_argument("reference_fasta", default = None, type = str,
                        help = "Path to the reference genome sequence")

    parser.add_argument("output_file", default = None, type = str,
                       help = "Path to output file")

    parser.add_argument("--read_num", type = int, default = None,
                        help = "Read number to use for polising")

    parser.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")

    args = parser.parse_args()

    main(args)
    
