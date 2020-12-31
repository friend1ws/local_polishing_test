#! /usr/bin/env python

from collections import Counter
import pysam
import parasail

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'W': 'W', 'S': 'S', 'M': 'K', 'K': 'M', 'R': 'Y', 'Y': 'R',
                  'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N'}

    return("".join(complement.get(base, base) for base in reversed(seq)))


def get_read_segments(input_bam, output_file, tchr, tstart, tend, tread_num = None, mapq_thres = 40):

    read_num = 0
    with pysam.AlignmentFile(input_bam, 'rb') as bam_h, open(output_file, 'w') as hout:

        for read in bam_h.fetch(tchr, tstart, tend):

            if read.is_secondary or read.is_supplementary: continue
            if read.mapping_quality < mapq_thres: continue        

            query_name = read.query_name
            query_strand = '-' if read.is_reverse else '+'
            query_length = read.infer_read_length()

            reference_name = read.reference_name
            reference_start = read.reference_start + 1
            reference_end = read.reference_end

            cigartuples = read.cigartuples
            left_hard_clipping_size, right_hard_clipping_size = 0, 0
            left_soft_clipping_size, right_soft_clipping_size = 0, 0
            if cigartuples[0][0] == 5: left_hard_clipping_size = cigartuples[0][1]
            if cigartuples[-1][0] == 5: right_hard_clipping_size = cigartuples[-1][1]
            if cigartuples[0][0] == 4: left_soft_clipping_size = cigartuples[0][1]
            if cigartuples[-1][0] == 4: right_soft_clipping_size = cigartuples[-1][1]

            """
            if query_strand == '+':
                query_start = read.query_alignment_start + 1 
                query_end = read.query_alignment_end
            else:
                query_start = query_length - read.query_alignment_end + 1
                query_end = query_length - read.query_alignment_start

            query_pos_cur = query_start - 1 if query_strand == '+' else query_end
            query_pos_check = query_pos_cur
            reference_pos_cur = reference_start - 1
            reference_pos_check = reference_start - 1
            sstart, send = None, None
            """

            query_pos_cur = read.query_alignment_start
            reference_pos_cur = reference_start - 1 
            reference_pos_check = reference_start - 1
            sstart, send = None, None
            # if query_name == "b6eec70e-7769-4214-a2b6-a5742bb865a8":
            #     import pdb; pdb.set_trace()

            for cigar in cigartuples:
                if cigar[0] == 0:

                    if tstart >= reference_pos_cur and tstart < reference_pos_cur + cigar[1]:
                        sstart = query_pos_cur + (tstart - reference_pos_cur)
                    if tend >= reference_pos_cur and tend < reference_pos_cur + cigar[1]:
                        send = query_pos_cur + (tend - reference_pos_cur)

                    query_pos_cur = query_pos_cur + cigar[1] # if query_strand == '+' else query_pos_cur - cigar[1]
                    reference_pos_cur = reference_pos_cur + cigar[1] 

                elif cigar[0] == 1:
                    query_pos_cur = query_pos_cur + cigar[1] # if query_strand == '+' else query_pos_cur - cigar[1]
                elif cigar[0] == 2: 

                    if tstart >= reference_pos_cur and tstart < reference_pos_cur + cigar[1]:
                        sstart = query_pos_cur
                    if tend >= reference_pos_cur and tend < reference_pos_cur + cigar[1]:
                        send = query_pos_cur + cigar[1]

                    reference_pos_cur = reference_pos_cur + cigar[1]


            if sstart is None or send is None: continue
            print(">%s\n%s" % (query_name, read.query_sequence[sstart:send]), file = hout)
            read_num = read_num + 1
            if tread_num is not None and tread_num == read_num: break

    return(read_num)

def get_consensus_from_mafft_result(input_file):

    id2seq = {}
    with open(input_file, 'r') as hin:
        for line in hin:
            line = line.rstrip('\n')
            if line.startswith('>'):
                tid = line
                id2seq[tid] = ''
            else:
                id2seq[tid] = id2seq[tid] + line


    ind2bases = {}
    for tid in id2seq:
        seq = id2seq[tid]
        for i in range(len(seq)):
            if i not in ind2bases: ind2bases[i] = []
            ind2bases[i].append(seq[i])


    seq_len = len(list(ind2bases))
    consensus = ''
    for i in range(seq_len):
        # import pdb; pdb.set_trace()
        mycounter = Counter(ind2bases[i] )
        consensus = consensus + mycounter.most_common()[0][0]

    consensus = consensus.replace('-', '').upper()

    return(consensus)


def summarize_bwa_alignment(input_sam, output_file):

    with pysam.AlignmentFile(input_sam, 'r') as sam_h, open(output_file, 'w') as hout:
        for read in sam_h.fetch():

            if read.is_unmapped: continue
            if read.is_secondary: continue

            query_name = read.query_name
            key = query_name

            query_strand = '-' if read.is_reverse else '+'
            reference_name = read.reference_name
            reference_start = read.reference_start + 1
            reference_end = read.reference_end
            mapping_quality = read.mapping_quality
            query_length = read.infer_read_length()

            cigartuples = read.cigartuples
            left_hard_clipping_size, right_hard_clipping_size = 0, 0
            total_match_size, total_ins_size, total_del_size = 0, 0, 0
            if cigartuples[0][0] == 5: left_hard_clipping_size = cigartuples[0][1]
            if cigartuples[-1][0] == 5: right_hard_clipping_size = cigartuples[-1][1]

            for tcigar in cigartuples:
                if tcigar[0] == 0: total_match_size = total_match_size + tcigar[1]
                if tcigar[0] == 1: total_ins_size = total_ins_size + tcigar[1]
                if tcigar[0] == 2: total_del_size = total_del_size + tcigar[1]

            total_mismatch_size = 0
            for ttag in read.get_tags():
                if ttag[0] == "NM": total_mismatch_size = total_mismatch_size + ttag[1]

            total_match_size = total_match_size - total_mismatch_size

            error_ratio = float(total_mismatch_size + total_ins_size + total_del_size) / \
                float(total_match_size + total_mismatch_size + total_ins_size + total_del_size)
 
            print("%s\t%d\t%d\t%d\t%d\t%f" % 
                (key, total_match_size, total_mismatch_size, total_ins_size, total_del_size, error_ratio),
                file = hout)



def generate_paf_file(query_fasta, target_fasta, output_file):

    user_matrix = parasail.matrix_create("ACGT", 2, -2)

    with open(target_fasta, 'r') as hin:
        for line in hin:
            if line.startswith('>'): 
                tid = line.rstrip('\n').lstrip('>')
            else:
                tseq = line.rstrip('\n')

    with open(query_fasta, 'r') as hin, open(output_file, 'w') as hout:
        for line in hin:
            if line.startswith('>'):
                qid = line.rstrip('\n').lstrip('>')
            else:
                qseq = line.rstrip('\n')
                
                res = parasail.ssw(qseq, tseq, 3, 1, user_matrix)
                print("%s\t%d\t%d\t%d\t+\t%s\t%d\t%d\t%d\t*\t*\t60" %
                    (qid, len(qseq), res.read_begin1, res.read_end1,
                    tid, len(tseq), res.ref_begin1, res.ref_end1), file = hout)
       
                """ 
                print("%s\t%d\t%d\t%d\t+\t%s\t%d\t%d\t%d\t*\t*\t60" %
                    (qid, len(qseq), 0, len(qseq),
                    tid, len(tseq), 0, len(tseq)), file = hout)
                """


if __name__ == "__main__":

    import sys
    print(get_consensus_from_mafft_result(sys.argv[1]))

