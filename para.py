#! /usr/bin/env python

import os.path as op
from func import reverse_complement
import parasail


def read(sFile):
    """
    read a sequence file
    @param  sFile   sequence file
    """
    def read_one_fasta(f):
        """
        read a fasta file
        @param  f   file handler
        """
        sId = ''
        sSeq = ''
        for l in f:
            if l.startswith('>'):
                if sSeq:
                    yield sId, sSeq, ''
                sId = l.strip()[1:].split()[0]
                sSeq = ''
            else:
                sSeq += l.strip()

        yield sId, sSeq, ''

    def read_one_fastq(f):
        """
        read a fastq file
        @param  f   file handler
        """
        sId = ''
        sSeq = ''
        s3 = ''
        sQual = ''
        """
        for l in f:
            sId = l.strip()[1:].split()[0]
            sSeq = f.next()
            s3 = f.next()
            sQual = f.next()
            yield sId, sSeq, sQual
        """

        sID = f.readline()
        while sId:
            sId = sId.rstrip('\n').split()[0]
            sSeq = f.readline().rstrip('\n')
            sId2 = f.readline().rstrip('\n')
            sQual = f.readline().rstrip('\n')
            yield sId, sSeq, sQual
            sID = f.readline()


        f.close()

# test if fasta or fastq
    bFasta = True
    ext = op.splitext(sFile)[1][1:].strip().lower()
    if ext == 'gz' or ext == 'gzip':
        with gzip.open(sFile, 'r') as f:
            # l = f.next()
            l = f.readline().rstrip('\n')
            if l.startswith('>'):
                bFasta = True
            elif l.startswith('@'):
                bFasta = False
            else:
                print('file format cannot be recognized', file = sys.stderr)
                sys.exit()
    else:
        with open(sFile, 'r') as f:
            # l = f.next()
            l = f.readline().rstrip('\n')
            if l.startswith('>'):
                bFasta = True
            elif l.startswith('@'):
                bFasta = False
            else:
                print('file format cannot be recognized', file = sys.stderr)
                sys.exit()

# read
    if ext == 'gz' or ext == 'gzip':
        with gzip.open(sFile, 'r') as f:
            if bFasta == True:
                for sId,sSeq,sQual in read_one_fasta(f):
                    yield sId, sSeq, sQual
            else:
                for sId,sSeq,sQual in read_one_fastq(f):
                    yield sId, sSeq, sQual
    else:
        with open(sFile, 'r') as f:
            if bFasta == True:
                for sId,sSeq,sQual in read_one_fasta(f):
                    yield sId, sSeq, sQual
            else:
                for sId,sSeq,sQual in read_one_fastq(f):
                    yield sId, sSeq, sQual


def ssw_check_parasail(query, target):

    user_matrix = parasail.matrix_create("ACGT", 2, -2)

    alignment_info = {}
    for sQId, sQSeq, sQQual in read(query):
    
        sQSeq_r = reverse_complement(sQSeq)

        for sTId, sTSeq, STQual in read(target):

            res = parasail.ssw(sQSeq, sTSeq, 3, 1, user_matrix)
            res_r = parasail.ssw(sQSeq_r, sTSeq, 3, 1, user_matrix)

            if res.score1 > res_r.score1:
                score = res.score1
                qstart, qend = res.read_begin1 + 1., res.read_end1 + 1
                tstart, tend = res.ref_begin1 + 1, res.ref_end1 + 1
                strand = '+'
            else:
                score = res_r.score1
                qstart, qend = len(sQSeq) - res_r.read_end1, len(sQSeq) - res_r.read_begin1
                tstart, tend = res_r.ref_begin1 + 1, res_r.ref_end1 + 1
                strand = '-'

            alignment_info[sTId] = [score, int(qstart), int(qend), int(tstart), int(tend), strand]

    return(alignment_info)


if __name__ == "__main__":

    import sys
    import pdb; pdb.set_trace()
    print(ssw_check_parasail(sys.argv[1], sys.argv[2]))   
