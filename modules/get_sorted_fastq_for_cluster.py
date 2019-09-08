from __future__ import print_function
import os,sys
import argparse
import pysam

import signal
from multiprocessing import Pool
import multiprocessing as mp

import operator
import functools
import errno
from time import time
from collections import deque
import sys
import itertools
import math

from modules import help_functions

D = {chr(i) : min( 10**( - (ord(chr(i)) - 33)/10.0 ), 0.79433)  for i in range(128)}
D_no_min = {chr(i) : 10**( - (ord(chr(i)) - 33)/10.0 )  for i in range(128)}

def expected_number_of_erroneous_kmers(quality_string, k):
    prob_error = [D[char_] for char_ in quality_string]
    window = deque([ (1.0 - p_e) for p_e in prob_error[:k]])
    # print(window)
    qurrent_prob_no_error = functools.reduce(operator.mul, window, 1)
    # print(qurrent_prob_no_error)
    sum_of_expectations = qurrent_prob_no_error # initialization 
    for p_e in prob_error[k:]:
        p_to_leave = window.popleft()
        # print(window)
        # print(p_to_leave, "!" in quality_string)
        qurrent_prob_no_error *= ((1.0 -p_e)/(p_to_leave))
        # print(qurrent_prob_no_error)
        sum_of_expectations += qurrent_prob_no_error
        window.append(1.0 -p_e)
    return len(quality_string) - k + 1 - sum_of_expectations 


def reverse_complement(string):
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}
    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)

def batch(iterable, n=1):
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx:min(ndx + n, l)]


def calc_score_new(d):
    for key,value in d.items():
        l, k, q_threshold = value

    read_array = []
    error_rates = []
    for i, (acc, seq, qual) in enumerate(l):
        if i % 10000 == 0:
            print(i, "reads processed.")

        # skip very short reads or degenerate reads
        seq_hpol_comp = ''.join(ch for ch, _ in itertools.groupby(seq))
        if len(seq) < 2*k or len(seq_hpol_comp) < k:
            continue

        poisson_mean = sum([ qual.count(char_) * D_no_min[char_] for char_ in set(qual)])
        error_rate = poisson_mean/float(len(qual))
        if 10*-math.log(error_rate, 10) <= q_threshold:
            # print("Filtered read with:", 10*-math.log(error_rate, 10), error_rate)
            continue

        error_rates.append(error_rate)
        exp_errors_in_kmers = expected_number_of_erroneous_kmers(qual, k)
        p_no_error_in_kmers = 1.0 - exp_errors_in_kmers/ float((len(seq) - k +1))
        score =  p_no_error_in_kmers  * (len(seq) - k +1)
        read_array.append((acc, seq, qual, score) )
    return {key : (read_array, error_rates)}


def fastq_parallel(args):
    k = args.k
    q_threshold = args.quality_threshold
    error_rates = []
    reads = [ (acc,seq, qual) for acc, (seq, qual) in help_functions.readfq(open(args.fastq, 'r'))]
    start = time()
    read_chunk_size = int( len(reads)/args.nr_cores ) + 1
    read_batches = [b for b in batch(reads, read_chunk_size)]
    del reads
    ####### parallelize alignment #########
    # pool = Pool(processes=mp.cpu_count())
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.signal(signal.SIGINT, original_sigint_handler)
    mp.set_start_method('spawn')
    print(mp.get_context())
    print("Environment set:", mp.get_context())
    print("Using {0} cores.".format(args.nr_cores))
    start_multi = time()
    pool = Pool(processes=int(args.nr_cores))
    try:
        print([len(b) for b in read_batches])
        data = [ {i : (b,k, q_threshold)} for i, b in enumerate(read_batches)] #[ {i+1 :((cluster_batches[i], cluster_seq_origin_batches[i], read_batches[i], p_emp_probs, lowest_batch_index_db[i], i+1, args), {})} for i in range(len(read_batches))]
        res = pool.map_async(calc_score_new, data)
        score_results =res.get(999999999) # Without the timeout this blocking call ignores all signals.
    except KeyboardInterrupt:
        print("Caught KeyboardInterrupt, terminating workers")
        pool.terminate()
        sys.exit()
    else:
        pool.close()
    pool.join()

    print("Time elapesd multiprocessing:", time() - start_multi)
    read_array, error_rates = [], []

    for output_dict in score_results:
        for k, v in output_dict.items():
            r_a, err_rates = v
            print("Batch index", k)
            for item in r_a:
                read_array.append(item)
            for item2 in err_rates:
                error_rates.append(item2)

    read_array.sort(key=lambda x: x[3], reverse=True)
    error_rates.sort()
    return read_array, error_rates


def fastq_single_core(args):
    k = args.k
    q_threshold = args.quality_threshold
    error_rates = []
    read_array = []
    for i, (acc, (seq, qual)) in enumerate(help_functions.readfq(open(args.fastq, 'r'))):
        if i % 10000 == 0:
            print(i, "reads processed.")

        # skip very short reads or degenerate reads
        seq_hpol_comp = ''.join(ch for ch, _ in itertools.groupby(seq))
        if len(seq) < 2*k or len(seq_hpol_comp) < args.k:
            continue
        ########################
    
        exp_errors_in_kmers = expected_number_of_erroneous_kmers(qual, k)
        p_no_error_in_kmers = 1.0 - exp_errors_in_kmers/ float((len(seq) - k +1))
        score =  p_no_error_in_kmers  * (len(seq) - k +1)
        
        ## For (inferred) average error rate only, based on quality values
        ### These values are used in evaluations in the paper only, and are not used in clustering
        poisson_mean = sum([ qual.count(char_) * D_no_min[char_] for char_ in set(qual)])
        error_rate = poisson_mean/float(len(qual))
        if 10*-math.log(error_rate, 10) <= q_threshold:
            # print("Filtered read with:", 10*-math.log(error_rate, 10), error_rate)
            continue
        error_rates.append(error_rate)
        ##############################################
        
        read_array.append((acc, seq, qual, score) )

    read_array.sort(key=lambda x: x[3], reverse=True)
    return read_array, error_rates


def isoseq(args):
    k = args.k
    error_rates = []
    flnc_file = pysam.AlignmentFile(args.flnc, "rb", check_sq=False)
    ccs_file = pysam.AlignmentFile(args.ccs, "rb", check_sq=False)
    flnc_dict = {}
    for read in flnc_file.fetch(until_eof=True):
        
        # while quality values are not implemented in unpolished.flnc
        flnc_dict[read.qname] = read.seq
        # If quality values gets implemented, use this one-liner insetead..
        # flnc_dict[read.qname] = (read.seq, read.qual)
    
    read_array = []
    for read in ccs_file.fetch(until_eof=True):
        if read.qname in flnc_dict:
            # while quality values are not implemented in unpolished.flnc
            seq = flnc_dict[read.qname]
            full_seq = read.seq
            full_seq_rc = reverse_complement(full_seq)
            
            if seq in full_seq:
                start_index = full_seq.index(seq)
                stop_index = start_index + len(seq)
                qual = read.qual[start_index: stop_index]

            elif seq in full_seq_rc:
                qual = read.qual[::-1]
                start_index = full_seq_rc.index(seq)
                stop_index = start_index + len(seq)
                qual = qual[start_index: stop_index]

            else:
                print("Bug, flnc not in ccs file")
                sys.exit()

            assert len(qual) == len(seq)

            poisson_mean = sum([ qual.count(char_) * D_no_min[char_] for char_ in set(qual)])
            error_rate = poisson_mean/float(len(qual))
            error_rates.append(error_rate)
            exp_errors_in_kmers = expected_number_of_erroneous_kmers(qual, k)
            p_no_error_in_kmers = 1.0 - exp_errors_in_kmers/ float((len(seq) - k +1))
            score =  p_no_error_in_kmers  * (len(seq) - k +1)

            read_array.append((read.qname, seq, qual, score) )

            # If quality values gets implemented, simply use the code below and remove everythin above..
            # seq, qual = flnc_dict[read.qname][0], flnc_dict[read.qname][1]
            # p_no_error_in_kmers_appr =  get_p_no_error_in_kmers_approximate(qual,k)
            # score = p_no_error_in_kmers_appr * len(seq)
            # read_array.append((read.qname, seq, qual, score) )
    read_array.sort(key=lambda x: x[3], reverse=True)
    return read_array, error_rates


def main(args):
    start = time()
    logfile = open(os.path.join(args.outfolder, "logfile.txt"), 'w')
    if os.path.isfile(args.outfile) and args.use_old_sorted_file:
        print("Using already existing sorted file in specified directory, in not intended, specify different outfolder or delete the current file.")
        return args.outfile

    elif args.fastq:
        if args.nr_cores > 1: 
            read_array, error_rates = fastq_parallel(args)
        else:
            read_array, error_rates = fastq_single_core(args)

    elif args.flnc and args.ccs:
        read_array, error_rates = isoseq(args)
    else:
        print("Wrong input format")
        sys.exit()


    reads_sorted_outfile = open(args.outfile, "w")
    for i, (acc, seq, qual, score) in enumerate(read_array):
        reads_sorted_outfile.write("@{0}\n{1}\n+\n{2}\n".format(acc + "_{0}".format(score), seq, qual))
    reads_sorted_outfile.close()
    print(len(read_array), "reads passed quality critera (avg phred Q val over {0} and length > 2*k) and will be clustered.".format(args.quality_threshold))
    error_rates.sort()
    min_e = error_rates[0]
    max_e = error_rates[-1]
    median_e = error_rates[int(len(error_rates)/2)]
    mean_e = sum(error_rates)/len(error_rates)
    logfile.write("Lowest read error rate:{0}\n".format(min_e))
    logfile.write("Highest read error rate:{0}\n".format(max_e))
    logfile.write("Median read error rate:{0}\n".format(median_e))
    logfile.write("Mean read error rate:{0}\n".format(mean_e))
    logfile.write("\n")
    logfile.close()
    print("Sorted all reads in {0} seconds.".format(time() - start) )
    return reads_sorted_outfile.name


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('--fastq', type=str,  default=False, help='Path to consensus fastq file(s)')
    parser.add_argument('--flnc', type=str, default=False, help='The flnc reads generated by the isoseq3 algorithm (BAM file)')
    parser.add_argument('--ccs', type=str, default=False, help='Path to lima demultiplexed BAM file')
    parser.add_argument('--outfile', type=str,  default=None, help='A fasta file with transcripts that are shared between samples and have perfect illumina support.')
    parser.add_argument('--k', type=int, default=15, help='kmer size')
    
    args = parser.parse_args()

    if (args.fastq and (args.flnc or args.ccs)):
        print("Either (1) only a fastq file, or (2) a ccs and a flnc file should be specified. ")
        sys.exit()

    if (args.flnc != False and args.ccs == False ) or (args.flnc == False and args.ccs != False ):
        print("qt-clust needs both the ccs.bam file produced by ccs and the flnc file produced by isoseq3 cluster. ")
        sys.exit()


    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()
    path_, file_prefix = os.path.split(args.outfile)
    help_functions.mkdir_p(path_)

    main(args)