from __future__ import print_function
from functools import reduce
import os,sys
import argparse

import pickle
# import pysam
from collections import defaultdict, Counter
import math
from collections import deque
import itertools
from operator import mul

import signal
from multiprocessing import Pool
import multiprocessing as mp

import errno
from time import time
import re
import copy

try:
    import parasail
except:
    pass



'''
    Below code taken from https://github.com/lh3/readfq/blob/master/readfq.py
'''

def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].replace(" ", "_"), [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, (''.join(seqs), None) # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, (seq, ''.join(seqs)); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, (seq, None) # yield a fasta record instead
                break


def cigar_to_seq(cigar, query, ref):
    cigar_tuples = []
    result = re.split(r'[=DXSMI]+', cigar)
    i = 0
    for length in result[:-1]:
        i += len(length)
        type_ = cigar[i]
        i += 1
        cigar_tuples.append((int(length), type_ ))

    r_index = 0
    q_index = 0
    q_aln = []
    r_aln = []
    for length_ , type_ in cigar_tuples:
        if type_ == "=" or type_ == "X":
            q_aln.append(query[q_index : q_index + length_])
            r_aln.append(ref[r_index : r_index + length_])

            r_index += length_
            q_index += length_
        
        elif  type_ == "I":
            # insertion w.r.t. reference
            r_aln.append('-' * length_)
            q_aln.append(query[q_index: q_index + length_])
            #  only query index change
            q_index += length_

        elif type_ == 'D':
            # deletion w.r.t. reference
            r_aln.append(ref[r_index: r_index + length_])
            q_aln.append('-' * length_)
            #  only ref index change
            r_index += length_
        
        else:
            print("error")
            print(cigar)
            sys.exit()

    return  "".join([s for s in q_aln]), "".join([s for s in r_aln])


def homopolymer_compress(seq):
    corr = [ n1  for n1,n2 in zip(seq[:-1], seq[1:]) if n1 != n2 ]
    #last base corner case
    if seq[-1] != seq[-2]:
        corr.append(seq[-1])
    return "".join([nt for nt in corr])


def get_kmer_minimizers(seq, k_size, w_size):
    # kmers = [seq[i:i+k_size] for i in range(len(seq)-k_size) ]
    w = w_size - k_size
    window_kmers = deque([seq[i:i+k_size] for i in range(w)])
    curr_min = min(window_kmers)
    minimizers = [ (curr_min, list(window_kmers).index(curr_min)) ]

    for i in range(w,len(seq) - k_size):
        new_kmer = seq[i:i+k_size]
        # updateing window
        discarded_kmer = window_kmers.popleft()
        window_kmers.append(new_kmer)

        # we have discarded previous windows minimizer, look for new minimizer brute force
        if curr_min == discarded_kmer: 
            curr_min = min(window_kmers)
            minimizers.append( (curr_min, list(window_kmers).index(curr_min) + i - w ) )

        # Previous minimizer still in window, we only need to compare with the recently added kmer 
        elif new_kmer < curr_min:
            curr_min = new_kmer
            minimizers.append( (curr_min, i) )

    return minimizers



def get_all_hits_new(minimizers, H, Clusters, read_cl_id):
    hit_clusters_ids = defaultdict(int)
    hit_clusters_hit_index = defaultdict(list)
    hit_clusters_hit_positions = defaultdict(list)
    for i, (m, pos) in enumerate(minimizers): # iterating over minimizers from upstream to downstream in read
        if m in H:
            for cl_id in H[m]: 
                hit_clusters_ids[cl_id] += 1
                hit_clusters_hit_index[cl_id].append(i) # index of the minimizer among the coordinate sorted minimizers in the read
                hit_clusters_hit_positions[cl_id].append(pos) # positions of the minimizer among the coordinate sorted minimizers in the read

    if read_cl_id in hit_clusters_ids:
        del hit_clusters_ids[read_cl_id]
        del hit_clusters_hit_index[read_cl_id]
        del hit_clusters_hit_positions[read_cl_id]

    return hit_clusters_ids, hit_clusters_hit_index, hit_clusters_hit_positions




def get_best_cluster(read_cl_id, compressed_seq_len, hit_clusters_ids, hit_clusters_hit_positions, minimizers, nummber_of_minimizers, hit_clusters_hit_index, cluster_seq_origin, p_emp_probs, args):
    """
        Tally up total covered (mapped region) and compare it with total unmapped region. What is counted as a consecutive mapped block depends on minimizer qualitity
        threshold as mapped-read-legnth/total-read-length is chosen to classify a read as belonging to the cluster
        
        Return: An integer >= 0 denoting the cluster ID that this read was assigned to. In not assigend to any previous cluster, return -1. 
                [Also returing mapped ratio and nr shared minimizers to best read for logging purposes.]
    """
    best_cluster_id = -1
    nr_shared_kmers = 0
    mapped_ratio = 0.0
    if hit_clusters_ids:
        top_matches = sorted(hit_clusters_ids.items(), key=lambda x: x[1],  reverse=True)
        top_hits = top_matches[0][1]
        nr_shared_kmers = top_hits
        if top_hits < args.min_shared:
            pass
        else:
            for tm in top_matches:
                cl_id = tm[0]
                nm_hits = tm[1]
                if nm_hits < args.min_fraction * top_hits or nm_hits < args.min_shared:
                    break

                cl_size = len(hit_clusters_ids)
                minimizer_hit_positions = hit_clusters_hit_positions[cl_id]
                minimizer_hit_indices = hit_clusters_hit_index[cl_id]
                assert len(minimizer_hit_indices) == len(minimizer_hit_positions)
                _, _, _, _, _, error_rate_c = cluster_seq_origin[cl_id]
                _, _, _, _, _, error_rate_read = cluster_seq_origin[read_cl_id]
                p_error_in_kmers_emp =  1.0 - p_shared_minimizer_empirical(error_rate_read, error_rate_c, p_emp_probs)
                minimizer_error_probabilities = [p_error_in_kmers_emp]*nummber_of_minimizers
                total_mapped = 0
                prev_mpos = 0
                prob_all_errors_since_last_hit = [reduce(mul, minimizer_error_probabilities[: minimizer_hit_indices[0]], 1)] +  [ reduce(mul, minimizer_error_probabilities[hit_idx1+1: hit_idx2], 1) for hit_idx1, hit_idx2 in zip(minimizer_hit_indices[:-1], minimizer_hit_indices[1:]) ] + [reduce(mul, minimizer_error_probabilities[minimizer_hit_indices[-1]+1 : ], 1)]

                assert len(prob_all_errors_since_last_hit) == len(minimizer_hit_positions) + 1
                for i in range(len(minimizer_hit_indices)):
                    if prob_all_errors_since_last_hit[i] < args.min_prob_no_hits:
                        pass
                    else:
                        if i == 0:
                            total_mapped += minimizer_hit_positions[i]
                        else:
                            total_mapped += minimizer_hit_positions[i] - minimizer_hit_positions[i-1]
                if prob_all_errors_since_last_hit[-1] < args.min_prob_no_hits:
                    pass
                else:
                    total_mapped += compressed_seq_len - minimizer_hit_positions[-1]

                mapped_ratio = total_mapped /float(compressed_seq_len) 
                
                if mapped_ratio > args.mapped_threshold:
                    is_covered = True
                    return cl_id, nm_hits, mapped_ratio

    return  best_cluster_id, nr_shared_kmers, mapped_ratio 


def parasail_block_alignment(s1, s2, k, match_id, x_acc = "", y_acc = "", match_score = 2, mismatch_penalty = -2, opening_penalty = 5, gap_ext = 1, ends_discrepancy_threshold = 0):
    user_matrix = parasail.matrix_create("ACGT", match_score, mismatch_penalty)
    result = parasail.sg_trace_scan_16(s1, s2, opening_penalty, gap_ext, user_matrix)
    if result.saturated:
        print("SATURATED!")
        result = parasail.sg_trace_scan_32(s1, s2, opening_penalty, gap_ext, user_matrix)
    if sys.version_info[0] < 3:
        cigar_string = str(result.cigar.decode).decode('utf-8')
    else:
        cigar_string = str(result.cigar.decode, 'utf-8')
    
    s1_alignment, s2_alignment = cigar_to_seq(cigar_string, s1, s2)

    # Rolling window of matching blocks
    # k=15
    # match_id = int(k*0.8)  1.0 - math.ceil(window_fraction)
    match_vector = [ 1 if n1 == n2 else 0 for n1, n2 in zip(s1_alignment, s2_alignment) ]
    # print("".join([str(m) for m in match_vector]))
    
    match_window = deque(match_vector[:k]) # initialization
    current_match_count = sum(match_window)
    aligned_region = []
    if current_match_count >= match_id:
        aligned_region.append(1)
    else:
        aligned_region.append(0)


    for new_m_state in match_vector[k:]:
        prev_m_state = match_window.popleft()
        current_match_count = current_match_count - prev_m_state + new_m_state 
        match_window.append(new_m_state)
        
        if current_match_count >= match_id:
            aligned_region.append(1)
        else:        
            aligned_region.append(0)

    # print("".join([str(m) for m in aligned_region]))
    # print("Aligned ratio (tot aligned/len(seq1):", sum(aligned_region)/float(len(s1)))
    alignment_ratio = sum(aligned_region)/float(len(s1))
    return (s1, s2, (s1_alignment, s2_alignment, alignment_ratio))


def get_best_cluster_block_align(read_cl_id, cluster_seq_origin, hit_clusters_ids, phred_char_to_p, args):
    best_cluster_id = -1
    top_matches = sorted(hit_clusters_ids.items(), key=lambda x: x[1],  reverse=True)
    _, _, seq, r_qual, _, _ = cluster_seq_origin[read_cl_id]
    # print(top_matches)
    top_hits = top_matches[0][1]
    aln_counter = 0
    alignment_ratio = 0.0
    for tm in top_matches:
        cl_id = tm[0]
        nm_hits = tm[1]
        if nm_hits < top_hits:
            break
        aln_counter +=1
        _, _, c_seq, c_qual, _, _ = cluster_seq_origin[cl_id]

        poisson_mean = sum([ r_qual.count(char_) * phred_char_to_p[char_] for char_ in set(r_qual)])
        poisson_mean2 = sum([ c_qual.count(char_) * phred_char_to_p[char_] for char_ in set(c_qual)])

        error_rate_sum = poisson_mean/float(len(seq)) + poisson_mean2/float(len(c_seq))  # k = max(int(mean_plus_two_stdvs_q2 + mean_plus_two_stdvs_q1) + 1 + int(len(seq)*args.variant_rate) , 40)
        if error_rate_sum <= 0.01:
            gap_opening_penalty = 5
        elif  0.01 < error_rate_sum <= 0.04:
            gap_opening_penalty = 4
        elif  0.04 < error_rate_sum <= 0.1:
            gap_opening_penalty = 3
        elif  0.1 < error_rate_sum:
            gap_opening_penalty = 2

        match_id_tailored = math.floor((1.0 - error_rate_sum) * args.k)
        (s1, s2, (s1_alignment, s2_alignment, alignment_ratio)) = parasail_block_alignment(seq, c_seq, args.k, match_id_tailored, opening_penalty = gap_opening_penalty,  )
        # print("Expected errors:", poisson_mean, poisson_mean2)
        if alignment_ratio >= args.aligned_threshold: #args.mapped_threshold:
            return cl_id, nm_hits,  error_rate_sum, s1_alignment, s2_alignment, alignment_ratio

    return  best_cluster_id, 0,  -1, -1, -1, alignment_ratio


def reads_to_clusters(Cluster, cluster_seq_origin, sorted_reads, p_emp_probs, args, minimizer_database = {}, task = "clustering", new_cluster_allowed = True):
    """
        Iterates throughreads in sorted order (w.r.t. score) and:
            1. homopolymenr compresses the read
            2. Finds the homopolymenr compressed error rate (if not computed in previous pass if more than 1 core specified to the program)
            3. Finds all the representatives with shared minimizers (in "get_all_hits_new")
            4. Finds the best of the hits using mapping approach
            5. If no hit is found in 4. tries to align to representative with th most shared minimizers.
            6. Adds current read to representative, or makes it a new representative of a new cluster.
                6''. If new representative, and add the minimizers to the miniizer database. 
    """
    print("USING w:{0}, k:{1}".format(args.w, args.k))


    if task  == "clustering":
        init_clusters = copy.deepcopy(Cluster)
        init_cluster_seq_origin = copy.deepcopy(cluster_seq_origin)

    phred_char_to_p = {chr(i) : min( 10**( - (ord(chr(i)) - 33)/10.0 ), 0.5)  for i in range(128)} # PHRED encoded quality character to prob of error. Need this locally if multiprocessing
    cluster_to_new_cluster_id = {}
    ## logging counters 
    aln_passed_criteria = 0
    mapped_passed_criteria = 0
    aln_called = 0
    # total_reads = 0
    ###################

    for i, (read_cl_id, acc, seq, qual, score) in enumerate(sorted_reads):
        
        ################################################################################
        ############  Just for develop purposes, print some info to std out ############
        if i%5000 == 0 and i > 0: 
            inv_map = {}
            for k, v in cluster_to_new_cluster_id.items():
                inv_map.setdefault(v, set()).add(k)
            cl_tmp = sorted( [ 1 + sum([len(Cluster[cl_id]) for cl_id in c ]) for c in inv_map.values() ], reverse= True)
            cl_tmp_nontrivial = [cl_size_tmp for cl_size_tmp in cl_tmp if cl_size_tmp > 1]
            print("Processing read", i, "seq length:", len(seq), "nr non-trivial clusters:", len(cl_tmp_nontrivial), "kmers stored:", len(minimizer_database))
            print("clust distr:", [c_len for c_len in cl_tmp if c_len > 100] )
            depth = [len(nr_cl) for kmer, nr_cl in  sorted(minimizer_database.items(), key=lambda x: len(x[1]), reverse= True)[:50]]
            print("Depth of minimizer_database:", sum(depth)/float(len(depth)), depth)
        ################################################################################
        ################################################################################

        if task  == "read_assignment":
            if read_cl_id in cluster_seq_origin: # last step read assignment
                # print("HERE")
                Cluster[read_cl_id].append(acc)
                continue


        # homopolymenr compress read
        seq_hpol_comp = ''.join(ch for ch, _ in itertools.groupby(seq))
        assert read_cl_id not in cluster_seq_origin

        indices = [i for i, (n1,n2) in enumerate(zip(seq[:-1],seq[1:])) if n1 != n2] # indicies we want to take quality values from to get quality string of homopolymer compressed read 
        indices.append(len(seq) - 1)
        qualcomp = ''.join([qual[i] for i in indices])
        assert len(seq_hpol_comp) == len(qualcomp)

        # compute the average error rate after compression
        poisson_mean = sum([ qualcomp.count(char_) * phred_char_to_p[char_] for char_ in set(qualcomp)])
        h_pol_compr_error_rate = poisson_mean/float(len(qualcomp))
        cluster_seq_origin[read_cl_id] = (read_cl_id, acc, seq, qual, score, h_pol_compr_error_rate) # adding homopolymenr compressed error rate to info tuple of cluster origin sequence
        
        # get minimizers
        if len(seq_hpol_comp) < args.k:
            print( "skipping read of length:", len(seq), "homopolymer compressed:", len(seq_hpol_comp), seq)
            continue


        minimizers = get_kmer_minimizers(seq_hpol_comp, args.k, args.w)
        hit_clusters_ids, hit_clusters_hit_index, hit_clusters_hit_positions = get_all_hits_new(minimizers, minimizer_database, Cluster, read_cl_id)
        best_cluster_id_m, nr_shared_kmers_m, mapped_ratio = get_best_cluster(read_cl_id, len(seq_hpol_comp), hit_clusters_ids, hit_clusters_hit_positions, minimizers, len(minimizers), hit_clusters_hit_index, cluster_seq_origin, p_emp_probs, args)
        
        if best_cluster_id_m >= 0:
            mapped_passed_criteria += 1

        if best_cluster_id_m < 0 and nr_shared_kmers_m >= args.min_shared:
            aln_called += 1
            best_cluster_id_a, nr_shared_kmers_a, error_rate_sum, s1_alignment, s2_alignment, alignment_ratio = get_best_cluster_block_align(read_cl_id, cluster_seq_origin, hit_clusters_ids, phred_char_to_p, args)
            if best_cluster_id_a >= 0:
                aln_passed_criteria += 1

        else:
            best_cluster_id_a = -1


        best_cluster_id = max(best_cluster_id_m,best_cluster_id_a)
        if read_cl_id == 5432:
            print(5432, best_cluster_id, acc, best_cluster_id_a, nr_shared_kmers_a, best_cluster_id_m)
        
        if best_cluster_id >= 0:
            # cluster_to_new_cluster_id[read_cl_id] = best_cluster_id
            Cluster[best_cluster_id].append(acc) 
            del cluster_seq_origin[read_cl_id]

        elif new_cluster_allowed:  # New cluster, adding representative minimixers
            # if task == "read_assignment":
            #     print("THIS SHOULD NOT HAPPEN")
            #     print(acc, seq)
                # sys.exit()
            Cluster[read_cl_id] = [acc]
            for m, pos in minimizers:
                if m in minimizer_database:
                    minimizer_database[m].add(read_cl_id)
                else:
                    minimizer_database[m] = set()
                    minimizer_database[m].add(read_cl_id)     

    
    # new no graph approach ####
    # for read_cl_id in cluster_to_new_cluster_id:
    #     new_cl_id = cluster_to_new_cluster_id[read_cl_id]
    #     # add all read acc to new origin
    #     all_reads = Cluster[read_cl_id]
    #     for read_acc in all_reads:
    #         Cluster[new_cl_id].append(read_acc)
    #     del Cluster[read_cl_id]
    #     # delete old origins
    #     del cluster_seq_origin[read_cl_id]
    ##########################


    # print("PASS")
    # print("Total number of reads iterated through:{0}".format(i+1))

    # print("Passed mapping criteria:{0}".format(mapped_passed_criteria))
    # print("Passed alignment criteria in this process:{0}".format(aln_passed_criteria))
    # print("Total calls to alignment mudule in this process:{0}".format(aln_called))

    # print("Percent passed mapping criteria:{0}".format( round(100*mapped_passed_criteria/float(i+1), 2) ))
    # print("Percent passed alignment criteria total:{0}".format( round(100*aln_passed_criteria/float(i+1), 2) ))    
    # if aln_called > 0:
    #     print("Percent passed alignment criteria out of number of calls to the alignment module:{0}".format(round(100*aln_passed_criteria/float(aln_called), 2) )) 

    if task  == "read_assignment":
        return Cluster, cluster_seq_origin, minimizer_database

    else:
        new_cluster_seq_origins = {k : v for k, v in cluster_seq_origin.items() if k not in init_cluster_seq_origin}
        new_clusters = {k : v for k,v in Cluster.items() if k not in init_clusters}
        if new_cluster_allowed:
            assert len(new_clusters) == len(new_cluster_seq_origins)
        else:
            assert len(new_clusters) == 0

        return new_clusters, new_cluster_seq_origins, minimizer_database


def p_shared_minimizer_empirical(error_rate_read, error_rate_center, p_emp_probs):
    e1 = round(error_rate_read, 2)
    # print(e1)
    if e1 > 0.15:
        e1 = 0.15
    if e1 < 0.01:
        e1 = 0.01
    e2 = round(error_rate_center, 2)
    # print(e2)
    if e2 > 0.15:
        e2 = 0.15
    if e2 < 0.01:
        e2 = 0.01
    # print(e2)
    p_kmer_shared = p_emp_probs[(e1,e2)]
    return p_kmer_shared

def reads_to_clusters_helper(arguments):
    args, kwargs = arguments
    return reads_to_clusters(*args, **kwargs)

def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z

def merge_dicts(*dict_args):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result

def batch(lst, n=1):
    l = len(lst)
    for ndx in range(0, l, n):
        yield lst[ndx:min(ndx + n, l)]


def parallelize(function, data_chunked, nr_cores):
    pool = Pool(processes=nr_cores-1)
    try:
        res = pool.map_async(function, data_chunked )
        results = res.get(999999999) # Without the timeout this blocking call ignores all signals.
    except KeyboardInterrupt:
        print("Caught KeyboardInterrupt, terminating workers")
        pool.terminate()
        sys.exit()
    else:
        # print("Normal termination")
        pool.close()
    pool.join()
    return results


def cluster_seqs(read_array, p_emp_probs, args):
    all_start = time()
    # split sorted reads into batches
    print("Using {0} batches (#cores - 1).".format(args.nr_cores -1))

    # process the first batch single core
    
    # top_read_batch = read_array[0:10000] # TODO: make 10k a parameter with default = 100,000 
    # top_new_clusters, top_new_seq_origins, top_minimizer_database = reads_to_clusters({}, {}, top_read_batch, p_emp_probs, args, minimizer_database = {})
    
    # Parallelize over args.nr_cores cores
    chunk_size = int((len(read_array) )/ (args.nr_cores - 1)) + 1
    print("Chunk sizes:", chunk_size)
    # print("Number of reference clusters:", len(top_new_clusters))
    # cluster_batches = []
    # cluster_seq_origin_batches = []
    read_batches = [batch for batch in batch(read_array, chunk_size)]

    # print( [info for info in read_array if "m141129_125831_42161_c100698142550000001823143403261592_s1_p0/27236/30_955_CCS_strand=+;fiveseen=1;polyAseen=1;threeseen=1;fiveend=30;polyAend=955;threeend=981;primer=1;chimera=0_900.8398439504809" == info[1] ])
    # if not read_batches:
    #     return top_new_clusters, top_new_seq_origins

    # print(read_batches)
    # for sublist in batch(read_array[3000:], chunk_size):
    #     tmp_clust = {}
    #     tmp_clust_origin = {}
    #     for i, acc, seq, qual, score in sublist:
    #         tmp_clust[i] = [acc]
    #         tmp_clust_origin[i] = (i, acc, seq, qual, score)
        # cluster_batches.append(tmp_clust)
        # cluster_seq_origin_batches.append(tmp_clust_origin)
    # results = parallelize(reads_to_clusters_helper, [((top_new_clusters, top_new_seq_origins, read_batches[i], p_emp_probs, args), {"minimizer_database" : top_minimizer_database}) for i in range(len(read_batches))], args.nr_cores)

    results = parallelize(reads_to_clusters_helper, [(({}, {}, read_batches[i], p_emp_probs, args), {"minimizer_database" : {}}) for i in range(len(read_batches))], args.nr_cores)
    # all_cl, all_repr, all_minimizer_databases = [top_new_clusters],[top_new_seq_origins], [minimizer_database]
    all_repr = [] # all_repr = [top_new_seq_origins]
    all_cl = []
    for new_clusters, new_cluster_origins, min_db in results: 
        all_cl.append(new_clusters)
        all_repr.append(new_cluster_origins)

    all_representatives = merge_dicts(*all_repr)
    all_representatives_sorted = sorted(list(all_representatives.values()), key = lambda x: x[4], reverse = True)
    all_representatives_sorted = [ (read_cl_id, acc, seq, qual, score) for (read_cl_id, acc, seq, qual, score, h_pol_compr_error_rate) in  all_representatives_sorted ]     #(read_cl_id, acc, seq, qual, score, h_pol_compr_error_rate)
    print("Clusters created in parallel (Iteration 1):", len(all_representatives_sorted))

    # chunk_size = int((len(all_representatives_sorted) )/ (args.nr_cores - 1)) + 1
    # print("Chunk sizes:", chunk_size)
    # read_batches_it2 =  [batch for batch in batch(all_representatives_sorted, chunk_size)] 
    # results = parallelize(reads_to_clusters_helper, [(({}, {}, read_batches_it2[i], p_emp_probs, args), {"minimizer_database" : {}}) for i in range(len(read_batches_it2))], args.nr_cores)
    # # all_cl, all_repr, all_minimizer_databases = [top_new_clusters],[top_new_seq_origins], [minimizer_database]
    # all_repr = [] #all_repr = [top_new_seq_origins]
    # all_cl = []
    # for new_clusters, new_cluster_origins, min_db in results: 
    #     all_cl.append(new_clusters)
    #     all_repr.append(new_cluster_origins)

    # all_representatives = merge_dicts(*all_repr)
    # all_representatives_sorted = sorted(list(all_representatives.values()), key = lambda x: x[4], reverse = True)
    # all_representatives_sorted = [ (read_cl_id, acc, seq, qual, score) for (read_cl_id, acc, seq, qual, score, h_pol_compr_error_rate) in  all_representatives_sorted ]     #(read_cl_id, acc, seq, qual, score, h_pol_compr_error_rate)
    # print("Clusters created in parallel (Iteration 2):", len(all_representatives_sorted))


    # TODO: MAYBE REMOVE top_representative before iterating on a single core the second time!
    # Process all representiatives on a single core to create final set of representatives

    all_clusters, all_seq_origins, minimizer_database = reads_to_clusters({}, {}, all_representatives_sorted, p_emp_probs, args, minimizer_database = {})
    print("Final clusters generated new:", len(all_seq_origins))
    assert len(all_clusters) == len(all_seq_origins)
    # Parallelize assignemt of reads to representatives

    read_batches = [batch for batch in batch(read_array, chunk_size)]
    print()
    print("READ ASSIGNMENT")
    print()
    all_clusters = {k: [] for k in all_clusters.keys()} # remove all accessions before read assignment
    results = parallelize(reads_to_clusters_helper, [((all_clusters, all_seq_origins, read_batches[i], p_emp_probs, args), {"minimizer_database" : minimizer_database, "task" : "read_assignment"}) for i in range(len(read_batches))], args.nr_cores)
    # final_clusters = representative_assignment(all_seq_origins_empty_clusters, all_seq_origins, read_array_batches[i], p_emp_probs, args, minimizer_database = minimizer_database)
    # all_cl = [all_clusters] #, [all_seq_origins]
    final_clusters = {k: [] for k in all_clusters.keys()}
    print("final cluster length", len(final_clusters))
    for clusters, cluster_origins, min_db in results: 
        print("Cluster length", len(clusters))
        for cl_id in clusters:
            if cl_id in final_clusters:
                final_clusters[cl_id].extend(clusters[cl_id]) 
            else:
                final_clusters[cl_id] = clusters[cl_id]
    print("final cluster length last", len(final_clusters))
    print("TOTAL TIME NEW:", time() - all_start)
    # return final_clusters, all_seq_origins

    old_start = time()
    old_clusters, old_seq_origins, top_minimizer_database = reads_to_clusters({}, {}, read_array, p_emp_probs, args, minimizer_database = {})
    print("Final clusters generated old:", len(old_clusters), len(old_seq_origins))
    print("TOTAL TIME OLD:", time() - old_start)

    sys.exit()






    # ## OLD CODE

    # read_batches = [read_array[i:len(read_array):args.nr_cores] for i in range(args.nr_cores)]
    
    # cluster_batches = []
    # cluster_seq_origin_batches = []
    # for batch in read_batches:
    #     tmp_clust = {}
    #     tmp_clust_origin = {}
    #     for i, acc, seq, qual, score in batch:
    #         tmp_clust[i] = [acc]
    #         tmp_clust_origin[i] = (i, acc, seq, qual, score)
    #     cluster_batches.append(tmp_clust)
    #     cluster_seq_origin_batches.append(tmp_clust_origin)

    # del read_array

    # # do clustering

    # ####### parallelize alignment #########
    # # pool = Pool(processes=mp.cpu_count())
    # original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    # signal.signal(signal.SIGINT, original_sigint_handler)
    # it = 1
    # while True:
    #     pool = Pool(processes=int(args.nr_cores/it))
    #     print("Iteration:", it)
    #     if args.nr_cores == 1:
    #         # print([len(cluster_batches[0][i]) for i in cluster_batches[0].keys()])
    #         Cluster, cluster_seq_origin = reads_to_clusters(cluster_batches[0], cluster_seq_origin_batches[0], read_batches[0], p_emp_probs, args)
    #         # print([len(Cluster[cl]) for cl in Cluster])
    #         assert len(Cluster) == len(cluster_seq_origin)
    #         break
    #     #     sys.exit()
    #     try:
    #         print([len(b) for b in read_batches])
            
    #         res = pool.map_async(reads_to_clusters_helper, [ ((cluster_batches[i], cluster_seq_origin_batches[i], read_batches[i], p_emp_probs, args), {}) for i in range(len(read_batches))] )
    #         cluster_results =res.get(999999999) # Without the timeout this blocking call ignores all signals.
    #     except KeyboardInterrupt:
    #         print("Caught KeyboardInterrupt, terminating workers")
    #         pool.terminate()
    #         sys.exit()
    #     else:
    #         # print("Normal termination")
    #         pool.close()
    #     pool.join()

    #     read_batches = []
    #     # H_batches = []
    #     cluster_batches = []
    #     cluster_seq_origin_batches = []
    #     if len(list(cluster_results)) == 1 :
    #         Cluster, cluster_seq_origin = cluster_results[0]
    #         break
    #     else:
    #         for i in range(0, len(cluster_results), 2): # merge read_batches, this is easy since by construction, all clusters have unique IDs
    #             # H = defaultdict(set)
    #             new_clusters1, cluster_seq_origin1 = cluster_results[i]
    #             assert len(new_clusters1) == len(cluster_seq_origin1)
    #             new_clusters2, cluster_seq_origin2 = cluster_results[i+1]

    #             cluster_seq_origin =  merge_two_dicts(cluster_seq_origin1, cluster_seq_origin2)
    #             Cluster =  merge_two_dicts(new_clusters1, new_clusters2)
    #             # for k in H1.keys():
    #             #     H[k].update(H1[k])
    #             # for k in H2.keys():
    #             #     H[k].update(H2[k])

    #             read_batches.append( [ (i, acc, seq, qual, score) for i, (i, acc, seq, qual, score, error_rate) in sorted(cluster_seq_origin.items(), key=lambda x: x[1][4], reverse=True)] )
                
    #             #### DIFF AFTER BUGFIX1 -- the iteration > 1 bug ####
    #             # H_batches.append(H)
    #             # H_batches.append({})
    #             #####################################################
    #             cluster_batches.append(Cluster)
    #             cluster_seq_origin_batches.append(cluster_seq_origin)
    #     it += 1

    # return Cluster, cluster_seq_origin



