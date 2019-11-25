from __future__ import print_function
from functools import reduce
import os,sys
import argparse
from collections import defaultdict
import math
from collections import deque
import itertools
from operator import mul
import re

import parasail

from modules import help_functions


def get_kmer_minimizers(seq, k_size, w_size):
    # kmers = [seq[i:i+k_size] for i in range(len(seq)-k_size) ]
    w = w_size - k_size
    window_kmers = deque([seq[i:i+k_size] for i in range(w +1)])
    curr_min = min(window_kmers)
    minimizers = [ (curr_min, list(window_kmers).index(curr_min)) ]

    for i in range(w+1,len(seq) - k_size + 1):
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



def get_all_hits(minimizers, minimizer_database, read_cl_id):
    """
        Get all representatives ID's that shares matches with the minimizers in the read.  
    """
    hit_clusters_ids = defaultdict(int)
    hit_clusters_hit_index = defaultdict(list)
    hit_clusters_hit_positions = defaultdict(list)
    for i, (m, pos) in enumerate(minimizers): # iterating over minimizers from upstream to downstream in read
        if m in minimizer_database:
            for cl_id in minimizer_database[m]: 
                hit_clusters_ids[cl_id] += 1
                hit_clusters_hit_index[cl_id].append(i) # index of the minimizer among the coordinate sorted minimizers in the read
                hit_clusters_hit_positions[cl_id].append(pos) # positions of the minimizer among the coordinate sorted minimizers in the read

    if read_cl_id in hit_clusters_ids:
        del hit_clusters_ids[read_cl_id]
        del hit_clusters_hit_index[read_cl_id]
        del hit_clusters_hit_positions[read_cl_id]

    return hit_clusters_ids, hit_clusters_hit_index, hit_clusters_hit_positions




def get_best_cluster(read_cl_id, compressed_seq_len, hit_clusters_ids, hit_clusters_hit_positions, minimizers, nummber_of_minimizers, hit_clusters_hit_index, representatives, p_emp_probs, args):
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
        top_matches = sorted(hit_clusters_hit_positions.items(), key=lambda x: (len(x[1]), sum(x[1]), representatives[x[0]][2]),  reverse=True) #sorted(hit_clusters_ids.items(), key=lambda x: x[1],  reverse=True)
        top_hits = len(top_matches[0][1])
        nr_shared_kmers = top_hits
        if top_hits < args.min_shared:
            pass
        else:
            for tm in top_matches:
                cl_id = tm[0]
                nm_hits = len(tm[1])
                if nm_hits < args.min_fraction * top_hits or nm_hits < args.min_shared:
                    break

                # cl_size = len(hit_clusters_ids)
                minimizer_hit_positions = hit_clusters_hit_positions[cl_id]
                minimizer_hit_indices = hit_clusters_hit_index[cl_id]
                assert len(minimizer_hit_indices) == len(minimizer_hit_positions)
                _, _, _, _, _, _, error_rate_c = representatives[cl_id]
                _, _, _, _, _, _, error_rate_read = representatives[read_cl_id]
                p_error_in_kmers_emp =  1.0 - p_shared_minimizer_empirical(error_rate_read, error_rate_c, p_emp_probs)
                minimizer_error_probabilities = [p_error_in_kmers_emp]*nummber_of_minimizers
                total_mapped = 0
                # prev_mpos = 0
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
                    # is_covered = True
                    return cl_id, nm_hits, mapped_ratio

    return  best_cluster_id, nr_shared_kmers, mapped_ratio 


def parasail_block_alignment(s1, s2, k, match_id, match_score = 2, mismatch_penalty = -2, opening_penalty = 5, gap_ext = 1):
    user_matrix = parasail.matrix_create("ACGT", match_score, mismatch_penalty)
    result = parasail.sg_trace_scan_16(s1, s2, opening_penalty, gap_ext, user_matrix)
    if result.saturated:
        print("SATURATED!",len(s1), len(s2))
        result = parasail.sg_trace_scan_32(s1, s2, opening_penalty, gap_ext, user_matrix)
        print("computed 32 bit instead")

    # difference in how to obtain string from parasail between python v2 and v3... 
    if sys.version_info[0] < 3:
        cigar_string = str(result.cigar.decode).decode('utf-8')
    else:
        cigar_string = str(result.cigar.decode, 'utf-8')
    
    s1_alignment, s2_alignment = help_functions.cigar_to_seq(cigar_string, s1, s2)

    # Rolling window of matching blocks
    match_vector = [ 1 if n1 == n2 else 0 for n1, n2 in zip(s1_alignment, s2_alignment) ]    
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


def get_best_cluster_block_align(read_cl_id, representatives, hit_clusters_ids, hit_clusters_hit_positions, phred_char_to_p, args):
    best_cluster_id = -1
    top_matches = sorted(hit_clusters_hit_positions.items(), key=lambda x: (len(x[1]), sum(x[1]), representatives[x[0]][2]),  reverse=True) #sorted(hit_clusters_ids.items(), key=lambda x: x[1],  reverse=True)
    _, _, _, seq, r_qual, _, _ = representatives[read_cl_id]
    # print(top_matches)
    top_hits = len(top_matches[0][1])
    alignment_ratio = 0.0
    for tm in top_matches:
        cl_id = tm[0]
        nm_hits = len(tm[1])
        if nm_hits < top_hits:
            break
        _, _, _, c_seq, c_qual, _, _ = representatives[cl_id]

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

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def reads_to_clusters(clusters, representatives, sorted_reads, p_emp_probs, minimizer_database, new_batch_index, args):
    """
        Iterates throughreads in sorted order (w.r.t. score) and:
            1. homopolymenr compresses the read and obtain minimizers
            2. Finds the homopolymenr compressed error rate (if not computed in previous pass if more than 1 core specified to the program)
            3. Finds all the representatives that shares minimizers with the read
            4. Finds the best of the hits using mapping approach
            5. If no hit is found in 4. tries to align to representative with th most shared minimizers.
            6. Adds current read to representative, or makes it a new representative of a new cluster.
            7. If new representative: add the minimizers to the minimizer database
            8. Assign the actual reads to their new cluster and their new cluster representative (since all reads were initialized as their own representatives to deal with multiprocessing) 
    """

    ## For multiprocessing only
    prev_b_indices = [ prev_batch_index for (read_cl_id, prev_batch_index, acc, seq, qual, score) in sorted_reads ]
    lowest_batch_index = max(1, min(prev_b_indices))
    skip_count = prev_b_indices.count(lowest_batch_index)
    print("Saved: {0} iterations.".format(skip_count) )
    ###################################
    
    ## logging counters 
    aln_passed_criteria = 0
    mapped_passed_criteria = 0
    aln_called = 0
    ###################
    
    phred_char_to_p = {chr(i) : min( 10**( - (ord(chr(i)) - 33)/10.0 ), 0.79433)  for i in range(128)} # PHRED encoded quality character to prob of error. Need this locally if multiprocessing
    cluster_to_new_cluster_id = {}

    if args.print_output:
        eprint("Iteration\tNrClusters\tMinDbSize\tCurrReadId\tClusterSizes")

    for i, (read_cl_id, prev_batch_index, acc, seq, qual, score) in enumerate(sorted_reads):

        ## This if statement is only active in parallelization code 
        ## to keep track of already processed reads in previous iteration
        if prev_batch_index == lowest_batch_index:
            lst = list(representatives[read_cl_id])
            lst[1] = new_batch_index
            t = tuple(lst)
            representatives[read_cl_id] =  t # just updated batch index
            continue
        ##############################################################
        
        ################################################################################
        ############  Just for develop purposes, print some info to std out ############
        if i % args.print_output == 0: 
            inv_map = {}
            for k, v in cluster_to_new_cluster_id.items():
                inv_map.setdefault(v, set()).add(k)
            cl_tmp = sorted( [ 1 + sum([len(clusters[cl_id]) for cl_id in c ]) for c in inv_map.values() ], reverse= True)
            cl_tmp_nontrivial = [cl_size_tmp for cl_size_tmp in cl_tmp if cl_size_tmp > 1]
            eprint("{0}\t{1}\t{2}\t{3}\t{4}".format(i, len(cl_tmp_nontrivial), len(minimizer_database), "_".join(acc.split("_")[:-1]), ",".join([str(s_) for s_ in sorted(cl_tmp_nontrivial, reverse=True)])))
            # print("Processing read", i+1 , "seq length:", len(seq), "nr non-trivial clusters:", len(cl_tmp_nontrivial), "kmers stored:", len(minimizer_database))
            # print("Non trivial cluster sizes:", sorted(cl_tmp_nontrivial, reverse=True))

            # print("clust distr:", [c_len for c_len in cl_tmp if c_len > 100] )
            # depth = [len(nr_cl) for kmer, nr_cl in  sorted(minimizer_database.items(), key=lambda x: len(x[1]), reverse= True) if len(nr_cl) > 1 ]
            # print("Minimizer database depths:", depth)
            # print("Nr trivial clusters :", cl_sizes)
        ################################################################################
        ################################################################################

        # 1. homopolymenr compress read and obtain minimizers

        seq_hpol_comp = ''.join(ch for ch, _ in itertools.groupby(seq))
        if len(seq_hpol_comp) < args.k:
            print( "skipping read of length:", len(seq), "homopolymer compressed:", len(seq_hpol_comp), seq)
            continue 
        minimizers = get_kmer_minimizers(seq_hpol_comp, args.k, args.w)
        
        # 2. Find the homopolymer compressed error rate (else statement is the only one active in single core mode)

        if len(representatives[read_cl_id]) == 7: # we have already computed homopolymenr compressed error rate in previous iteration (if isONclust is called with multiple cores):
            lst = list(representatives[read_cl_id])
            lst[1] = new_batch_index
            t = tuple(lst)
            representatives[read_cl_id] =  t # just updated batch index
        else:
            # indices = [i for i, (n1,n2) in enumerate(zip(seq[:-1],seq[1:])) if n1 != n2] # indicies we want to take quality values from to get quality string of homopolymer compressed read 
            # indices.append(len(seq) - 1)
            # qualcomp = ''.join([qual[i] for i in indices])
            # assert len(seq_hpol_comp) == len(qualcomp)
            all_read_hpol_lengths = [len([c for c in g]) for ch, g in itertools.groupby(seq)]
            # print(all_read_hpol_lengths)
            qualcomp = []
            start = 0
            for h_len in all_read_hpol_lengths:
                q_max = min(qual[start: start + h_len], key = lambda x: phred_char_to_p[x])
                qualcomp.append(q_max)
                # if h_len > 2:
                #     print(qual[start: start + h_len], q_max)
                start += h_len
            qualcomp = "".join([q for q in qualcomp])
            assert len(seq_hpol_comp) == len(qualcomp)
            # print(qualcomp)
            # assert len(qualcomp) == len(qualcomp2)
            # print(qualcomp)
            # print(qualcomp2)
            # print()

            # compute the average error rate after compression
            poisson_mean = sum([ qualcomp.count(char_) * phred_char_to_p[char_] for char_ in set(qualcomp)])
            h_pol_compr_error_rate = poisson_mean/float(len(qualcomp))
            representatives[read_cl_id] = (read_cl_id, new_batch_index, acc, seq, qual, score, h_pol_compr_error_rate) # adding homopolymenr compressed error rate to info tuple of cluster origin sequence


        # 3. Find all the representatives with shared minimizers (this is the time consuming function for noisy and large datasets)

        hit_clusters_ids, hit_clusters_hit_index, hit_clusters_hit_positions = get_all_hits(minimizers, minimizer_database, read_cl_id)

        
        # 4. Finds the best of the hits using mapping approach

        best_cluster_id_m, nr_shared_kmers_m, mapped_ratio = get_best_cluster(read_cl_id, len(seq_hpol_comp), hit_clusters_ids, hit_clusters_hit_positions, minimizers, len(minimizers), hit_clusters_hit_index, representatives, p_emp_probs, args)
        
        
        # 5. If step 4 is unsuccessfull we try to align the read to the representative(s) with the most shared minimizers.

        if best_cluster_id_m >= 0:
            mapped_passed_criteria += 1

        if best_cluster_id_m < 0 and nr_shared_kmers_m >= args.min_shared:
            aln_called += 1
            best_cluster_id_a, nr_shared_kmers_a, error_rate_sum, s1_alignment, s2_alignment, alignment_ratio = get_best_cluster_block_align(read_cl_id, representatives, hit_clusters_ids, hit_clusters_hit_positions, phred_char_to_p, args)
            if best_cluster_id_a >= 0:
                aln_passed_criteria += 1

        else:
            best_cluster_id_a = -1


        # 6. Adds current read to representative, or makes it a new representative of a new cluster.
        
        best_cluster_id = max(best_cluster_id_m,best_cluster_id_a)
        if best_cluster_id >= 0:
            cluster_to_new_cluster_id[read_cl_id] = best_cluster_id

        # 7. If new representative: add the minimizers to the minimizer database. 

        else :  # Stays in current cluser, adding representative minimixers
            for m, pos in minimizers:
                if m in minimizer_database:
                    minimizer_database[m].add(read_cl_id)
                else:
                    minimizer_database[m] = set()
                    minimizer_database[m].add(read_cl_id)     
    
    
    # 8. Since all reads were initialized as their own representatives we need to reassign reads to their new representative, (this approach was implemented to deal with iterative assigment in the multiprocessing version)
    for read_cl_id in cluster_to_new_cluster_id:
        new_cl_id = cluster_to_new_cluster_id[read_cl_id]
        all_reads = clusters[read_cl_id]
        for read_acc in all_reads:
            clusters[new_cl_id].append(read_acc)
        del clusters[read_cl_id]
        # delete old origins
        del representatives[read_cl_id]
    ##########################

    print("Total number of reads iterated through:{0}".format(len(sorted_reads)))
    print("Passed mapping criteria:{0}".format(mapped_passed_criteria))
    print("Passed alignment criteria in this process:{0}".format(aln_passed_criteria))
    print("Total calls to alignment mudule in this process:{0}".format(aln_called))

    return { new_batch_index : (clusters, representatives, minimizer_database, new_batch_index)}


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





