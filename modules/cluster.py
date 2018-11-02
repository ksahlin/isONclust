
#BONUS: homopolymenr compress and do minimizers only
# store all k-mers of assigned reads in a set S
# Store a structure H with {k-mer: set(cluster ids)}
# for each read, assign the read to the cluster with the most matching k-mers if over a threshold, otherwise assign the read to a new cluster
# BONUS sensitivity: we can keep 3 different k-mer sizes in parallell to customize the length given the quality or the read and the highest sequence cuality in the set 
# If n is the number of reads, m is the number of clusters and k is the number of k-menrs in the read, the complexity is O(nk)


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


# def get_all_hits(minimizers, H, Clusters, read_cl_id):
#     hit_clusters_ids = defaultdict(int)
#     hit_clusters_start = defaultdict(int)
#     hit_clusters_stop = defaultdict(int)
#     hit_clusters_min_start_index = defaultdict(int)
#     hit_clusters_min_stop_index = defaultdict(int)

#     hit_clusters_prevpos = defaultdict(int)
#     hit_clusters_deltamax = defaultdict(int)
#     hit_clusters_hit_index = defaultdict(list)
#     hit_clusters_hit_positions = defaultdict(list)
#     # hit_cl_id_array = np.zeros((50000,), dtype=int)
#     nr_rep = 0
#     for i, (m, pos) in enumerate(minimizers):
#         if m in H:
#             if len(H[m]) > 100:
#                 nr_rep += 1
#                 # print("repetitive minimizer, masking")
#             for cl_id in H[m]: 
#                 hit_clusters_ids[cl_id] += 1
#                 hit_clusters_hit_index[cl_id].append(i)
#                 hit_clusters_hit_positions[cl_id].append(pos)

#                 if hit_clusters_ids[cl_id] == 1:
#                     hit_clusters_start[cl_id] = pos
#                     hit_clusters_min_start_index[cl_id] = i
#                     hit_clusters_prevpos[cl_id] = pos
#                     hit_clusters_deltamax[cl_id] = 0

#                 else:
#                     if pos - hit_clusters_prevpos[cl_id] > hit_clusters_deltamax[cl_id]:
#                         hit_clusters_deltamax[cl_id] =  pos - hit_clusters_prevpos[cl_id]

#                 hit_clusters_stop[cl_id] = pos
#                 hit_clusters_prevpos[cl_id] = pos
#                 hit_clusters_min_stop_index[cl_id] = i
    
#     if read_cl_id in hit_clusters_ids:
#         del hit_clusters_ids[read_cl_id]
#         del hit_clusters_min_start_index[read_cl_id]
#         del hit_clusters_min_stop_index[read_cl_id]
#         del hit_clusters_hit_index[read_cl_id]
#         del hit_clusters_hit_positions[read_cl_id]

#     return hit_clusters_ids, hit_clusters_min_start_index, hit_clusters_min_stop_index, hit_clusters_hit_index, hit_clusters_hit_positions


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
        # best_cluster_id, nr_shared_kmers = max(hit_clusters_ids.items(), key=lambda x: x[1])
        # print("here")
        top_matches = sorted(hit_clusters_ids.items(), key=lambda x: x[1],  reverse=True)
        top_hits = top_matches[0][1]
        nr_shared_kmers = top_hits
        if top_hits < args.min_shared:
            # print(top_hits)
            pass

        else:
            for tm in top_matches:
                cl_id = tm[0]
                nm_hits = tm[1]
                # print("ok")
                if nm_hits < args.min_fraction * top_hits or nm_hits < args.min_shared:
                    # print(nm_hits, args.min_fraction * top_hits, top_hits)
                    break

                cl_size = len(hit_clusters_ids)
                minimizer_hit_positions = hit_clusters_hit_positions[cl_id]
                minimizer_hit_indices = hit_clusters_hit_index[cl_id]
                assert len(minimizer_hit_indices) == len(minimizer_hit_positions)

                _, _, _, _, _, error_rate_c = cluster_seq_origin[cl_id]
                _, _, _, _, _, error_rate_read = cluster_seq_origin[read_cl_id]
                p_error_in_kmers_emp =  1.0 - p_shared_minimizer_empirical(error_rate_read, error_rate_c, p_emp_probs)
                # c_seq_quals_comp = ''.join([ c_qual[i] for i, (n1,n2) in enumerate(zip(c_seq[:-1], c_seq[1:])) if n1 != n2]) + c_qual[-1]
                # p_error_in_kmers_emp = 1.0 - p_shared_minimizer_empirical(seq_quals_comp, c_seq_quals_comp, p_emp_probs)
                # print("map:", p_error_in_kmers_emp)
                minimizer_error_probabilities = [p_error_in_kmers_emp]*nummber_of_minimizers
                # p_non_hits = [minimizer_error_probabilities[pos] for pos in minimizer_hit_indices]
                # assert len(p_non_hits) == len(minimizer_hit_positions)
                # check minimizera across read
                # is_covered = True
                total_mapped = 0
                prev_mpos = 0
                prob_all_errors_since_last_hit = [reduce(mul, minimizer_error_probabilities[: minimizer_hit_indices[0]], 1)] +  [ reduce(mul, minimizer_error_probabilities[hit_idx1+1: hit_idx2], 1) for hit_idx1, hit_idx2 in zip(minimizer_hit_indices[:-1], minimizer_hit_indices[1:]) ] + [reduce(mul, minimizer_error_probabilities[minimizer_hit_indices[-1]+1 : ], 1)]
                # print(minimizer_hit_indices)
                # print(minimizer_hit_positions)
                # print("P minimizer not shared:", minimizer_error_probabilities[0])
                # print("new:", p_error_in_kmers_emp, "old:",minimizer_error_probabilities[0])
                # print(prob_all_errors_since_last_hit)
                # print()



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
                # print("mapped ratio:", mapped_ratio)
                # print(total_mapped, seq_len, total_mapped /float(seq_len), nm_hits, num_min )
                    # if i == 0:
                    #     mapped_segment = m_pos
                    #     if m_pos - prev_mpos - 100 > args.max_minimizer_gap: # allow 100 extra in beginning
                    #         if prob_all_errors_since_last_hit[i] < 0.1:
                    #             # print("Beginning Breaking, but",  m_pos - prev_mpos, nm_hits, num_min, "p-error:", prob_all_errors_since_last_hit[i], minimizer_hit_indices[i], minimizer_error_probabilities[ : minimizer_hit_indices[i]] )#, minimizer_error_probabilities)
                    #             # is_covered = False
                    #             # break
                    #             pass
                    #         else:
                    #             total_mapped += mapped_segment
    
                    # elif m_pos - prev_mpos > args.max_minimizer_gap:
                    #     if prob_all_errors_since_last_hit[i] < 0.1:
                    #         # print("Breaking, but", m_pos - prev_mpos, nm_hits, num_min, "p-error:", prob_all_errors_since_last_hit[i],  minimizer_hit_indices[i-1:i+1], minimizer_error_probabilities[minimizer_hit_indices[i-1] + 1: minimizer_hit_indices[i] ] )
                    #         # is_covered = False
                    #         # break
                    #         pass
                    #     else:
                    #         total_mapped += mapped_segment
                    # else:
                    #     # just extend the mapped segment here
                    #     mapped_segment += m_pos - prev_mpos

                    # prev_mpos = m_pos
                # if is_covered:
                #     if seq_len - minimizer_hit_positions[-1] - 50 > args.max_minimizer_gap: # allow 50 extra in end
                #         if  prob_all_errors_since_last_hit[-1] < 0.1:
                #             # print("End Breaking, but", seq_len - prev_mpos, nm_hits, num_min, "p-error:", prob_all_errors_since_last_hit[-1], minimizer_hit_indices[i], minimizer_error_probabilities[minimizer_hit_indices[-1]: ] ) #, minimizer_error_probabilities)
                #             is_covered = False
                #             break
                
                if mapped_ratio > args.mapped_threshold:
                    is_covered = True
                    return cl_id, nm_hits, mapped_ratio
                    # break
                # else:
                #     # print(total_mapped, seq_len, total_mapped/float(seq_len)) #, prob_all_errors_since_last_hit) #, minimizer_hit_indices, minimizer_error_probabilities)
                #     is_covered = False
                #     break


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
        # print(seq)
        # print(c_seq)
        # print(gap_opening_penalty, match_id_tailored, error_rate_sum, len(seq), acc )
        (s1, s2, (s1_alignment, s2_alignment, alignment_ratio)) = parasail_block_alignment(seq, c_seq, args.k, match_id_tailored, opening_penalty = gap_opening_penalty,  )
        # print("Expected errors:", poisson_mean, poisson_mean2)
        if alignment_ratio >= args.aligned_threshold: #args.mapped_threshold:
            return cl_id, nm_hits,  error_rate_sum, s1_alignment, s2_alignment, alignment_ratio
        # else:
        #     return -1, 0,              mean_plus_two_stdvs_q1, mean_plus_two_stdvs_q2, lower_exp_identity, avg_id1, avg_id2, error_rate_sum, s1_alignment, s2_alignment, alignment_ratio

    return  best_cluster_id, 0,  -1, -1, -1, alignment_ratio


# def get_fraction_minimizers_destroyed(qual_string, minimizers, k_size):
#     minimizer_quals = [ qual_string[pos : pos +k_size] for m, pos in minimizers]
#     p_error_in_minimizers = [ 1 - reduce(mul, [ (1 - 10**( - (ord(char_) - 33)/10.0 ))**m.count(char_) for char_ in set(m)], 1) for m in minimizer_quals]
#     # print("HERE", p_error_in_minimizers, minimizer_quals)
#     # approximate this Poisson-Binomial distribution with Poisson, the mean will be expected number of of minimizers destroyed
#     # poisson_mean = sum(p_error_in_minimizers)
#     # mean_plus_two_stdvs = poisson_mean + 2*math.sqrt(poisson_mean) 
#     # fraction_minimizers_destroyed = min(0.9,  mean_plus_two_stdvs / float(len(minimizers))) # set upper limit to 0.9
#     return p_error_in_minimizers




# def archive_reads(H, Cluster, cluster_seq_origin, archived_reads, reads, mask_threshold, cluster_threshold, masked, next_pass_cleanup = False):
#     # get kmers for ech individual cluster

#     #TODO: Store this instead of recomputing..
#     kmers_in_cluster = defaultdict(set)
#     for k in H:
#         for cl_id in H[k]:
#             kmers_in_cluster[cl_id].add(k)

#     size_sorted_clusters = sorted(Cluster.items(), key=lambda x: len(x[1]) )
#     removed_clusters = []
#     saved_counter = 0 
#     for cl_id, cl_set in size_sorted_clusters:
#         cl_size = len(cl_set)
#         if cl_size > cluster_threshold:
#             break
#             # for acc in cl_set:
#             #     del reads[acc]
#         else:
#             if len(masked.intersection(kmers_in_cluster[cl_id]) ) >= 1 or next_pass_cleanup: # only remove the cluster if it has at least 1 high complexity minimizers
#                 removed_clusters.append(cl_id)
#                 for acc in cl_set:
#                     archived_reads[acc] = reads[acc]
#                     del reads[acc]
#             else:
#                 saved_counter +=1 
#     print("{0} clusters saved from removal".format(saved_counter) )

#     # remove the clusters and kmers of the archived reads
#     for cl_id in removed_clusters:
#         del Cluster[cl_id]
#         del cluster_seq_origin[cl_id]
#         del kmers_in_cluster[cl_id]

#     # Recreate kmer dict from remaining clusters
#     H = defaultdict(set)
#     for cl_id in kmers_in_cluster:
#         for k in kmers_in_cluster[cl_id]:
#             H[k].add(cl_id)

#     size_sorted_clusters_after_cleanup = sorted(H.items(), key=lambda x: len(x[1]) )
#     current_complexity_metric = len(size_sorted_clusters_after_cleanup[-100][1]) if len(size_sorted_clusters_after_cleanup) >=100 else 0
#     # current_complexity_metric = sum([len(cl_id) for k, cl_id in size_sorted_clusters_after_cleanup[-100:])/ float(len(size_sorted_clusters_after_cleanup[-100:]))
#     mask_threshold  = current_complexity_metric + 200
#     current_masked = 0
#     masked = set()
#     for k in H:
#         if len(H[k]) > mask_threshold:
#             current_masked += 1
#             masked.add(k)

#     print(len(archived_reads), "archived reads.", "New mask threshold:", mask_threshold, "Current kmers masked:", current_masked)
#     return H, archived_reads, current_masked, mask_threshold


def reads_to_clusters(Cluster, cluster_seq_origin, sorted_reads, p_emp_probs, args):
    """
        Iterates throughreads in sorted order (w.r.t. score) and:
            1. homopolymenr compresses the read
            2. Finds the homopolymenr compressed error rate (if not computed in previous pass if more than 1 core specified to the program)
            3. Finds all the representatives with shared minimizers (in "get_all_hits_new")
            4. Finds the best of the hits using mapping approach
            5. If no hit is found in 4. tries to align to representative with th most shared minimizers.
            6. Adds current read to representative, or makes it a new representative of a new cluster.
                6''. If new representative, and add the minimizers to the miniizer database H. 
    """
    H = {}
    print("USING w:{0}, k:{1}".format(args.w, args.k))
    phred_char_to_p = {chr(i) : min( 10**( - (ord(chr(i)) - 33)/10.0 ), 0.5)  for i in range(128)} # PHRED encoded quality character to prob of error. Need this locally if multiprocessing
    cluster_to_new_cluster_id = {}
    ## logging counters 
    aln_passed_criteria = 0
    mapped_passed_criteria = 0
    aln_called = 0
    # total_reads = 0
    ###################

    for i, (read_cl_id, acc, seq, qual, score) in enumerate(sorted_reads):
        # print(read_cl_id)
        # print(len(seq))
        
        ################################################################################
        ############  Just for develop purposes, print some info to std out ############
        if i%5000 == 0 and i > 0: 
            inv_map = {}
            for k, v in cluster_to_new_cluster_id.items():
                inv_map.setdefault(v, set()).add(k)
            cl_tmp = sorted( [ 1 + sum([len(Cluster[cl_id]) for cl_id in c ]) for c in inv_map.values() ], reverse= True)
            cl_tmp_nontrivial = [cl_size_tmp for cl_size_tmp in cl_tmp if cl_size_tmp > 1]
            print("Processing read", i, "seq length:", len(seq), "nr non-trivial clusters:", len(cl_tmp_nontrivial), "kmers stored:", len(H))
            print("clust distr:", [c_len for c_len in cl_tmp if c_len > 100] )
            depth = [len(nr_cl) for kmer, nr_cl in  sorted(H.items(), key=lambda x: len(x[1]), reverse= True)[:50]]
            print("Depth of H:", sum(depth)/float(len(depth)), depth)
        ################################################################################
        ################################################################################

        # homopolymenr compress read
        seq_hpol_comp = ''.join(ch for ch, _ in itertools.groupby(seq))
        
        if len(cluster_seq_origin[read_cl_id]) == 6: # we have already computed homopolymenr compressed error rate in previous iteration (if qtclust is called with multiple cores):
            pass
        else:
            # qualcomp = ''.join([ qual[i] for i, (n1,n2) in enumerate(zip(seq[:-1],seq[1:])) if n1 != n2]) + qual[-1]
            indices = [i for i, (n1,n2) in enumerate(zip(seq[:-1],seq[1:])) if n1 != n2] # indicies we want to take quality values from to get quality string of homopolymer compressed read 
            indices.append(len(seq) - 1)
            qualcomp = ''.join([qual[i] for i in indices])
            # assert qualcomp == qualcomp2
            # print(qualcomp == qualcomp2)
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
        # p_error_in_minimizers = get_fraction_minimizers_destroyed(qualcomp, minimizers, args.k)
        # p_error_in_kmers_appr = 1.0 - get_p_no_error_in_kmers_approximate(qualcomp,args.k)
        # print("P minimizernt shared:", p_error_in_kmers_appr)
        # p_error_in_minimizers2 = [p_error_in_kmers_appr]*len(minimizers)


        # hit_clusters_ids, hit_clusters_min_start_index, hit_clusters_min_stop_index, hit_clusters_hit_index, hit_clusters_hit_positions = get_all_hits(minimizers, H, Cluster, read_cl_id)
        hit_clusters_ids, hit_clusters_hit_index, hit_clusters_hit_positions = get_all_hits_new(minimizers, H, Cluster, read_cl_id)
        best_cluster_id_m, nr_shared_kmers_m, mapped_ratio = get_best_cluster(read_cl_id, len(seq_hpol_comp), hit_clusters_ids, hit_clusters_hit_positions, minimizers, len(minimizers), hit_clusters_hit_index, cluster_seq_origin, p_emp_probs, args)
        
        # if best_cluster_id_m > 0:
        #     best_cluster_id_a, nr_shared_kmers_a, mean_plus_two_stdvs_q1, mean_plus_two_stdvs_q2, lower_exp_identity, avg_id1, avg_id2, k, s1_alignment, s2_alignment, ed = get_best_cluster_align(seq, qual, cluster_seq_origin, hit_clusters_ids, args)
        #     print("alignment:",len(minimizers), nr_shared_kmers_m, best_cluster_id_a, ed, k )
        #     print(s1_alignment)
        #     print(s2_alignment)
        #     print()
        #     print()
        #     best_cluster_id_a, nr_shared_kmers_a, mean_plus_two_stdvs_q1, mean_plus_two_stdvs_q2, lower_exp_identity, avg_id1, avg_id2, k, s1_alignment, s2_alignment, ed = get_best_cluster_block_align(seq, qual, cluster_seq_origin, hit_clusters_ids, args)
        #     print("alignment:",len(minimizers), nr_shared_kmers_m, best_cluster_id_a, ed, k )
        #     print(s1_alignment)
        #     print(s2_alignment)
        if best_cluster_id_m >= 0:
            mapped_passed_criteria += 1

        if best_cluster_id_m < 0 and nr_shared_kmers_m >= args.min_shared:
            aln_called += 1
            # best_cluster_id_a, nr_shared_kmers_a, mean_plus_two_stdvs_q1, mean_plus_two_stdvs_q2, lower_exp_identity, avg_id1, avg_id2, k, s1_alignment, s2_alignment, ed = get_best_cluster_align(seq, qual, cluster_seq_origin, hit_clusters_ids, args)
            
            # print("mapping:",len(minimizers), nr_shared_kmers_m, "id:", best_cluster_id_m )
            # print("alignment:",len(minimizers), nr_shared_kmers_m, best_cluster_id_a, ed, k )
            # print(s1_alignment)
            # print(s2_alignment)
            # print()
            # print()
            best_cluster_id_a, nr_shared_kmers_a, error_rate_sum, s1_alignment, s2_alignment, alignment_ratio = get_best_cluster_block_align(read_cl_id, cluster_seq_origin, hit_clusters_ids, phred_char_to_p, args)
            if best_cluster_id_a >= 0:
                aln_passed_criteria += 1
            # print("alignment block:",len(minimizers), nr_shared_kmers_m, best_cluster_id_a, ed, k )
            # print("aln ratio", alignment_ratio)
            # print(s1_alignment)
            # print(s2_alignment)
            # print()
            # if len(minimizers) == 69 and nr_shared_kmers_m == 51:
            #     sys.exit()
        else:
            best_cluster_id_a = -1
        # if best_cluster_id_m == -1 and best_cluster_id != -1:
        #     print("mapping not assigned", len(minimizers), nr_shared_kmers_m, ed, p_error_in_kmers_appr)
        #     FN +=1
        #     tmp_out.write(">{0}\n{1}\n>{2}\n{3}\n".format(acc+ "FN", seq, acc+"_origin", "".join([c for c in s2_alignment if c != "-"]) ))
        #     print(mean_plus_two_stdvs_q1, mean_plus_two_stdvs_q2, lower_exp_identity, avg_id1, avg_id2, k)
        #     print(s1_alignment)
        #     print(s2_alignment)
        # if best_cluster_id_m != -1 and best_cluster_id == -1:
        #     print("Alignment not assigned", len(minimizers), nr_shared_kmers_m, ed, p_error_in_kmers_appr)
        #     FP +=1
        #     tmp_out.write(">{0}\n{1}\n>{2}\n{3}\n".format(acc+ "FP", seq, acc+"_origin", "".join([c for c in s2_alignment if c != "-"]) ))
        #     print(mean_plus_two_stdvs_q1, mean_plus_two_stdvs_q2, lower_exp_identity, avg_id1, avg_id2, k)
        #     print(s1_alignment)
        #     print(s2_alignment)
        
        # if best_cluster_id_m != best_cluster_id_a and best_cluster_id_m >= 0 and best_cluster_id_a >= 0:
        #     print("Different best!")
        #     print("mapping not assigned", len(minimizers), nr_shared_kmers_m, ed, p_error_in_kmers_appr)
        #     print(s1_alignment)
        #     print(s2_alignment)

        best_cluster_id = max(best_cluster_id_m,best_cluster_id_a)

        # if i > 160000:
        #     if best_cluster_id > 0:
        #         (c_id, acc, c_seq, qual, score) = cluster_seq_origin[best_cluster_id]
        #         c_seq_hpol_comp = ''.join(ch for ch, _ in itertools.groupby(c_seq))
        #         c_minimizers = get_kmer_minimizers(c_seq_hpol_comp, args.k, args.w)
        #         print(best_cluster_id_m, best_cluster_id_a, nr_shared_kmers_m, len(minimizers), len(c_minimizers), len(seq), len(seq_hpol_comp), len(c_seq_hpol_comp), i)
        #         print(minimizers)
        #         print(c_minimizers)
        #         print()

        
        if best_cluster_id >= 0:
            cluster_to_new_cluster_id[read_cl_id] = best_cluster_id

        else :  # Stays in current cluser, adding representative minimixers
            # print("here")
            # new_read_cluster_id = best_cluster_id if best_cluster_id >=0 else read_cl_id
            for m, pos in minimizers:
                if m in H:
                    H[m].add(read_cl_id)
                else:
                    H[m] = set()
                    H[m].add(read_cl_id)     

    
    # new no graph approach ####
    for read_cl_id in cluster_to_new_cluster_id:
        new_cl_id = cluster_to_new_cluster_id[read_cl_id]
        # (r_cl_id, acc, seq, qual, score, error_rate) = cluster_seq_origin[read_cl_id]
        # add all read acc to new origin
        all_reads = Cluster[read_cl_id]
        for read_acc in all_reads:
            Cluster[new_cl_id].append(read_acc)
        del Cluster[read_cl_id]
        # delete old origins
        del cluster_seq_origin[read_cl_id]
    ##########################


    print("PASS")
    print("Total number of reads iterated through:{0}".format(i+1))

    print("Passed mapping criteria:{0}".format(mapped_passed_criteria))
    print("Passed alignment criteria:{0}".format(aln_passed_criteria))
    print("Total calls to alignment mudule:{0}".format(aln_called))

    print("Percent passed mapping criteria:{0}".format( round(100*mapped_passed_criteria/float(i+1), 2) ))
    print("Percent passed alignment criteria total:{0}".format( round(100*aln_passed_criteria/float(i+1), 2) ))    
    print("Percent passed alignment criteria out of number of calls to this module:{0}".format( round(100*aln_passed_criteria/float(aln_called), 2) )) 

    print("total alingments called in core:", aln_called)
    print("total alingments passed:", aln_passed_criteria)
    return Cluster, cluster_seq_origin


# def get_p_no_error_in_kmers_approximate(qual_string, k):
#     poisson_mean = sum([ qual_string.count(char_) * 10**( - (ord(char_) - 33)/10.0 ) for char_ in set(qual_string)])
#     error_rate = poisson_mean/float(len(qual_string))
#     return (1.0 - error_rate)**k #1.0 - min(error_rate * k, 1.0)


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

# def p_shared_minimizer_empirical(seq_quals_comp, c_seq_quals_comp, p_emp_probs):
#     poisson_mean1 = sum([ seq_quals_comp.count(char_) * 10**( - (ord(char_) - 33)/10.0 ) for char_ in set(seq_quals_comp)])
#     error_rate1 = poisson_mean1/float(len(seq_quals_comp))
#     e1 = round(error_rate1, 2)
#     # print(e1)
#     if e1 > 0.15:
#         e1 = 0.15
#     if e1 < 0.01:
#         e1 = 0.01
#     poisson_mean2 = sum([ c_seq_quals_comp.count(char_) * 10**( - (ord(char_) - 33)/10.0 ) for char_ in set(c_seq_quals_comp)])
#     error_rate2 = poisson_mean2/float(len(c_seq_quals_comp))
#     e2 = round(error_rate2, 2)
#     # print(e2)
#     if e2 > 0.15:
#         e2 = 0.15
#     if e2 < 0.01:
#         e2 = 0.01
#     # print(e2)
#     p_kmer_shared = p_emp_probs[(e1,e2)]
#     return p_kmer_shared


def reads_to_clusters_helper(arguments):
    args, kwargs = arguments
    return reads_to_clusters(*args, **kwargs)

def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z

def next_power_of_2(x):  
    return 1 if x == 0 else 2**(x - 1).bit_length()


def cluster_seqs(read_array, p_emp_probs, args):
    # split sorted reads into batches
    num_batches = next_power_of_2(args.nr_cores)
    print("Using {0} batches.".format(num_batches))
    read_batches = [read_array[i:len(read_array):num_batches] for i in range(num_batches)]
    
    cluster_batches = []
    cluster_seq_origin_batches = []
    for batch in read_batches:
        tmp_clust = {}
        tmp_clust_origin = {}
        for i, acc, seq, qual, score in batch:
            tmp_clust[i] = [acc]
            tmp_clust_origin[i] = (i, acc, seq, qual, score)
        cluster_batches.append(tmp_clust)
        cluster_seq_origin_batches.append(tmp_clust_origin)

    # H_batches = [{} for i in range(num_batches) ]
    del read_array

    # do clustering

    ####### parallelize alignment #########
    # pool = Pool(processes=mp.cpu_count())
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.signal(signal.SIGINT, original_sigint_handler)
    it = 1
    while True:
        pool = Pool(processes=int(num_batches/it))
        print("Iteration:", it)
        if num_batches == 1:
            # print([len(cluster_batches[0][i]) for i in cluster_batches[0].keys()])
            Cluster, cluster_seq_origin = reads_to_clusters(cluster_batches[0], cluster_seq_origin_batches[0], read_batches[0], p_emp_probs, args)
            # print([len(Cluster[cl]) for cl in Cluster])
            assert len(Cluster) == len(cluster_seq_origin)
            break
        #     sys.exit()
        try:
            print([len(b) for b in read_batches])
            
            res = pool.map_async(reads_to_clusters_helper, [ ((cluster_batches[i], cluster_seq_origin_batches[i], read_batches[i], p_emp_probs, args), {}) for i in range(len(read_batches))] )
            cluster_results =res.get(999999999) # Without the timeout this blocking call ignores all signals.
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            pool.terminate()
            sys.exit()
        else:
            # print("Normal termination")
            pool.close()
        pool.join()

        read_batches = []
        # H_batches = []
        cluster_batches = []
        cluster_seq_origin_batches = []
        if len(list(cluster_results)) == 1 :
            Cluster, cluster_seq_origin = cluster_results[0]
            break
        else:
            for i in range(0, len(cluster_results), 2): # merge read_batches, this is easy since by construction, all clusters have unique IDs
                # H = defaultdict(set)
                new_clusters1, cluster_seq_origin1 = cluster_results[i]
                assert len(new_clusters1) == len(cluster_seq_origin1)
                new_clusters2, cluster_seq_origin2 = cluster_results[i+1]

                cluster_seq_origin =  merge_two_dicts(cluster_seq_origin1, cluster_seq_origin2)
                Cluster =  merge_two_dicts(new_clusters1, new_clusters2)
                # for k in H1.keys():
                #     H[k].update(H1[k])
                # for k in H2.keys():
                #     H[k].update(H2[k])

                read_batches.append( [ (i, acc, seq, qual, score) for i, (i, acc, seq, qual, score, error_rate) in sorted(cluster_seq_origin.items(), key=lambda x: x[1][4], reverse=True)] )
                
                #### DIFF AFTER BUGFIX1 -- the iteration > 1 bug ####
                # H_batches.append(H)
                # H_batches.append({})
                #####################################################
                cluster_batches.append(Cluster)
                cluster_seq_origin_batches.append(cluster_seq_origin)
        it += 1

    return Cluster, cluster_seq_origin



# def main(sorted_reads_fastq_file, args):
#     read_array = [ (i, acc, seq, qual, float(acc.split("_")[-1])) for i, (acc, (seq, qual)) in enumerate(readfq(open(sorted_reads_fastq_file, 'r')))]

#     clusters = cluster_seqs(read_array, args)

#     outfile = open(os.path.join(args.outfolder,  "pre_clusters.csv"), "w")
#     reads = {acc: (seq,qual) for (acc, (seq, qual)) in  readfq(open(sorted_reads, 'r'))}
#     output_index = 1
#     # nonclustered_outfile = open(os.path.join(reads_outfolder, "non_clustered.fa"), "w")
#     for c_id, all_read_acc in sorted(clusters.items(), key = lambda x: len(x[1]), reverse=True):
#         for r_acc in all_read_acc:
#             outfile.write("{0}\t{1}\n".format(c_id,r_acc))
#         if len(all_read_acc) > 1:
#             output_index += 1

#     print("Nr clusters larger than 1:", output_index) #, "Non-clustered reads:", len(archived_reads))



