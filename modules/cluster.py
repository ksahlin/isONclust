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



def get_all_hits_new(minimizers, H, Clusters, read_cl_id):
    hit_clusters_ids = defaultdict(int)
    hit_clusters_hit_index = defaultdict(list)
    hit_clusters_hit_positions = defaultdict(list)
    for i, (m, pos) in enumerate(minimizers): # iterating over minimizers from upstream to downstream in read
        if m in H:
            for cl_id in H[m]: 
                # if cl_id not in Clusters:
                #     continue
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
                # if cl_id not in cluster_seq_origin:
                #     # print()
                #     # print("GAAAAAAAAAH")
                #     # print()
                #     # return -1, nr_shared_kmers, mapped_ratio 
                #     continue

                nm_hits = tm[1]
                if nm_hits < args.min_fraction * top_hits or nm_hits < args.min_shared:
                    break

                cl_size = len(hit_clusters_ids)
                minimizer_hit_positions = hit_clusters_hit_positions[cl_id]
                minimizer_hit_indices = hit_clusters_hit_index[cl_id]
                assert len(minimizer_hit_indices) == len(minimizer_hit_positions)
                _, _, _, _, _, _, error_rate_c = cluster_seq_origin[cl_id]
                _, _, _, _, _, _, error_rate_read = cluster_seq_origin[read_cl_id]
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
    _, _, _, seq, r_qual, _, _ = cluster_seq_origin[read_cl_id]
    # print(top_matches)
    top_hits = top_matches[0][1]
    aln_counter = 0
    alignment_ratio = 0.0
    for tm in top_matches:
        cl_id = tm[0]
        nm_hits = tm[1]
        # if cl_id not in cluster_seq_origin:
        #     # print("GGAGGAGAGGAAGAG2222")
        #     # return  -1, 0,  -1, -1, -1, alignment_ratio
        #     continue

        if nm_hits < top_hits:
            break
        aln_counter +=1
        _, _, _, c_seq, c_qual, _, _ = cluster_seq_origin[cl_id]

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


def reads_to_clusters(Cluster, cluster_seq_origin, sorted_reads, p_emp_probs, H, new_batch_index, args):
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
    # H = {}
    prev_b_indices = [ prev_batch_index for (read_cl_id, prev_batch_index, acc, seq, qual, score) in sorted_reads ]
    lowest_batch_index = max(1, min(prev_b_indices))
    skip_count = prev_b_indices.count(lowest_batch_index)
    print("Saved: {0} iterations.".format(skip_count) )
    # print(sorted(prev_b_indices))
    print("USING w:{0}, k:{1}".format(args.w, args.k))
    # print("HERE:",len(cluster_seq_origin))
    phred_char_to_p = {chr(i) : min( 10**( - (ord(chr(i)) - 33)/10.0 ), 0.5)  for i in range(128)} # PHRED encoded quality character to prob of error. Need this locally if multiprocessing
    cluster_to_new_cluster_id = {}
    ## logging counters 
    aln_passed_criteria = 0
    mapped_passed_criteria = 0
    aln_called = 0
    # total_reads = 0
    ###################

    for i, (read_cl_id, prev_batch_index, acc, seq, qual, score) in enumerate(sorted_reads):

        if prev_batch_index == lowest_batch_index:
            lst = list(cluster_seq_origin[read_cl_id])
            lst[1] = new_batch_index
            t = tuple(lst)
            cluster_seq_origin[read_cl_id] =  t # just updated batch index
            continue
        
        ################################################################################
        ############  Just for develop purposes, print some info to std out ############
        if i%10000 == 0 and i > 0: 
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
        
        if len(cluster_seq_origin[read_cl_id]) == 7: # we have already computed homopolymenr compressed error rate in previous iteration (if qtclust is called with multiple cores):
            lst = list(cluster_seq_origin[read_cl_id])
            lst[1] = new_batch_index
            t = tuple(lst)
            cluster_seq_origin[read_cl_id] =  t # just updated batch index
        else:
            indices = [i for i, (n1,n2) in enumerate(zip(seq[:-1],seq[1:])) if n1 != n2] # indicies we want to take quality values from to get quality string of homopolymer compressed read 
            indices.append(len(seq) - 1)
            qualcomp = ''.join([qual[i] for i in indices])
            assert len(seq_hpol_comp) == len(qualcomp)

            # compute the average error rate after compression
            poisson_mean = sum([ qualcomp.count(char_) * phred_char_to_p[char_] for char_ in set(qualcomp)])
            h_pol_compr_error_rate = poisson_mean/float(len(qualcomp))
            cluster_seq_origin[read_cl_id] = (read_cl_id, new_batch_index, acc, seq, qual, score, h_pol_compr_error_rate) # adding homopolymenr compressed error rate to info tuple of cluster origin sequence
        
        # get minimizers
        if len(seq_hpol_comp) < args.k:
            print( "skipping read of length:", len(seq), "homopolymer compressed:", len(seq_hpol_comp), seq)
            continue


        minimizers = get_kmer_minimizers(seq_hpol_comp, args.k, args.w)
        hit_clusters_ids, hit_clusters_hit_index, hit_clusters_hit_positions = get_all_hits_new(minimizers, H, Cluster, read_cl_id)
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
        
        if best_cluster_id >= 0:
            cluster_to_new_cluster_id[read_cl_id] = best_cluster_id

        else :  # Stays in current cluser, adding representative minimixers
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


    # print("PASS")
    # print("Total number of reads iterated through:{0}".format(i+1))

    # print("Passed mapping criteria:{0}".format(mapped_passed_criteria))
    # print("Passed alignment criteria in this process:{0}".format(aln_passed_criteria))
    # print("Total calls to alignment mudule in this process:{0}".format(aln_called))

    # print("Percent passed mapping criteria:{0}".format( round(100*mapped_passed_criteria/float(i+1), 2) ))
    # print("Percent passed alignment criteria total:{0}".format( round(100*aln_passed_criteria/float(i+1), 2) ))    
    # if aln_called > 0:
    #     print("Percent passed alignment criteria out of number of calls to the alignment module:{0}".format(round(100*aln_passed_criteria/float(aln_called), 2) )) 
    # print("PIckled:", get_pickled_memory((Cluster, cluster_seq_origin, H, new_batch_index)))
    # print("PIckled2:", get_pickled_memory({ new_batch_index : (Cluster, cluster_seq_origin, H, new_batch_index)}))
    return { new_batch_index : (Cluster, cluster_seq_origin, H, new_batch_index)}


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

# def reads_to_clusters_helper(arguments):
#     args, kwargs = arguments
#     return reads_to_clusters(*args, **kwargs)

def reads_to_clusters_helper2(arguments):
    for k,v in arguments.items():
        args, kwargs = v
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


def next_power_of_2(x):  
    return 1 if x == 0 else 2**(x - 1).bit_length()


def batch_list(lst, nr_cores=1, batch_type = "nr_reads" , merge_consecutive = False ):
    if merge_consecutive:
        batch_id = 2
        batch = []
        for info in lst:
            if info[1] <= batch_id:
                batch.append(info)
            else: # first sequence in new batch
                yield batch
                batch_id += 2
                batch = []
                batch.append(info)
        yield batch

    else:
        if batch_type == "nr_reads":
            l = len(lst)
            chunk_size = int(l/nr_cores) + 1
            for ndx in range(0, l, chunk_size):
                yield lst[ndx:min(ndx + chunk_size, l)]
        
        elif batch_type == "total_nt":
            tot_length = sum([len(seq) for i, b_i, acc, seq, qual, score in lst] )
            nt_chunk_size = int(tot_length/nr_cores) + 1

            batch = []
            curr_size = 0
            for info in lst:
                curr_size += len(info[3])
                batch.append(info)
                if curr_size >= nt_chunk_size:
                    yield batch
                    batch = []
                    curr_size = 0
            yield batch   
        
        elif batch_type == "weighted":
            tot_length = sum([math.sqrt(len(seq)) for i, b_i, acc, seq, qual, score in lst] )
            nt_chunk_size = int(tot_length/nr_cores) + 1
            batch = []
            curr_size = 0
            for info in lst:
                curr_size += math.sqrt(len(info[3]))
                batch.append(info)
                if curr_size >= nt_chunk_size:
                    yield batch
                    batch = []
                    curr_size = 0
            yield batch               




# def batch_wrt_total_nucl_length(lst, nr_cores=1):
#     tot_length = sum([len(seq) for i, b_i, acc, seq, qual, score in lst] )
#     nt_chunk_size = int(tot_length/nr_cores) + 1

#     batch = []
#     curr_size = 0
#     for info in lst:
#         curr_size += len(info[3])
#         batch.append(info)
#         if curr_size >= nt_chunk_size:
#             yield batch
#             batch = []
#             curr_size = 0
#     yield batch

    # l = len(lst)
    # for ndx in range(0, l, n):
    #     yield lst[ndx:min(ndx + n, l)]


#### TMP MEMORY CHECK
import pickle
import sys
# from collections import OrderedDict

def get_pickled_memory(data):
    return sum(sys.getsizeof(pickle.dumps(d)) for d in data)
    # print("pikled size:", sum(sys.getsizeof(pickle.dumps(d)) for d in data))
    # numprocs = 10
    # a = ['a' for i in range(1000000)]
    # b = [a+[] for i in range(100)]

    # data1 = [b+[] for i in range(numprocs)]
    # data2 = [data1+[]] + ['1' for i in range(numprocs-1)]
    # data3 = [['1'] for i in range(numprocs)]
    # sizes = OrderedDict()
    # for idx, data in enumerate((data1, data2, data3)):
    #     sizes['data{}'.format(idx+1)] = sum(sys.getsizeof(pickle.dumps(d))
    #                                             for d in data)

    # for k, v in sizes.items():
    #     print("{}: {}".format(k, v))

def mkdir_p(path):
    try:
        os.makedirs(path)
        print("creating", path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def print_intermediate_results(clusters, cluster_seq_origin, args, iter_nr):
    path = args.outfolder +"/{0}".format(iter_nr) 
    mkdir_p( path )
    outfile = open(os.path.join(path, "pre_clusters.csv"), "w")
    nontrivial_cluster_index = 0
    for c_id, all_read_acc in sorted(clusters.items(), key = lambda x: len(x[1]), reverse=True):
        for r_acc in all_read_acc:
            outfile.write("{0}\t{1}\n".format(c_id, "_".join([item for item in r_acc.split("_")[:-1]]) ))
        if len(all_read_acc) > 1:
            nontrivial_cluster_index += 1
    print("Nr clusters larger than 1:", nontrivial_cluster_index) #, "Non-clustered reads:", len(archived_reads))
    print("Nr clusters (all):", len(clusters)) #, "Non-clustered reads:", len(archived_reads))


    origins_outfile = open(os.path.join(path,  "cluster_origins.csv"), "w")
    for cl_id, all_read_acc in sorted(clusters.items(), key = lambda x: len(x[1]), reverse=True):
        read_cl_id, b_i, acc, c_seq, c_qual, score, error_rate = cluster_seq_origin[cl_id]
        origins_outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(read_cl_id, acc, c_seq, c_qual, score, error_rate)) 
    outfile.close()
    origins_outfile.close()






def parallel_clustering(read_array, p_emp_probs, args):

    num_batches = args.nr_cores 
    # prev_nr_repr = len(read_array)

    read_batches = [batch for batch in batch_list(read_array, num_batches, batch_type = "weighted" )]
    print("Using total nucleotide batch sizes:", [sum([len(seq) for i, b_i, acc, seq, qual, score in b]) for b in read_batches] )
    print("Using nr reads batch sizes:", [len(b)  for b in read_batches] )
    cluster_batches = []
    cluster_seq_origin_batches = []
    lowest_batch_index_db = []
    for batch in read_batches:
        tmp_clust = {}
        tmp_clust_origin = {}
        for i, b_i, acc, seq, qual, score in batch:
            tmp_clust[i] = [acc]
            tmp_clust_origin[i] = (i, b_i, acc, seq, qual, score)
        cluster_batches.append(tmp_clust)
        cluster_seq_origin_batches.append(tmp_clust_origin)
        lowest_batch_index_db.append({})
    del read_array

    ####### parallelize alignment #########
    # pool = Pool(processes=mp.cpu_count())
    original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.signal(signal.SIGINT, original_sigint_handler)
    try:
        mp.set_start_method('spawn')
        print("Environment set:", mp.get_context())
    except RuntimeError:
        print("Environment already set:", mp.get_context())
    it = 1
    while True:
        # Structure up batches
        print()
        print("ITERATION", it)
        print("Using {0} batches.".format(num_batches))

        if len(read_batches) == 1:
            start_cluster = time()

            data2 = {i+1 :((cluster_batches[0], cluster_seq_origin_batches[0], read_batches[0], p_emp_probs, lowest_batch_index_db[0], 1, args), {})} 
            result = reads_to_clusters_helper2(data2) # { new_batch_index : (Cluster, cluster_seq_origin, H, new_batch_index)}
            Cluster, cluster_seq_origin, H, new_batch_index = result[1]
            print("Time elapesd clustering last iteration single core:", time() - start_cluster)
            return Cluster, cluster_seq_origin


        start_multi = time()
        pool = Pool(processes=int(num_batches))
        try:
            print([len(b) for b in read_batches])

            ############################################
            # data1 = [ ((cluster_batches[i], cluster_seq_origin_batches[i], read_batches[i], p_emp_probs, lowest_batch_index_db[i], i+1, args), {}) for i in range(len(read_batches))]
            # print("Size:{0}Mb".format( [round( get_pickled_memory(d)/float(1000000), 2) for d in data1 ] ))
            data2 = [ {i+1 :((cluster_batches[i], cluster_seq_origin_batches[i], read_batches[i], p_emp_probs, lowest_batch_index_db[i], i+1, args), {})} for i in range(len(read_batches))]
            # data2 = [ {i+1 :(([0]*1000000,0), {})} for i in range(2)]
            # data2 = [ (([0]*2000000000,0), {}) for i in range(2)]
            print("Size:{0}Mb".format( [round( get_pickled_memory(d)/float(1000000), 2) for d in data2 ] ))
            # sys.exit()
            # for data in [ ((cluster_batches[i], cluster_seq_origin_batches[i], read_batches[i], p_emp_probs, args), {}) for i in range(len(read_batches))]:
            #     print("Size:{0}Mb".format(round( get_pickled_memory(data)/float(1000000), 2) ))
            #     data2 = {"k" :data, 'k2': []}
            #     print("Size:{0}Mb".format(round( get_pickled_memory(data2)/float(1000000), 7) ))
            #     print("Size:{0}Mb".format( [round( get_pickled_memory(d)/float(1000000), 2) for d in data for s in d] ))
            #     print()
            # sys.exit()
            ############################################

            res = pool.map_async(reads_to_clusters_helper2, data2)
            cluster_results =res.get(999999999) # Without the timeout this blocking call ignores all signals.
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            pool.terminate()
            sys.exit()
        else:
            pool.close()
        pool.join()

        print("Time elapesd multiprocessing:", time() - start_multi)

        start_joining = time()
        all_repr = [] # all_repr = [top_new_seq_origins]
        all_cl = []
        all_H = {}
        for output_dict in cluster_results:
            print("New batch")
            for k, v in output_dict.items():
                new_clusters, new_cluster_origins, H_new, batch_index = v
                print("Batch index", k)
            # for new_clusters, new_cluster_origins, H_new, batch_index in cluster_results: 
                all_cl.append(new_clusters)
                all_repr.append(new_cluster_origins)
                all_H[batch_index] =  H_new

        all_clusters = merge_dicts(*all_cl)
        all_representatives = merge_dicts(*all_repr)
        read_array =  [ (i, b_index, acc, seq, qual, score) for i, (i, b_index, acc, seq, qual, score, error_rate) in sorted(all_representatives.items(), key=lambda x: x[1][5], reverse=True)] 
        new_nr_repr = len(read_array)
        print("number of representatives left to cluster:", new_nr_repr)
        print("Time elapesd joining clusters:", time() - start_joining)

        # Determine new number of batches
        if num_batches == 1:
            return all_clusters, all_representatives
        else:
            print_intermediate_results(all_clusters, all_representatives, args, it)

        # elif new_nr_repr/float(prev_nr_repr) > 0.8:
        #     num_batches = 1
        # else:
        #     num_batches = max(1, int(round(num_batches*(new_nr_repr/float(prev_nr_repr)))) )
        # prev_nr_repr = new_nr_repr
        it += 1
        read_batches = [batch for batch in batch_list(read_array, num_batches, batch_type = "weighted", merge_consecutive = True)]
        num_batches = len(read_batches)
        print("Batches after pairwise consecutive merge:", num_batches)
        print("Using total nucleotide batch sizes:", [sum([len(seq) for i, b_i, acc, seq, qual, score in b]) for b in read_batches] )
        print("Using nr reads batch sizes:", [len(b)  for b in read_batches] )
        cluster_batches = []
        cluster_seq_origin_batches = []
        lowest_batch_index_db = []
        for batch in read_batches:
            tmp_clust = {}
            tmp_clust_origin = {}
            lowest_batch_index = min( [ prev_batch_index for (read_cl_id, prev_batch_index, acc, seq, qual, score) in batch ] )
            for i, b_i, acc, seq, qual, score in batch:
                tmp_clust[i] = all_clusters[i]
                tmp_clust_origin[i] = all_representatives[i]
            cluster_batches.append(tmp_clust)
            cluster_seq_origin_batches.append(tmp_clust_origin)
            lowest_batch_index_db.append( all_H[lowest_batch_index])

        del all_H


def single_clustering(read_array, p_emp_probs, args):

    num_batches = args.nr_cores 
    # prev_nr_repr = len(read_array)
    start_cluster = time()
    read_batches = [batch for batch in batch_list(read_array, num_batches, batch_type = "weighted" )]
    print("Using total nucleotide batch sizes:", [sum([len(seq) for i, b_i, acc, seq, qual, score in b]) for b in read_batches] )
    print("Using nr reads batch sizes:", [len(b)  for b in read_batches] )
    cluster_batches = []
    cluster_seq_origin_batches = []
    lowest_batch_index_db = []
    for batch in read_batches:
        tmp_clust = {}
        tmp_clust_origin = {}
        for i, b_i, acc, seq, qual, score in batch:
            tmp_clust[i] = [acc]
            tmp_clust_origin[i] = (i, b_i, acc, seq, qual, score)
        cluster_batches.append(tmp_clust)
        cluster_seq_origin_batches.append(tmp_clust_origin)
        lowest_batch_index_db.append({})
    del read_array

    print("Using 1 batch.")

    data2 = {i+1 :((cluster_batches[0], cluster_seq_origin_batches[0], read_batches[0], p_emp_probs, lowest_batch_index_db[0], 1, args), {})} 
    result = reads_to_clusters_helper2(data2) # { new_batch_index : (Cluster, cluster_seq_origin, H, new_batch_index)}
    Cluster, cluster_seq_origin, H, new_batch_index = result[1]

    print("Time elapesd clustering:", time() - start_cluster)
    return Cluster, cluster_seq_origin


def cluster_seqs(read_array, p_emp_probs, args):
    if args.nr_cores > 1:
        all_clusters, all_representatives = parallel_clustering(read_array, p_emp_probs, args)
    else:
        all_clusters, all_representatives = single_clustering(read_array, p_emp_probs, args)

    return all_clusters, all_representatives




