import os,sys
import random
import itertools
from collections import deque

def get_kmer_minimizers(seq, k_size, w_size):
    # print(seq, k_size, w_size)
    # kmers = [seq[i:i+k_size] for i in range(len(seq)-k_size) ]
    w = min(len(seq), w_size) - k_size +1
    # print("w", w, k_size, w_size)
    window_kmers = deque([seq[i:i+k_size] for i in range(w)])
    # print(w, seq, k_size, w_size,window_kmers)

    curr_min = min(window_kmers)
    minimizers = [ (curr_min, list(window_kmers).index(curr_min)) ]

    for i in range(w, len(seq) - k_size):
        new_kmer = seq[i:i+k_size]
        # updateing window
        discarded_kmer = window_kmers.popleft()
        window_kmers.append(new_kmer)

        # we have discarded previous windows minimizer, look for new minimizer brute force
        if curr_min == discarded_kmer: 
            curr_min = min(window_kmers)
            minimizers.append( (curr_min, list(window_kmers).index(curr_min) + i - w + 1 ) )

        # Previous minimizer still in window, we only need to compare with the recently added kmer 
        elif new_kmer < curr_min:
            curr_min = new_kmer
            minimizers.append( (curr_min, i) )

    return minimizers


def get_minimizers_random_hash(array, w):
    window_kmers = deque([array[i] for i in range(w)])
    curr_min = min(window_kmers)
    minimizers = [ (curr_min, list(window_kmers).index(curr_min)) ]

    for i in range(w,len(array)):
        new_kmer = array[i]
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


def calc_p_share_random_hash(k_size, w_size, read1, read2, transcript):
    # create random hash table
    H = {}
    read1_h = []
    for i in range(len(read1)-k_size+1):
        if hash(read1[i:i+k_size]) in H:
            read1_h.append( H[ hash(read1[i:i+k_size]) ] )
        else:
            H[ hash(read1[i:i+k_size]) ] = random.getrandbits(128)
            read1_h.append( H[ hash(read1[i:i+k_size]) ] )

    read2_h = []
    for i in range(len(read2)-k_size+1):
        if hash(read2[i:i+k_size]) in H:
            read2_h.append( H[ hash(read2[i:i+k_size]) ] )
        else:
            H[ hash(read2[i:i+k_size]) ] = random.getrandbits(128)
            read2_h.append( H[ hash(read2[i:i+k_size]) ] )


    m1 = get_minimizers_random_hash(read1_h, w_size - k_size + 1 )
    m2 = get_minimizers_random_hash(read2_h, w_size - k_size + 1 )
    # t = get_minimizers_random_hash(t, w_size)
    # M1 = {}
    # for m, i in m1:
    #     if m in M1:
    #         M1[m].append(i)
    #     else:
    #         M1[m] = [i]
    M2 = {}
    for m, i in m2:
        if m in M2:
            M2[m].append(i)
        else:
            M2[m] = [i]

    shared_w_m2 = []
    for m, i1 in m1:
        if m in M2:
            for i2 in M2[m]:
                if -500 < i1 - i2 < 500:
                    # print(i1, i2)
                    shared_w_m2.append(1)
                    break
    return len(shared_w_m2)/float(len(m1))



def calc_p_share(k_size, w_size, read1, read2, t):
    m1 = get_kmer_minimizers(read1, k_size, w_size)
    m2 = get_kmer_minimizers(read2, k_size, w_size)
    t = get_kmer_minimizers(t, k_size, w_size)
    # print(m1, m2)
    M1 = {}
    for m, i in m1:
        if m in M1:
            M1[m].append(i)
        else:
            M1[m] = [i]
    M2 = {}
    for m, i in m2:
        if m in M2:
            M2[m].append(i)
        else:
            M2[m] = [i]

    # M1 = {m : i for m,i in m1.itmes()}
    # M2 = {m : i for m,i in m2.itmes()}

    shared_w_m2 = []
    for m, i1 in m1:
        if m in M2:
            for i2 in M2[m]:
                if -500 < i1 - i2 < 500:
                    # print(i1, i2)
                    shared_w_m2.append(1)
                    break
    # if len(shared_w_m2)/float(len(m1)) == 1.0:
    #     # print(len(shared_w_m2)/float(len(m1)))
    #     print(m1, m2,k_size, w_size,read1, read2 )
    # print(shared_w_m2)
    return len(shared_w_m2)/float(len(m1))


def calc_probs( t ):
    e1, e2  = t[0],t[1]
    pre_calc_probabilities = []
    for k in range(10, 31):
        for w in range(k, 101,5):
            probs_repl = []
            # probs2_repl = []
            # tmp = []
            for repl in range(1,1000):
                transcript = [random.choice("ACGT") for i in range(1000)] #initialization
                transcript = "".join([n for n in transcript])
                read1 = [n for n in transcript if random.random() > e1/2.0] #  delations
                read1 = [n + (random.choice("ACGT") if random.random() < e1/2.0 else "") for n in read1] # insertions
                read1 = "".join([n for n in read1])
                r1_comp = ''.join(ch for ch, _ in itertools.groupby(read1))

                read2 = [n for n in transcript if random.random() > e2/2.0] #  delations
                read2 = [n + (random.choice("ACGT") if random.random() < e2/2.0 else "") for n in read2] # insertions
                read2 = "".join([n for n in read2])
                r2_comp = ''.join(ch for ch, _ in itertools.groupby(read2))
                # print(e2,len(transcript), len(read1), len(read2), len(r1_comp), len(r2_comp) )
                p_shared_comp = calc_p_share(k, w, r1_comp, r2_comp, transcript)
                probs_repl.append(p_shared_comp)
            print(k,w, e1, e2, "avg:", sum(probs_repl)/len(probs_repl), "min:", min(probs_repl), "max", max(probs_repl) )
            # pre_calc_probabilities.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(k,w, sum(probs_repl)/len(probs_repl), e1, e2, "LH_HC" ))
            pre_calc_probabilities.append( (k,w, sum(probs_repl)/len(probs_repl), e1, e2) )
    return pre_calc_probabilities

from multiprocessing import Pool
p = Pool(2)
# error_rates = [0.15, 0.14, 0.13, 0.12, 0.11, 0.10, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01]
error_rates = [0.005, 0.002]
# results=[]
# for  e1, e2 in itertools.combinations_with_replacement(error_rates, 2):
#     r = calc_probs((e1, e2)) 
#     results.append(r)
results = p.map(calc_probs, [(e1, e2) for  e1, e2 in itertools.combinations_with_replacement(error_rates, 2)])
results_flattened = [item for sublist in results for item in sublist]
print(results_flattened)

outfile = open(sys.argv[1], "w")
outfile.write("L = {0}\n".format(results_flattened))
outfile.write("def read_empirical_p():\n")
outfile.write("    return L\n")
outfile.close()


########## OLD
# ####################################################
# ####################################################

# pre_calc_probabilities = open(sys.argv[1], "w")
# pre_calc_probabilities.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format("k", "w", "p", "e1", "e2", "hash" ))
# error_rates = [0.15, 0.14, 0.13, 0.12, 0.11, 0.10, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01]
# for  e1, e2 in itertools.combinations_with_replacement(error_rates, 2):
#     for k in range(10, 31):
#         for w in range(k, 101,5):
#             probs_repl = []
#             for repl in range(1,1000):
#                 transcript = [random.choice("ACGT") for i in range(1000)] #initialization
#                 transcript = "".join([n for n in transcript])
#                 read1 = [n for n in transcript if random.random() > e1/2.0] #  delations
#                 read1 = [n + (random.choice("ACGT") if random.random() < e1/2.0 else "") for n in read1] # insertions                
#                 read1 = "".join([n for n in read1])
#                 r1_comp = ''.join(ch for ch, _ in itertools.groupby(read1))

#                 read2 = [n for n in transcript if random.random() > e2/2.0] #  delations
#                 read2 = [n + (random.choice("ACGT") if random.random() < e2/2.0 else "") for n in read2] # insertions
#                 read2 = "".join([n for n in read2])
#                 r2_comp = ''.join(ch for ch, _ in itertools.groupby(read2))
#                 # p_shared = calc_p_share(k, w, read1, read2, transcript)
#                 # p_shared_rh = calc_p_share_random_hash(k, w, read1, read2, transcript)
#                 p_shared_comp = calc_p_share(k, w, r1_comp, r2_comp, transcript)
#                 probs_repl.append(p_shared_comp)
#             print(k,w, e1, e2, "avg:", sum(probs_repl)/len(probs_repl), "min:", min(probs_repl), "max", max(probs_repl) )
#             # pre_calc_probabilities.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(k,w, p_shared, e1, e2, "LH" ))
#             pre_calc_probabilities.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(k,w, sum(probs_repl)/len(probs_repl), e1, e2, "LH_HC" ))
#             # pre_calc_probabilities.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(k,w, p_shared_rh, e1, e2, "RH"))
# pre_calc_probabilities.close()
# ####################################################
# ####################################################


# ### PARAMETER COMBINATION FOR PLOT IN PAPER ########
# ####################################################

# transcript = [random.choice("ACGT") for i in range(100000)] #initialization
# transcript = "".join([n for n in transcript])
# pre_calc_probabilities = open(sys.argv[1], "w")
# pre_calc_probabilities.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format("k", "w", "P_minimizer_shared", "error_rate", "error_rate2", "hash" ))
# for e1 in [0.10, 0.05, 0.02]: 
#     read1 = [n for n in transcript if random.random() > e1/2.0] #  delations
#     read1 = [n + (random.choice("ACGT") if random.random() < e1/2.0 else "") for n in read1] # insertions
#     read1 = "".join([n for n in read1])
#     r1_comp = ''.join(ch for ch, _ in itertools.groupby(read1))
#     for e2 in [0.10, 0.05, 0.02]:
#         read2 = [n for n in transcript if random.random() > e2/2.0] #  delations
#         read2 = [n + (random.choice("ACGT") if random.random() < e2/2.0 else "") for n in read2] # insertions
#         read2 = "".join([n for n in read2])
#         r2_comp = ''.join(ch for ch, _ in itertools.groupby(read2))
#         for k in [10,15,20]:
#             for w in [20, 30, 50, 100]: #range(k, 100, 5):  # range(k+1, 100): 
#                 p_shared = calc_p_share(k, w, read1, read2, transcript)
#                 p_shared_rh = calc_p_share_random_hash(k, w, read1, read2, transcript)
#                 p_shared_comp = calc_p_share(k, w, r1_comp, r2_comp, transcript)
#                 print(k,w, round(p_shared, 3), round(p_shared_comp, 3),round(p_shared_rh, 3), e1, e2)
#                 pre_calc_probabilities.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(k,w, p_shared, e1, e2, "LH" ))
#                 pre_calc_probabilities.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(k,w, p_shared_comp, e1, e2, "LH_HC" ))
#                 pre_calc_probabilities.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(k,w, p_shared_rh, e1, e2, "RH"))

# pre_calc_probabilities.close()
# #################################################################
# #################################################################

