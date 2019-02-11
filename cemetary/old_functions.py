def homopolymer_compress(seq):
    corr = [ n1  for n1,n2 in zip(seq[:-1], seq[1:]) if n1 != n2 ]
    #last base corner case
    if seq[-1] != seq[-2]:
        corr.append(seq[-1])
    return "".join([nt for nt in corr])



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
calling code:
# print("Size:{0}Mb".format( [round( get_pickled_memory(d)/float(1000000), 2) for d in data ] ))
def get_pickled_memory(data):
    return sum(sys.getsizeof(pickle.dumps(d)) for d in data)
    # print("pikled size:", sum(sys.getsizeof(pickle.dumps(d)) for d in data))
    # numprocs = 10
    # a = ['a' for i in range(1000000)]
    # b = [a+[] for i in range(100)]

    # data1 = [b+[] for i in range(numprocs)]
    # data = [data1+[]] + ['1' for i in range(numprocs-1)]
    # data3 = [['1'] for i in range(numprocs)]
    # sizes = OrderedDict()
    # for idx, data in enumerate((data1, data, data3)):
    #     sizes['data{}'.format(idx+1)] = sum(sys.getsizeof(pickle.dumps(d))
    #                                             for d in data)

    # for k, v in sizes.items():
    #     print("{}: {}".format(k, v))


# code for writing statistics for the paper, called in reads_to_clusters function
    # print("PASS")
    # print("Total number of reads iterated through:{0}".format(i+1))

    # print("Passed mapping criteria:{0}".format(mapped_passed_criteria))
    # print("Passed alignment criteria in this process:{0}".format(aln_passed_criteria))
    # print("Total calls to alignment mudule in this process:{0}".format(aln_called))

    # print("Percent passed mapping criteria:{0}".format( round(100*mapped_passed_criteria/float(i+1), 2) ))
    # print("Percent passed alignment criteria total:{0}".format( round(100*aln_passed_criteria/float(i+1), 2) ))    
    # if aln_called > 0:
    #     print("Percent passed alignment criteria out of number of calls to the alignment module:{0}".format(round(100*aln_passed_criteria/float(aln_called), 2) )) 
    # print("PIckled:", get_pickled_memory((Cluster, cluster_seq_origin, minimizer_database, new_batch_index)))
    # print("PIckled2:", get_pickled_memory({ new_batch_index : (Cluster, cluster_seq_origin, minimizer_database, new_batch_index)}))



