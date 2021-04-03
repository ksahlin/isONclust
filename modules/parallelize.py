
import os
import math
import signal
from multiprocessing import Pool
import multiprocessing as mp
from time import time


from modules import help_functions
from modules import cluster


def reads_to_clusters_helper(arguments):
    for k,v in arguments.items():
        args, kwargs = v
    return cluster.reads_to_clusters(*args, **kwargs)


def merge_dicts(*dict_args):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result


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
        if batch:
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
            if batch:
                yield batch
        
        elif batch_type == "read_lengths_squared":
            tot_length = sum([ math.pow(len(seq),2) for i, b_i, acc, seq, qual, score in lst] )
            nt_chunk_size = int(tot_length/nr_cores) + 1
            batch = []
            curr_size = 0
            for info in lst:
                curr_size +=  math.pow(len(info[3]),2) 
                batch.append(info)
                if curr_size >= nt_chunk_size:
                    yield batch
                    batch = []
                    curr_size = 0
            if batch:
                yield batch



def print_intermediate_results(clusters, cluster_seq_origin, args, iter_nr):
    path = args.outfolder +"/{0}".format(iter_nr) 
    help_functions.mkdir_p( path )
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
    read_batches = [batch for batch in batch_list(read_array, num_batches, batch_type = args.batch_type )]
    print("Using total nucleotide batch sizes:", [sum([len(seq) for i, b_i, acc, seq, qual, score in b]) for b in read_batches] )
    print("Nr reads in batches:", [len(b)  for b in read_batches] )
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

            data = {i+1 :((cluster_batches[0], cluster_seq_origin_batches[0], read_batches[0], p_emp_probs, lowest_batch_index_db[0], 1, args), {})} 
            result = reads_to_clusters_helper(data) # { new_batch_index : (Cluster, cluster_seq_origin, H, new_batch_index)}
            Cluster, cluster_seq_origin, _, _ = result[1]
            print("Time elapesd clustering last iteration single core:", time() - start_cluster)
            return Cluster, cluster_seq_origin


        start_multi = time()
        pool = Pool(processes=int(num_batches))
        try:
            # print([len(b) for b in read_batches])
            data = [ {i+1 :((cluster_batches[i], cluster_seq_origin_batches[i], read_batches[i], p_emp_probs, lowest_batch_index_db[i], i+1, args), {})} for i in range(len(read_batches))]
            res = pool.map_async(reads_to_clusters_helper, data)
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
        all_minimizer_databases = {}
        for output_dict in cluster_results:
            print("New batch")
            for k, v in output_dict.items():
                new_clusters, new_representatives, minimizer_database_new, batch_index = v
                print("Batch index", k)
            # for new_clusters, new_representatives, minimizer_database_new, batch_index in cluster_results: 
                all_cl.append(new_clusters)
                all_repr.append(new_representatives)
                all_minimizer_databases[batch_index] =  minimizer_database_new

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

        it += 1
        read_batches = [batch for batch in batch_list(read_array, num_batches, batch_type = args.batch_type, merge_consecutive = True)]
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
            lowest_batch_index_db.append( all_minimizer_databases[lowest_batch_index])

        del all_minimizer_databases

