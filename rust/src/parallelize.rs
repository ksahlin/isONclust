//! `modules/parallelize.py` -- the `--t > 1` path.
//!
//! **This is a different algorithm, not a parallelised version of the same one.**
//! Reads are cut into batches, each batch is clustered independently, the
//! surviving representatives are pooled and re-batched, and the whole thing
//! repeats until one batch remains. `--t` therefore changes the answer: on the
//! smoke corpus `--t 1` gives 40 clusters and `--t 8` gives 35. See PORTING.md,
//! Finding 3.
//!
//! **The batches are independent, so execution order is free.** Each worker
//! receives its own cluster map, representative map, read slice and minimizer
//! database, and returns new ones; nothing is shared. `pool.map_async` preserves
//! result order, so running the batches sequentially gives byte-identical output
//! to running them in parallel. That is why this is written sequentially first:
//! the batching structure is what has to be exact, and threading it afterwards
//! is behaviour-neutral.

use crate::sweep::{self, OrderedClusters, ReadInfo, SweepParams, SweepRead};
use rustc_hash::FxHashMap;
use std::collections::HashMap;

/// `batch_list`'s three real modes.
///
/// Note the CLI's help advertises `"weighted"`, which **no branch implements** --
/// the third mode is spelled `read_lengths_squared`. Passing the documented
/// value yields no batches at all and the reference then dies with
/// `ValueError: Number of processes must be at least 1`. PORTING.md, Finding 15.
#[derive(Clone, Copy, PartialEq)]
pub enum BatchType {
    NrReads,
    TotalNt,
    ReadLengthsSquared,
}

impl BatchType {
    pub fn parse(s: &str) -> Option<BatchType> {
        match s {
            "nr_reads" => Some(BatchType::NrReads),
            "total_nt" => Some(BatchType::TotalNt),
            "read_lengths_squared" => Some(BatchType::ReadLengthsSquared),
            _ => None,
        }
    }
}

/// `batch_list(lst, nr_cores, batch_type)` -- the first-pass split.
pub fn batch_list(reads: &[SweepRead], nr_cores: usize, bt: BatchType) -> Vec<Vec<SweepRead>> {
    let mut out: Vec<Vec<SweepRead>> = Vec::new();
    match bt {
        BatchType::NrReads => {
            // `chunk_size = int(l / nr_cores) + 1`, so the last chunk can be
            // short and there can be fewer chunks than cores.
            let l = reads.len();
            let chunk = l / nr_cores.max(1) + 1;
            let mut ndx = 0usize;
            while ndx < l {
                out.push(reads[ndx..(ndx + chunk).min(l)].to_vec());
                ndx += chunk;
            }
        }
        BatchType::TotalNt | BatchType::ReadLengthsSquared => {
            let weight = |r: &SweepRead| -> f64 {
                let n = r.seq.len() as f64;
                if bt == BatchType::TotalNt {
                    n
                } else {
                    // math.pow(len, 2)
                    n * n
                }
            };
            // `tot_length` is a sum of ints for total_nt and of floats for the
            // squared variant; `int(tot/nr_cores) + 1` truncates either way.
            let total: f64 = reads.iter().map(weight).sum();
            let chunk = (total / nr_cores.max(1) as f64) as i64 + 1;
            let mut batch: Vec<SweepRead> = Vec::new();
            let mut curr = 0f64;
            for r in reads {
                curr += weight(r);
                batch.push(r.clone());
                if curr >= chunk as f64 {
                    out.push(std::mem::take(&mut batch));
                    curr = 0.0;
                }
            }
            if !batch.is_empty() {
                out.push(batch);
            }
        }
    }
    out
}

/// `batch_list(..., merge_consecutive=True)` -- the re-batching between rounds.
///
/// This walks the score-ordered representatives and closes a batch whenever a
/// read's *previous batch index* exceeds a threshold that starts at 2 and rises
/// by 2 each time. The reads are in score order, not batch order, so the group
/// boundaries fall wherever the batch indices happen to cross the threshold --
/// arbitrary, but deterministic, and it is what decides which minimizer database
/// each group inherits.
pub fn batch_list_merge_consecutive(reads: &[SweepRead]) -> Vec<Vec<SweepRead>> {
    let mut out: Vec<Vec<SweepRead>> = Vec::new();
    let mut batch_id: i64 = 2;
    let mut batch: Vec<SweepRead> = Vec::new();
    for r in reads {
        if r.prev_batch_index <= batch_id {
            batch.push(r.clone());
        } else {
            // The reference yields the batch even when it is empty.
            out.push(std::mem::take(&mut batch));
            batch_id += 2;
            batch.push(r.clone());
        }
    }
    if !batch.is_empty() {
        out.push(batch);
    }
    out
}

/// What a completed round hands back.
pub struct ParallelResult {
    pub clusters: OrderedClusters,
    pub representatives: FxHashMap<usize, ReadInfo>,
    /// One entry per merge iteration, in order: `(pre_clusters, cluster_origins)`.
    pub intermediates: Vec<(String, String)>,
    pub mapped_passed: usize,
    pub aln_passed: usize,
    pub aln_called: usize,
    /// Reads whose homopolymer-compressed length was under k. The reference
    /// prints one line per read to stdout; a count is more useful and stdout is
    /// not in the byte-identity contract.
    pub skipped_short: usize,
    pub times: sweep::StageTimes,
}

/// `parallel_clustering`.
pub fn parallel_clustering(
    read_array: &[SweepRead],
    nr_cores: usize,
    bt: BatchType,
    table: &crate::p_emp::Table,
    p: SweepParams,
    // `sorted.fastq`, re-read for the representatives' quality strings when an
    // intermediate is rendered; the clustering stage does not keep them. See
    // `crate::quals_for`.
    sorted_path: &std::path::Path,
) -> ParallelResult {
    let mut num_batches = nr_cores;
    let mut read_batches = batch_list(read_array, num_batches, bt);

    // Per batch: its own cluster map, representative map and (empty) database.
    let mut cluster_batches: Vec<OrderedClusters> = Vec::new();
    let mut rep_batches: Vec<FxHashMap<usize, ReadInfo>> = Vec::new();
    let mut db_batches: Vec<crate::cluster::MinimizerDatabase> = Vec::new();
    for batch in &read_batches {
        let mut c = OrderedClusters::default();
        let mut r = FxHashMap::default();
        for x in batch {
            c.insert(x.id, vec![x.id as u32]);
            r.insert(
                x.id,
                ReadInfo {
                    id: x.id,
                    batch_index: x.prev_batch_index,
                    acc: x.acc.clone(),
                    seq: x.seq.clone(),
                    err_per_base: x.err_per_base,
                    score: x.score,
                    error_rate: None,
                },
            );
        }
        cluster_batches.push(c);
        rep_batches.push(r);
        db_batches.push(crate::cluster::MinimizerDatabase::new());
    }

    let mut out = ParallelResult {
        clusters: OrderedClusters::default(),
        representatives: FxHashMap::default(),
        intermediates: Vec::new(),
        mapped_passed: 0,
        aln_passed: 0,
        aln_called: 0,
        skipped_short: 0,
        times: sweep::StageTimes::default(),
    };

    let mut it = 1usize;
    loop {
        // The one-batch case returns straight out of the loop, without a merge
        // iteration and without writing intermediates.
        if read_batches.len() == 1 {
            let res = sweep::reads_to_clusters(
                std::mem::take(&mut cluster_batches[0]),
                std::mem::take(&mut rep_batches[0]),
                &read_batches[0],
                std::mem::take(&mut db_batches[0]),
                1,
                table,
                p,
            );
            out.mapped_passed += res.mapped_passed;
            out.aln_passed += res.aln_passed;
            out.aln_called += res.aln_called;
            out.skipped_short += res.skipped_short;
            out.times.add(&res.times);
            out.clusters = res.clusters;
            out.representatives = res.representatives;
            return out;
        }

        // The pool. One thread per batch, matching the reference's
        // `Pool(processes=num_batches)`.
        //
        // This is behaviour-neutral, and the reason is worth stating rather than
        // hoping: each worker owns its cluster map, representative map, read
        // slice and minimizer database, and returns new ones -- nothing is
        // shared and nothing is mutated across batches. `map_async` preserves
        // result order and so does joining the handles in order, so the merge
        // below sees the same sequence whatever order the threads finish in.
        // Verified, not assumed: all 27 equivalence cases pass either way.
        let inputs: Vec<_> = (0..read_batches.len())
            .map(|i| {
                (
                    std::mem::take(&mut cluster_batches[i]),
                    std::mem::take(&mut rep_batches[i]),
                    &read_batches[i],
                    std::mem::take(&mut db_batches[i]),
                    i as i64 + 1,
                )
            })
            .collect();
        let results: Vec<sweep::SweepResult> = std::thread::scope(|scope| {
            let handles: Vec<_> = inputs
                .into_iter()
                .map(|(c, r, reads, db, idx)| {
                    scope.spawn(move || sweep::reads_to_clusters(c, r, reads, db, idx, table, p))
                })
                .collect();
            handles
                .into_iter()
                .map(|h| h.join().expect("a clustering worker panicked"))
                .collect()
        });
        for res in &results {
            out.mapped_passed += res.mapped_passed;
            out.aln_passed += res.aln_passed;
            out.aln_called += res.aln_called;
            out.skipped_short += res.skipped_short;
            out.times.add(&res.times);
        }

        // merge_dicts: later dicts win, but the batches are disjoint by read id.
        let mut all_clusters = OrderedClusters::default();
        let mut all_reps: FxHashMap<usize, ReadInfo> = FxHashMap::default();
        let mut dbs: HashMap<i64, crate::cluster::MinimizerDatabase> = HashMap::new();
        for res in results {
            for (id, accs) in res.clusters.iter() {
                all_clusters.insert(id, accs.clone());
            }
            all_reps.extend(res.representatives);
            dbs.insert(res.batch_index, res.db);
        }

        // The survivors, re-sorted by score descending.
        let mut survivors: Vec<SweepRead> = all_reps
            .values()
            .map(|r| SweepRead {
                id: r.id,
                prev_batch_index: r.batch_index,
                acc: r.acc.clone(),
                seq: r.seq.clone(),
                hp_error_rate: r.error_rate,
                err_per_base: r.err_per_base,
                score: r.score,
            })
            .collect();
        // `sorted(all_representatives.items(), key=lambda x: x[1][5], reverse=True)`
        // -- a dict of ints, so its own order is deterministic but not insertion
        // order. Ties are broken by id here to keep this reproducible; see the
        // note in PORTING.md.
        survivors.sort_by(|a, b| {
            b.score
                .partial_cmp(&a.score)
                .expect("scores are finite")
                .then(a.id.cmp(&b.id))
        });

        if num_batches == 1 {
            out.clusters = all_clusters;
            out.representatives = all_reps;
            return out;
        }
        out.intermediates.push(render_intermediate(
            &all_clusters,
            &all_reps,
            read_array,
            sorted_path,
        ));

        it += 1;
        let _ = it;
        read_batches = batch_list_merge_consecutive(&survivors);
        num_batches = read_batches.len();
        if num_batches == 0 {
            // The reference reaches Pool(processes=0) and dies with
            // "ValueError: Number of processes must be at least 1".
            out.clusters = all_clusters;
            out.representatives = all_reps;
            return out;
        }

        cluster_batches.clear();
        rep_batches.clear();
        db_batches.clear();
        for batch in &read_batches {
            let mut c = OrderedClusters::default();
            let mut r = FxHashMap::default();
            let lowest = batch
                .iter()
                .map(|x| x.prev_batch_index)
                .min()
                .expect("batches are non-empty here");
            for x in batch {
                if let Some(v) = all_clusters.map.get(&x.id) {
                    c.insert(x.id, v.clone());
                }
                if let Some(v) = all_reps.get(&x.id) {
                    r.insert(x.id, v.clone());
                }
            }
            cluster_batches.push(c);
            rep_batches.push(r);
            db_batches.push(dbs.remove(&lowest).unwrap_or_default());
        }
    }
}

/// `print_intermediate_results` -- the two files written per merge iteration.
///
/// Both sort by cluster size only, with ties falling through to insertion order.
/// `pre_clusters.csv` uses the **internal** cluster id, not a renumbered one, and
/// `cluster_origins.csv` writes the accession **with** its score suffix -- unlike
/// `final_cluster_origins.tsv`, which strips it.
fn render_intermediate(
    clusters: &OrderedClusters,
    reps: &FxHashMap<usize, ReadInfo>,
    reads: &[SweepRead],
    sorted_path: &std::path::Path,
) -> (String, String) {
    // Quality strings for this pass's representatives only.
    let want: rustc_hash::FxHashSet<usize> = reps.keys().copied().collect();
    let quals = crate::quals_for(sorted_path, &want).unwrap_or_default();
    let mut order: Vec<usize> = clusters.order.clone();
    order.sort_by(|a, b| clusters.map[b].len().cmp(&clusters.map[a].len()));

    let mut pre = String::new();
    let mut origins = String::new();
    for c_id in &order {
        for id in &clusters.map[c_id] {
            pre.push_str(&format!(
                "{}\t{}\n",
                c_id,
                crate::strip_score(&reads[*id as usize].acc)
            ));
        }
    }
    for c_id in &order {
        let r = &reps[c_id];
        origins.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            r.id,
            r.acc,
            String::from_utf8_lossy(&r.seq.to_bytes()),
            quals.get(c_id).map(String::as_str).unwrap_or(""),
            crate::pyfloat::repr(r.score),
            crate::pyfloat::repr(r.error_rate.unwrap_or(f64::NAN)),
        ));
    }
    (pre, origins)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn r(id: usize, b: i64, len: usize, score: f64) -> SweepRead {
        SweepRead {
            id,
            prev_batch_index: b,
            acc: format!("r{id}_{score}").into(),
            seq: crate::packed::PackedSeq::from_bytes(&vec![b'A'; len]).0,
            hp_error_rate: None,
            err_per_base: 0.0,
            score,
        }
    }

    #[test]
    fn batch_type_weighted_is_not_implemented() {
        // The CLI documents it; no branch handles it. Finding 15.
        assert!(BatchType::parse("weighted").is_none());
        assert!(BatchType::parse("total_nt").is_some());
        assert!(BatchType::parse("nr_reads").is_some());
        assert!(BatchType::parse("read_lengths_squared").is_some());
    }

    #[test]
    fn nr_reads_uses_chunk_size_len_over_cores_plus_one() {
        let reads: Vec<SweepRead> = (0..10).map(|i| r(i, 0, 100, 1.0)).collect();
        // 10/4 + 1 = 3, so 3+3+3+1
        let b = batch_list(&reads, 4, BatchType::NrReads);
        assert_eq!(
            b.iter().map(|x| x.len()).collect::<Vec<_>>(),
            vec![3, 3, 3, 1]
        );
    }

    #[test]
    fn total_nt_closes_a_batch_when_the_running_length_reaches_the_chunk() {
        let reads: Vec<SweepRead> = (0..8).map(|i| r(i, 0, 100, 1.0)).collect();
        // total 800, /4 = 200, +1 = 201; so batches close after 3 reads (300)
        let b = batch_list(&reads, 4, BatchType::TotalNt);
        assert_eq!(b.iter().map(|x| x.len()).collect::<Vec<_>>(), vec![3, 3, 2]);
    }

    #[test]
    fn merge_consecutive_closes_on_a_rising_threshold() {
        // score order, batch indices crossing 2 then 4
        let reads = vec![
            r(0, 1, 10, 9.0),
            r(1, 2, 10, 8.0),
            r(2, 3, 10, 7.0),
            r(3, 4, 10, 6.0),
        ];
        let b = batch_list_merge_consecutive(&reads);
        // ids 0,1 have index <= 2; then 2 opens a new batch (threshold -> 4),
        // and 3 has index 4 <= 4 so it joins it.
        assert_eq!(b.len(), 2);
        assert_eq!(b[0].iter().map(|x| x.id).collect::<Vec<_>>(), vec![0, 1]);
        assert_eq!(b[1].iter().map(|x| x.id).collect::<Vec<_>>(), vec![2, 3]);
    }
}
