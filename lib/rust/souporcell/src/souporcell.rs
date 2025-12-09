//! Souporcell source code
//!
//! MIT License
//!
//! Copyright (c) 2019 Haynes Heaton
//!
//! Permission is hereby granted, free of charge, to any person obtaining a copy
//! of this software and associated documentation files (the "Software"), to deal
//! in the Software without restriction, including without limitation the rights
//! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//! copies of the Software, and to permit persons to whom the Software is
//! furnished to do so, subject to the following conditions:
//!
//! The above copyright notice and this permission notice shall be included in all
//! copies or substantial portions of the Software.
use crate::{ProbMatrix, SouporcellChunkResult};
use anyhow::Result;
use itertools::izip;
use metric::{TxHashMap, TxHashSet};
use rand::rngs::SmallRng;
use rand::Rng;
use rust_htslib::bcf::{Read, Reader as VCFReader};
use std::f32;
use std::fs::File;
use std::io::{BufRead, BufReader};

pub(super) fn expectation_maximization(
    loci: usize,
    mut cluster_centers: ProbMatrix,
    cell_data: &[CellData],
    params: &Params,
    epoch: usize,
    thread_num: usize,
) -> SouporcellChunkResult {
    let mut sums: ProbMatrix = Vec::new();
    let mut denoms: ProbMatrix = Vec::new();
    for cluster in 0..params.num_clusters {
        sums.push(Vec::new());
        denoms.push(Vec::new());
        for _index in 0..loci {
            sums[cluster].push(1.0);
            denoms[cluster].push(2.0); // psuedocounts
        }
    }

    let log_prior: f32 = (1.0 / (params.num_clusters as f32)).ln();

    let mut _change = 1000.0;
    let mut iterations = 0;
    let mut total_log_loss = f32::NEG_INFINITY;
    let mut _total_log_loss_binom = f32::NEG_INFINITY;
    let mut final_log_probabilities = Vec::new();
    for _cell in 0..cell_data.len() {
        final_log_probabilities.push(Vec::new());
    }
    let log_loss_change_limit = 0.01 * (cell_data.len() as f32);
    let temp_steps = 9;
    let mut last_log_loss = f32::NEG_INFINITY;
    for temp_step in 0..temp_steps {
        let mut log_loss_change = 10000.0;
        while log_loss_change > log_loss_change_limit && iterations < 1000 {
            let mut log_binom_loss = 0.0;
            reset_sums_denoms(
                loci,
                &mut sums,
                &mut denoms,
                &cluster_centers,
                params.num_clusters,
            );
            for (celldex, cell) in cell_data.iter().enumerate() {
                let log_binoms = binomial_loss(cell, &cluster_centers, log_prior, celldex);
                log_binom_loss += log_sum_exp(&log_binoms);
                let mut temp =
                    (cell.total_alleles / (20.0 * 2.0f32.powf(temp_step as f32))).max(1.0);
                if temp_step == temp_steps - 1 {
                    temp = 1.0;
                }
                let probabilities = normalize_in_log_with_temp(&log_binoms, temp);
                update_centers_average(&mut sums, &mut denoms, cell, &probabilities);

                final_log_probabilities[celldex] = log_binoms; //log_probabilities;
            }

            total_log_loss = log_binom_loss;
            log_loss_change = log_binom_loss - last_log_loss; //log_loss - last_log_loss;
            last_log_loss = log_binom_loss; //log_loss;

            update_final(loci, &sums, &denoms, &mut cluster_centers);
            iterations += 1;
            eprintln!(
                "binomial\t{thread_num}\t{epoch}\t{iterations}\t{temp_step}\t{log_binom_loss}\t{log_loss_change}"
            );
        }
    }

    SouporcellChunkResult {
        log_loss: total_log_loss,
        log_probs: final_log_probabilities,
    }
}

#[allow(unused_variables)]
fn binomial_loss(
    cell_data: &CellData,
    cluster_centers: &[Vec<f32>],
    log_prior: f32,
    _cellnum: usize,
) -> Vec<f32> {
    let mut log_probabilities: Vec<f32> = Vec::new();
    let mut sum = 0.0;
    for (cluster, center) in cluster_centers.iter().enumerate() {
        log_probabilities.push(log_prior);
        for (locus_index, locus) in cell_data.loci.iter().enumerate() {
            log_probabilities[cluster] += cell_data.log_binomial_coefficient[locus_index]
                + (cell_data.alt_counts[locus_index] as f32) * center[*locus].ln()
                + (cell_data.ref_counts[locus_index] as f32) * (1.0 - center[*locus]).ln();
        }
        sum += log_probabilities[cluster];
    }

    log_probabilities
}

fn log_sum_exp(p: &[f32]) -> f32 {
    let max_p: f32 = p.iter().copied().fold(f32::NEG_INFINITY, f32::max);
    let sum_rst: f32 = p.iter().map(|x| (x - max_p).exp()).sum();
    max_p + sum_rst.ln()
}

fn normalize_in_log_with_temp(log_probs: &[f32], temp: f32) -> Vec<f32> {
    let mut normalized_probabilities: Vec<f32> = Vec::new();
    let mut new_log_probs: Vec<f32> = Vec::new();
    for log_prob in log_probs {
        new_log_probs.push(log_prob / temp);
    }
    let sum = log_sum_exp(&new_log_probs);
    for prob in log_probs {
        normalized_probabilities.push((prob - sum).exp());
    }
    normalized_probabilities
}

fn update_final(
    loci: usize,
    sums: &[Vec<f32>],
    denoms: &[Vec<f32>],
    cluster_centers: &mut [Vec<f32>],
) {
    for locus in 0..loci {
        for cluster in 0..sums.len() {
            let update = sums[cluster][locus] / denoms[cluster][locus];
            cluster_centers[cluster][locus] = update.clamp(0.01, 0.99); //max(0.0001, min(0.9999, update));
        }
    }
}

fn reset_sums_denoms(
    loci: usize,
    sums: &mut [Vec<f32>],
    denoms: &mut [Vec<f32>],
    _cluster_centers: &[Vec<f32>],
    num_clusters: usize,
) {
    for cluster in 0..num_clusters {
        for index in 0..loci {
            sums[cluster][index] = 1.0;
            denoms[cluster][index] = 2.0;
        }
    }
}

fn update_centers_average(
    sums: &mut [Vec<f32>],
    denoms: &mut [Vec<f32>],
    cell: &CellData,
    probabilities: &[f32],
) {
    for locus in 0..cell.loci.len() {
        for (cluster, _probability) in probabilities.iter().enumerate() {
            sums[cluster][cell.loci[locus]] +=
                probabilities[cluster] * (cell.alt_counts[locus] as f32);
            denoms[cluster][cell.loci[locus]] +=
                probabilities[cluster] * ((cell.alt_counts[locus] + cell.ref_counts[locus]) as f32);
        }
    }
}

pub(super) fn init_cluster_centers(
    loci_used: usize,
    cell_data: &[CellData],
    params: &Params,
    rng: &mut SmallRng,
    locus_to_index: &TxHashMap<usize, usize>,
) -> Result<ProbMatrix> {
    if params.known_genotypes.is_some() {
        return init_cluster_centers_known_genotypes(loci_used, params, rng, locus_to_index);
    }
    match params.initialization_strategy {
        ClusterInit::RandomUniform => init_cluster_centers_uniform(loci_used, params, rng),
        ClusterInit::RandomAssignment => {
            init_cluster_centers_random_assignment(loci_used, cell_data, params, rng)
        }
    }
}

fn init_cluster_centers_known_genotypes(
    loci: usize,
    params: &Params,
    _rng: &mut SmallRng,
    locus_to_index: &TxHashMap<usize, usize>,
) -> Result<ProbMatrix> {
    let mut centers: ProbMatrix = Vec::new();
    for cluster in 0..params.num_clusters {
        centers.push(Vec::new());
        for _ in 0..loci {
            centers[cluster].push(0.5);
        }
    }
    let mut vcf_reader = VCFReader::from_path(params.known_genotypes.as_ref().unwrap())?;
    let known_sample_to_vcf_sample_idx: TxHashMap<String, usize> = vcf_reader
        .header()
        .samples()
        .iter()
        .enumerate()
        .map(|(i, s)| (String::from_utf8(s.to_vec()).unwrap(), i))
        .collect();
    let mut locus_id: usize = 0;
    for record in vcf_reader.records() {
        let record = record?;
        if let Some(loci_index) = locus_to_index.get(&locus_id) {
            if !params.known_genotypes_sample_names.is_empty() {
                for (sample_index, known_sample) in
                    params.known_genotypes_sample_names.iter().enumerate()
                {
                    let samples_genotypes = record.genotypes().unwrap();
                    let sample_vcf_idx = known_sample_to_vcf_sample_idx.get(known_sample).unwrap();
                    let sample_genotype = samples_genotypes.get(*sample_vcf_idx);
                    let alleles: Vec<_> = sample_genotype
                        .iter()
                        .filter_map(|x| x.index()) // Filter missing genotypes
                        .map(|x| x.max(1)) // Collapse all alt alleles to 1
                        .collect();
                    // Skip for this sample-locus if no alleles (SNP not present) or more than 2 alleles (triploid or more)
                    if alleles.is_empty() || alleles.len() > 2 {
                        continue;
                    }
                    let weight = alleles.iter().sum::<u32>() as f32 / (alleles.len() as f32);

                    centers[sample_index][*loci_index] = weight.clamp(0.01, 0.99);
                }
            } else {
                panic!("currently requiring known_genotypes_sample_names if known_genotypes set");
            }
        }
        locus_id += 1;
    }
    Ok(centers)
}

fn init_cluster_centers_uniform(
    loci: usize,
    params: &Params,
    rng: &mut SmallRng,
) -> Result<ProbMatrix> {
    let mut centers: ProbMatrix = Vec::new();
    for cluster in 0..params.num_clusters {
        centers.push(Vec::new());
        for _ in 0..loci {
            centers[cluster].push(rng.random::<f32>().clamp(0.0001, 0.9999));
        }
    }
    Ok(centers)
}

fn init_cluster_centers_random_assignment(
    loci: usize,
    cell_data: &[CellData],
    params: &Params,
    rng: &mut SmallRng,
) -> Result<ProbMatrix> {
    let mut sums: ProbMatrix = Vec::new();
    let mut denoms: ProbMatrix = Vec::new();
    for cluster in 0..params.num_clusters {
        sums.push(Vec::new());
        denoms.push(Vec::new());
        for _ in 0..loci {
            sums[cluster].push(rng.random::<f32>() * 0.01);
            denoms[cluster].push(0.01);
        }
    }
    for cell in cell_data {
        let cluster = rng.random_range(0..params.num_clusters);
        for locus in 0..cell.loci.len() {
            let alt_c = cell.alt_counts[locus] as f32;
            let total = alt_c + (cell.ref_counts[locus] as f32);
            let locus_index = cell.loci[locus];
            sums[cluster][locus_index] += alt_c;
            denoms[cluster][locus_index] += total;
        }
    }
    for cluster in 0..params.num_clusters {
        for locus in 0..loci {
            sums[cluster][locus] =
                sums[cluster][locus] / denoms[cluster][locus] + (rng.random::<f32>() / 2.0 - 0.25);
            sums[cluster][locus] = sums[cluster][locus].clamp(0.0001, 0.9999);
        }
    }

    Ok(sums)
}

pub(super) fn load_cell_data(
    params: &Params,
) -> (
    usize,
    usize,
    Vec<CellData>,
    Vec<usize>,
    TxHashMap<usize, usize>,
) {
    let alt_reader = File::open(&params.alt_mtx).expect("cannot open alt mtx file");

    let alt_reader = BufReader::new(alt_reader);
    let ref_reader = File::open(&params.ref_mtx).expect("cannot open ref mtx file");

    let ref_reader = BufReader::new(ref_reader);
    let mut used_loci: TxHashSet<usize> = TxHashSet::default();
    let mut total_loci = 0;
    let mut total_cells = 0;
    let mut all_loci: TxHashSet<usize> = TxHashSet::default();
    let mut locus_cell_counts: TxHashMap<usize, [u32; 2]> = TxHashMap::default();
    let mut locus_umi_counts: TxHashMap<usize, [u32; 2]> = TxHashMap::default();
    let mut locus_counts: TxHashMap<usize, TxHashMap<usize, [u32; 2]>> = TxHashMap::default();
    for (line_number, (alt_line, ref_line)) in
        izip!(alt_reader.lines(), ref_reader.lines()).enumerate()
    {
        let alt_line = alt_line.expect("cannot read alt mtx");
        let ref_line = ref_line.expect("cannot read ref mtx");
        if line_number > 2 {
            let alt_tokens: Vec<&str> = alt_line.split_whitespace().collect();
            let ref_tokens: Vec<&str> = ref_line.split_whitespace().collect();
            let locus = alt_tokens[0].to_string().parse::<usize>().unwrap() - 1;
            all_loci.insert(locus);
            let cell = alt_tokens[1].to_string().parse::<usize>().unwrap() - 1;
            let ref_count = ref_tokens[2].to_string().parse::<u32>().unwrap();
            let alt_count = alt_tokens[2].to_string().parse::<u32>().unwrap();
            assert!(locus < total_loci);
            assert!(cell < total_cells);
            let cell_counts = locus_cell_counts.entry(locus).or_insert([0; 2]);
            let umi_counts = locus_umi_counts.entry(locus).or_insert([0; 2]);
            if ref_count > 0 {
                cell_counts[0] += 1;
                umi_counts[0] += ref_count;
            }
            if alt_count > 0 {
                cell_counts[1] += 1;
                umi_counts[1] += alt_count;
            }
            let cell_counts = locus_counts.entry(locus).or_default();
            cell_counts.insert(cell, [ref_count, alt_count]);
        } else if line_number == 2 {
            let tokens: Vec<&str> = alt_line.split_whitespace().collect();
            total_loci = tokens[0].to_string().parse::<usize>().unwrap();
            total_cells = tokens[1].to_string().parse::<usize>().unwrap();
        }
    }
    let mut all_loci2: Vec<usize> = Vec::new();
    for loci in all_loci {
        all_loci2.push(loci);
    }
    let mut all_loci = all_loci2;

    all_loci.sort();
    let mut index_to_locus: Vec<usize> = Vec::new();
    let mut locus_to_index: TxHashMap<usize, usize> = TxHashMap::default();
    let mut cell_data: Vec<CellData> = Vec::new();
    for _cell in 0..total_cells {
        cell_data.push(CellData::new());
    }
    let mut locus_index = 0;
    for locus in all_loci {
        let cell_counts = locus_cell_counts.get(&locus).unwrap();
        let umi_counts = locus_umi_counts.get(&locus).unwrap();
        if cell_counts[0] >= params.min_ref
            && cell_counts[1] >= params.min_alt
            && umi_counts[0] >= params.min_ref_umis
            && umi_counts[1] >= params.min_alt_umis
        {
            used_loci.insert(locus);
            index_to_locus.push(locus);
            locus_to_index.insert(locus, locus_index);
            for (cell, counts) in locus_counts.get(&locus).unwrap() {
                if counts[0] + counts[1] == 0 {
                    continue;
                }
                cell_data[*cell].alt_counts.push(counts[1]);
                cell_data[*cell].ref_counts.push(counts[0]);
                cell_data[*cell].loci.push(locus_index);
                cell_data[*cell]
                    .allele_fractions
                    .push((counts[1] as f32) / ((counts[0] + counts[1]) as f32));
                cell_data[*cell].log_binomial_coefficient.push(
                    statrs::function::factorial::ln_binomial(
                        (counts[1] + counts[0]) as u64,
                        counts[1] as u64,
                    ) as f32,
                );
                cell_data[*cell].total_alleles += (counts[0] + counts[1]) as f32;
                //println!("cell {} locus {} alt {} ref {} fraction {}",*cell, locus_index, counts[1], counts[0],
                //    (counts[1] as f32)/((counts[0] + counts[1]) as f32));
            }
            locus_index += 1;
        }
    }
    eprintln!("total loci used {}", used_loci.len());

    (
        used_loci.len(),
        total_cells,
        cell_data,
        index_to_locus,
        locus_to_index,
    )
}

pub(super) struct CellData {
    allele_fractions: Vec<f32>,
    log_binomial_coefficient: Vec<f32>,
    alt_counts: Vec<u32>,
    ref_counts: Vec<u32>,
    loci: Vec<usize>,
    total_alleles: f32,
}

impl CellData {
    fn new() -> CellData {
        CellData {
            allele_fractions: Vec::new(),
            log_binomial_coefficient: Vec::new(),
            alt_counts: Vec::new(),
            ref_counts: Vec::new(),
            loci: Vec::new(),
            total_alleles: 0.0,
        }
    }
}

#[derive(Clone)]
pub(super) struct Params {
    pub(super) ref_mtx: String,
    pub(super) alt_mtx: String,
    pub(super) num_clusters: usize,
    pub(super) min_alt: u32,
    pub(super) min_ref: u32,
    pub(super) min_alt_umis: u32,
    pub(super) min_ref_umis: u32,
    pub(super) known_genotypes: Option<String>,
    pub(super) known_genotypes_sample_names: Vec<String>,
    pub(super) initialization_strategy: ClusterInit,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            ref_mtx: String::default(),
            alt_mtx: String::default(),
            num_clusters: 0,
            min_alt: 10,
            min_ref: 10,
            min_alt_umis: 0,
            min_ref_umis: 0,
            known_genotypes: None,
            known_genotypes_sample_names: Vec::new(),
            initialization_strategy: ClusterInit::RandomAssignment,
        }
    }
}

#[derive(Clone)]
#[allow(dead_code)]
pub(super) enum ClusterInit {
    RandomUniform,
    RandomAssignment,
}
