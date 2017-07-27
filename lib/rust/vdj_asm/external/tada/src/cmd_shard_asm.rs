//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

//! Command for assembling MSP shard Kmers into short DeBruijn graph snippets
//!
//! Enumerates all Kmers observed in the MSP shard and the barcodes observed on each.  Filters
//! the kmers by an arbitrary rule.  Then labels each Kmer with all the observed extensions from it.
//! Note, the extension in this phase will be a super-set of to Kmers that pass the filtering.
//! This will be resolved later.

use std::path::{Path, PathBuf};
use rustc_serialize::json::Json;
use std::collections::BTreeMap;
use std::fs;
use std::env;

use martian::{self, JsonDict, MartianStage};
use utils;
use debruijn;
use shard_buf::PodReader;
use cmd_msp;

pub fn main_shard_asm(shard_chunks: Vec<PathBuf>, sedge_asm_out: &Path, sedge_bcs_out: &Path) {

    let mut bsps = Vec::new();
    for p in shard_chunks {
        let mut reader = PodReader::new(&p);
        let mut read_bsps = reader.read();
        bsps.append(&mut read_bsps);
    }

    info!("Generating kmers from {} bsps", bsps.len());
    let (valid_kmers, all_kmers) = utils::process_kmer_shard(&bsps);
    info!("Got valid kmers: {}", valid_kmers.len());

    info!("Getting kmer extensions");
    let kmer_extensions = debruijn::kmer_extensions(&valid_kmers, &all_kmers, &bsps, true);

    info!("Building sedges");
    let sedges = debruijn::build_sedges(&kmer_extensions, true);
    info!("Got sedges: {}", sedges.len());

    // BCs covered by each sedge
    let mut bc_vec = Vec::new();
    let mut sedge_bcs = utils::MultiVec::new();

    info!("Getting barcodes for sedges");
    for &(sedge, _) in sedges.iter() {
        debruijn::barcodes_for_sedge(&mut bc_vec, sedge, &valid_kmers);
        sedge_bcs.add(bc_vec.drain(..));
    }

    info!("Writing sedges: {:?}", sedge_asm_out);
    utils::write_obj(&sedges, sedge_asm_out).expect("write failed");

    info!("Writing sedge bcs: {:?}", sedge_bcs_out);
    utils::write_obj(&sedge_bcs, sedge_bcs_out).expect("write failed");
}


pub struct ShardAsmMartian;

impl MartianStage for ShardAsmMartian {
    fn split(&self, args: JsonDict) -> JsonDict {

        let _msp_chunks= args["chunks"].as_array().unwrap().clone();
        let mut msp_chunks : Vec<cmd_msp::MspChunk> = Vec::new();
        for chnk in _msp_chunks
        {
            let cc = martian::json_decode(chnk);
            msp_chunks.push(cc);
        }

        let nshards = msp_chunks[0].msp.len();

        let mut shards = Vec::new();
        for _ in 0 .. nshards {
            shards.push(Vec::new());
        }

        for msp_chunk in msp_chunks {
            for (i, shard_chunk) in msp_chunk.msp.into_iter().enumerate() {
                shards[i].push(Json::String(shard_chunk));
            }
        }

        let mut chunks : Vec<Json> = Vec::new();

        for shard in shards
        {
            let mut bt = BTreeMap::new();
            bt.insert("chunk".to_string(), Json::Array(shard));
            bt.insert("__mem_gb".to_string(), Json::F64(3.0));
            chunks.push(Json::Object(bt));
        }

        let chunk_array = Json::Array(chunks);
        let mut cc =  BTreeMap::new();
        cc.insert("chunks".to_string(), chunk_array);
        cc
    }

    fn main(&self, args: JsonDict, outs: JsonDict) -> JsonDict {
        let shard_chunks = args["chunk"].as_array().unwrap().iter().map(|v| PathBuf::from(v.as_string().unwrap())).collect();

        let sedge_asm = Path::new({outs.get("sedge_asm").unwrap().as_string().unwrap().clone()});
        let sedge_bcs = Path::new({outs.get("sedge_bcs").unwrap().as_string().unwrap().clone()});

        main_shard_asm(shard_chunks, sedge_asm, sedge_bcs);
        outs.clone()
    }

    fn join(&self, _: JsonDict, _: JsonDict, _: Vec<JsonDict>, chunk_outs: Vec<JsonDict>) -> JsonDict {

            let mut sedge_asms: Vec<Json> = Vec::new();
            let mut sedge_bcss : Vec<Json> = Vec::new();

            let pwd = env::current_dir().unwrap();

            for (idx, c) in chunk_outs.into_iter().enumerate() {

                //let _old_name =
                let old_name = c.get("sedge_asm").unwrap().as_string().unwrap();
                let mut new_name = pwd.clone();
                new_name.push(format!("chunk{}.sedge_asm", idx));
                fs::rename(old_name, new_name.clone());
                sedge_asms.push(Json::String(new_name.to_str().unwrap().to_string()));

                let old_name = c.get("sedge_bcs").unwrap().as_string().unwrap();
                let mut new_name = pwd.clone();
                new_name.push(format!("chunk{}.sedge_bcs", idx));
                fs::rename(old_name, new_name.clone());
                sedge_bcss.push(Json::String(new_name.to_str().unwrap().to_string()));
            }

            let mut obj = BTreeMap::new();
            obj.insert("sedge_asm".to_string(), Json::Array(sedge_asms));
            obj.insert("sedge_bcs".to_string(), Json::Array(sedge_bcss));
            obj
    }
}
