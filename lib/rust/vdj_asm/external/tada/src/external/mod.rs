//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

// use load_graph;
// use bwt;
// use std::ffi::{CString, CStr};
// use std::ptr;
// use libc;
// use bitenc;
// use std::slice;
// use std::collections::HashSet;
// use std::mem;
//
// #[repr(C)]
// pub struct Array {
// data: *const libc::c_void,
// len: libc::size_t,
// }
//
// impl Array {
// unsafe fn as_u32_slice(&self) -> &[u32] {
// assert!(!self.data.is_null());
// slice::from_raw_parts(self.data as *const u32, self.len as usize)
// }
//
// fn from_vec<T>(mut vec: Vec<T>) -> Array {
// Important to make length and capacity match
// A better solution is to track both length and capacity
// vec.shrink_to_fit();
// let array = Array { data: vec.as_ptr() as *const libc::c_void, len: vec.len() as libc::size_t };
// Take ownership of vec and prevent it from being destroyed.
// mem::forget(vec);
//
// array
// }
// }
//
// #[no_mangle]
// pub extern fn sa_step() -> u32 {
// load_graph::SA_STEP as u32
// }
//
// #[no_mangle]
// pub extern fn bwt_bucket_pos() -> u32 {
// load_graph::BWT_BUCKET_POS as u32
// }
//
// Reads an assembly graph from a binary file and creates a Graph object.
//
// # Arguments
//
// * `name` - filename
// #[no_mangle]
// pub extern fn make_graph(name: *const libc::c_char) -> *mut libc::c_void {
// unsafe {
// This ptr to string convertion is unsafe.
// let filename: String = CStr::from_ptr(name).to_string_lossy().into_owned();
// let g = match load_graph::load_graph(&filename) {
// Err(err) => panic!(err),
// Ok(g) => g,
// };
// let out_g = Box::new(g);
// Box::into_raw(out_g) as *mut libc::c_void
// }
// }
//
// #[no_mangle]
// pub extern fn delete_graph(g: *const libc::c_void) {
// unsafe {
// drop is technically not necessary. The graph will be dropped anyway
// because it will go out of scope.
// drop(Box::from_raw(g as *mut load_graph::Graph));
// }
// }
//
// #[no_mangle]
// pub extern fn write_graph(g: *const libc::c_void, name: *const libc::c_char) {
// unsafe {
// let filename: String = CStr::from_ptr(name).to_string_lossy().into_owned();
// let _ = load_graph::write_graph(&*(g as *const load_graph::Graph), &filename);
// }
// }
//
// #[no_mangle]
// pub extern fn read_graph(name: *const libc::c_char) -> *mut libc::c_void {
// unsafe {
// let filename: String = CStr::from_ptr(name).to_string_lossy().into_owned();
// let g = match load_graph::read_graph(&filename) {
// Err(_) => return ptr::null_mut() as *mut libc::c_void,
// Ok(g) => g,
// };
// let out_g = Box::new(g);
// Box::into_raw(out_g) as *mut libc::c_void
// }
// }
//
// #[no_mangle]
// pub extern fn num_edges(g : *const libc::c_void) -> u32 {
// unsafe {
// let graph : &load_graph::Graph = &*(g as *const load_graph::Graph);
// graph.num_edges() as u32
// }
// }
//
// #[no_mangle]
// pub extern fn edge_seq(g : *const libc::c_void, edge_idx: libc::c_uint) -> *const libc::c_char {
// unsafe {
// let graph : &load_graph::Graph = &*(g as *const load_graph::Graph);
// let s = graph.edge_seq(edge_idx as usize);
// let res = CString::new(s).unwrap();
// res.into_raw()
// }
// }
//
// #[no_mangle]
// pub extern fn merge_bwt_buckets(filenames: *const *const libc::c_char,
// nfiles: libc::c_uint, sa_step: libc::c_uint, outfile: *const libc::c_char) {
// unsafe {
// let values = slice::from_raw_parts(filenames, nfiles as usize);
// let names: Vec<String> = values.iter()
// .map(|&p| CStr::from_ptr(p))
// .map(|cs| cs.to_string_lossy().into_owned())
// .collect();
// let b = bwt::merge_bwt_buckets(&names, sa_step as usize);
// let outf: String = CStr::from_ptr(outfile).to_string_lossy().into_owned();
// bwt::write_bwt(&b, &outf);
// }
// }
//
// Compute a fraction of the BWT and write it to a file.
// #[no_mangle]
// pub extern fn compute_bwt_bucket(g : *const libc::c_void, val: libc::c_uint, filename: *const libc::c_char) {
// unsafe {
// let f: String = CStr::from_ptr(filename).to_string_lossy().into_owned();
// let graph : &load_graph::Graph = &*(g as *const load_graph::Graph);
// let b: bwt::BWT = bwt::compute_bwt_bucket(&graph.all_seq, load_graph::BWT_BUCKET_POS, val as usize);
// bwt::write_bwt(&b, &f);
// }
// }
//
// Compute a fraction of the BWT and write it to a file.
// #[no_mangle]
// pub extern fn compute_bwt(g : *const libc::c_void, filename: *const libc::c_char) {
// unsafe {
// let f: String = CStr::from_ptr(filename).to_string_lossy().into_owned();
// let graph : &load_graph::Graph = &*(g as *const load_graph::Graph);
// let (sa, bwt) = bwt::compute_bwt(&graph.all_seq, load_graph::BWT_BUCKET_POS, load_graph::SA_STEP);
// let b = bwt::BWT {bwt: bwt, sa: sa};
// bwt::write_bwt(&b, &f);
// }
// }
//
// #[no_mangle]
// pub extern fn create_indexed_graph(g: *const libc::c_void, filename: *const libc::c_char) -> *mut libc::c_void {
// unsafe {
// let f: String = CStr::from_ptr(filename).to_string_lossy().into_owned();
// let graph : &load_graph::Graph = &*(g as *const load_graph::Graph);
// let new_graph = load_graph::Graph::from_graph_and_index(graph, &f);
// let out_g = Box::new(new_graph);
// Box::into_raw(out_g) as *mut libc::c_void
// }
// }
//
// #[no_mangle]
// pub extern fn graph_contains_seq(g: *const libc::c_void, seq: *const libc::c_char) -> libc::c_uchar {
// unsafe {
// let graph : &load_graph::Graph = &*(g as *const load_graph::Graph);
// let dna: String = CStr::from_ptr(seq).to_string_lossy().into_owned();
// let bit_seq = bitenc::BitEnc::from_dna_string(&dna);
// graph.contains_seq(&bit_seq) as libc::c_uchar
// }
// }
//
// #[no_mangle]
// pub extern fn find_graph_edges(g: *const libc::c_void, seq: *const libc::c_char) -> Array {
// unsafe {
// let graph : &load_graph::Graph = &*(g as *const load_graph::Graph);
// let dna: String = CStr::from_ptr(seq).to_string_lossy().into_owned();
// let bit_seq = bitenc::BitEnc::from_dna_string(&dna);
// let edges = graph.find_edges(&bit_seq)
// .iter().cloned()
// .map(|x| x as u32)
// .collect::<Vec<u32>>();
// Array::from_vec(edges)
// }
// }
//
// struct GraphPaths {
// paths: Vec<load_graph::GraphPath>
// }
//
// impl GraphPaths {
// pub fn empty() -> Self {
// GraphPaths { paths: Vec::new() }
// }
// }
//
// #[no_mangle]
// pub extern fn num_paths(p: *const libc::c_void) -> u32 {
// unsafe {
// let paths : &GraphPaths = &*(p as *const GraphPaths);
// paths.paths.len() as u32
// }
// }
//
// #[no_mangle]
// pub extern fn path_seq(g: *const libc::c_void, p: *const libc::c_void,
// path_idx: *const libc::c_uint) -> *const libc::c_char {
// unsafe {
// let graph : &load_graph::Graph = &*(g as *const load_graph::Graph);
// let paths : &GraphPaths = &*(p as *const GraphPaths);
// let s = paths.paths[path_idx as usize].seq(graph);
// let res = CString::new(s).unwrap();
// res.into_raw()
// }
// }
//
// #[no_mangle]
// pub extern fn path_edges(p: *const libc::c_void, path_idx: *const libc::c_uint) -> Array {
// unsafe {
// let paths : &GraphPaths = &*(p as *const GraphPaths);
// let edges = paths.paths[path_idx as usize].edges.clone();
// Array::from_vec(edges)
// }
// }
//
// #[no_mangle]
// pub extern fn find_graph_paths(g: *const libc::c_void, start_edge: libc::c_uint,
// end_edge: libc::c_uint, max_hops: libc::c_uint, good_edges: *const libc::c_uint,
// nedges: libc::c_uint) -> *const libc::c_void {
//
// unsafe {
// let graph : &load_graph::Graph = &*(g as *const load_graph::Graph);
// if start_edge as usize > graph.num_edges() {
// println!("Invalid starting edge idx.");
// return Box::into_raw(Box::new(GraphPaths::empty())) as *const libc::c_void;
// }
// if end_edge as usize > graph.num_edges() {
// println!("Invalid ending edge idx.");
// return Box::into_raw(Box::new(GraphPaths::empty())) as *const libc::c_void;
// }
// let values = slice::from_raw_parts(good_edges, nedges as usize);
// let good_edge_set: HashSet<usize> = values.iter()
// .map(|&p| p as usize)
// .collect();
// let paths = GraphPaths{ paths: graph.find_paths(start_edge as usize, end_edge as usize, max_hops as usize,
// &good_edge_set)};
// Box::into_raw(Box::new(paths)) as *const libc::c_void
// }
// }
//
// #[no_mangle]
// pub extern fn graph_stats(g: *const libc::c_void,
// summary_json: *const libc::c_char, edge_len_tsv: *const libc::c_char,
// seqs_file: *const libc::c_char, gt_hits_tsv: *const libc::c_char) {
// unsafe {
// let graph : &load_graph::Graph = &*(g as *const load_graph::Graph);
// let f1: String = CStr::from_ptr(summary_json).to_string_lossy().into_owned();
// let f2: String = CStr::from_ptr(edge_len_tsv).to_string_lossy().into_owned();
// let f3: String = CStr::from_ptr(seqs_file).to_string_lossy().into_owned();
// let f4: String = CStr::from_ptr(gt_hits_tsv).to_string_lossy().into_owned();
//
// graph.graph_stats(&f1, &f2, &f3, &f4);
// }
//
// }
//
