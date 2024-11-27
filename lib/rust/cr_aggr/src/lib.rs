mod errors;
pub mod match_vdj_outs;

/// Trim and merge molecules
pub mod merge_molecules;

pub mod create_antigen_clonotype_clustermap;
pub mod parse_aggr_csv;
pub mod process_vdj_proto;
pub mod setup_vdj_aggr;
pub mod websummary;
pub mod write_aggr_ann;
pub mod write_contig_proto;
pub mod write_ws_json;

#[cfg(test)]
#[ctor::ctor]
fn init() {
    // this ensures insta knows where to find its snap tests
    let cwd = std::env::current_dir().unwrap();
    let workspace_root = cwd.parent().unwrap();
    std::env::set_var("INSTA_WORKSPACE_ROOT", workspace_root);
}
