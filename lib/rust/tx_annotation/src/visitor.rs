//!
//! This module defines a visitor strategy useful for configuring the `Aligner` and
//! performing operations at different levels including the computation of metrics.
//!
//! The stategy is inspired by the visitor implementation within the [rust compiler](https://github.com/rust-lang/rust/blob/c9bacb70f0b19d324a548bd7942692ab18d159a4/src/librustc/hir/intravisit.rs)
//! You can read more about the pattern [here](https://github.com/rust-unofficial/patterns/blob/master/patterns/visitor.md)
//!
//! It is recommended that you gain familiarity with the trait definied here before looking into
//! how various metrics are computed in the `ALIGN_AND_COUNT` stage.
//!
//! ## Barcode Annotation Graph
//! The key concept here is to think about all the aligned and annotated read data from a barcode
//! as forming a graph. We will call this graph a `Barcode Annotation Graph`. We are only
//! considering reads from the gene expression library for this trait.
//!
//! Consider the following annotation graph for a barcode which has a total of `k+1` reads:
//!
//! ```text
//!                     +------------------------+
//!                     |Barcode CCGGCCTCAAGGCGTC|
//!                     +------------+-----------+
//!                                  |
//!                                  |
//!   +---------------------+--------+------------------------+
//!   |                     |                                 |
//!   |                     |                                 |
//! +-+--+                +-+--+                            +-+--+
//! |Read|                |Read|       ................     |Read|
//! |[0] |                |[1] |                            |[k] |
//! +-+--+                +-+--+                            +----+
//!   |                     |                                 .
//!   +>Genomes mapped      +>Unmapped                        .
//!   +>Genes mapped                                          .
//!   +>Regions mapped                                        .
//!   +>Genome conf mapped
//!   +>Gene conf mapped
//!   +>Region conf mapped
//!   +>Duplicate Info
//!
//! ```
//!
//! ### Properties of the graph
//! - Barcode node is the root of the graph
//! - Each read (which is aligned and annotated) node is a child of the barcode node
//! - Each read node could have multiple child nodes depending on the annotation of that read.
//! - Invariants associated with the annotation
//!   - Every mapped read will:
//!     - Map to at least one genome in the reference.
//!     - Map to one of the `AnnotationRegion`s (exonic/intronic/intergenic)
//!     - Map to zero or more genes
//!   - Every confidently mapped read will:
//!     - Map to exactly one genome
//!     - Map to exactly one of the `AnnotationRegion`s (exonic/intronic/intergenic)
//!     - Map to zero or one gene
//!
//! ### Graph operations
//! There are two types of operations the are defined on the graph
//! 1. **Visit** a node: Every method in the `AnnotatedReadVisitor` trait is a
//!    visit to one of the nodes in the graph (`visit_*`). Each function is a hook to potentially
//!    override what happens when you visit that node. Each function's default implementation
//!    recursively visits the children via the corresponding `walk_` method. e.g., the
//!    `visit_barcode_read_annotations` method by default calls `walk_barcode_read_annotations`.
//!    Each overridden visit method has full control over what happens with its node, it can do
//!    its own traversal of the node's children, call `walk_*` to apply the default traversal
//!    algorithm, or prevent deeper traversal by doing nothing.
//! 2. **Walk** from a node: This corresponds to a traversal of the node's children. For the
//!    graph described above there are walks described at two levels: the barcode level and the
//!    read level. The `walk_` function at the barcode level visits each of the read nodes (by
//!    calling various `visit_*` functions) and the `walk_` function at the read level visits
//!    various annotation properties of the read.

use crate::mark_dups::DupInfo;
use crate::read::{AnnotationInfo, ReadAnnotations};
use crate::transcript::AnnotationRegion;
use cr_types::rna_read::RnaRead;
use cr_types::GenomeName;
use std::collections::HashSet;
use transcriptome::Gene;

/// An abstract visitor over the set of aligned and annotated reads from a barcode. For an
/// overview of the strategy, read the comments at the top of this module. This visitor
/// separates the details of the alignment/annotation algorithm from what we would like to
/// compute from the annotated reads (metrics, for example).
///
/// Every function in this trait is of the form `visit_*` which visits one of the nodes in the
/// barcode annotation graph (see comments at the top). The default implementation will `walk`
/// over all the children of the node. The implementor can override the default behavior.
///
/// TODO
/// - Visit UMI count reads
/// - Should we abstract over mutability of annotations so that the visitor can mutate the
///   annotations? This might be a cleaner strategy to tag the BAM records for example.
/// - Right now this is made for single end reads. Need to upgrade this for paired end.
pub trait AnnotatedReadVisitor: Sized {
    /// Visit one of the nodes under one barcode. The caller needs to guarantee that the visit
    /// happens in barcode sorted order.
    ///
    /// IMPORTANT: THIS FUNCTION IS THE ONLY ENTRY POINT TO THE GRAPH THAT THE `Aligner` USES
    ///
    /// By default, this function recursively visits all the child nodes of the read node by
    /// calling the `walk` function.
    /// Be mindful of the traversal needed when you are overriding this method.
    fn visit_read_annotation(&mut self, annotation: &ReadAnnotations) {
        walk_read_annotation(self, annotation);
    }

    /// Visit the set of genomes of a mapped read. A mapped read will map to at least one
    /// genome.
    ///
    /// In the default implementation, this function is called exactly once for every **mapped** read.
    /// See `walk_read_annotations` for the traversal logic.
    fn visit_mapped_read_genomes(
        &mut self,
        _annotation: &ReadAnnotations,
        _genomes: HashSet<&GenomeName>,
    ) {
    }

    /// Visit the set of (genome, gene) of a mapped read. It is possible that the set is empty,
    /// for example, a read could be intergenic.
    ///
    /// In the default implementation, this function is called exactly once for every **mapped** read.
    /// See `walk_read_annotations` for the traversal logic.
    fn visit_mapped_read_genes(
        &mut self,
        _annotation: &ReadAnnotations,
        _genes: HashSet<(&GenomeName, &Gene)>,
    ) {
    }

    /// Visit the set of (genome, region) of a mapped read. The region is one of exonic, intronic
    /// or intergenic. The set is guaranteed to be non empty.
    ///
    /// In the default implementation, this function is called **exactly once for every mapped
    /// read**.See `walk_read_annotations` for the traversal logic.
    fn visit_mapped_read_regions(
        &mut self,
        _annotation: &ReadAnnotations,
        _regions: HashSet<(&GenomeName, AnnotationRegion)>,
    ) {
    }

    /// Visit the genome corresponding to a confidently mapped read.
    ///
    /// In the default implementation, this function is called **exactly once for every
    /// confidently mapped read**. See `walk_read_annotations` for the traversal logic.
    fn visit_conf_mapped_read_genome(
        &mut self,
        _annotation: &ReadAnnotations,
        _genome: &GenomeName,
    ) {
    }

    /// Visit the (genome, gene) corresponding to a confidently mapped read, if it exists.
    ///
    /// In the default implementation, this function is called **exactly once for every
    /// confidently mapped read that confidently maps to a gene**. See `walk_read_annotations`
    /// for the traversal logic.
    fn visit_conf_mapped_read_gene(
        &mut self,
        _annotation: &ReadAnnotations,
        _genome: &GenomeName,
        _gene: &Gene,
    ) {
    }

    /// Visit the (genome, region) corresponding to a confidently mapped read.
    ///
    /// In the default implementation, this function is called **exactly once for every
    /// confidently mapped read**. See `walk_read_annotations` for the traversal logic.
    fn visit_conf_mapped_read_region(
        &mut self,
        _annotation: &ReadAnnotations,
        _genome: &GenomeName,
        _region: AnnotationRegion,
    ) {
    }

    /// Visit a read that is confidently mapped antisense.
    ///
    /// In the default implementation, this function is called **exactly once for every
    /// confidently mapped read that is mapped antisense**. See `walk_read_annotations`
    /// for the traversal logic.
    fn visit_antisense_conf_mapped_read(&mut self, _annotation: &ReadAnnotations) {}

    /// In the default implementation, this function is called exactly once for each read that
    /// comes in the aligner stage.
    fn visit_every_read(&mut self, _read: &RnaRead) {}

    /// In the default implementation, this function is called exactly once for every unmapped read
    /// processed. See `walk_read_annotations` for the traversal logic.
    fn visit_unmapped_read(&mut self, _annotation: &ReadAnnotations) {}

    /// In the default implementation, this function is called exactly once for every mapped read
    /// processed. See `walk_read_annotations` for the traversal logic.
    fn visit_mapped_read(&mut self, _annotation: &ReadAnnotations) {}

    /// In the implementation, this function is called exactly once for every feature read
    /// processed. See `walk_read_annotations` for the traversal logic.
    fn visit_feature_read(&mut self, _annotation: &ReadAnnotations) {}

    /// Visit the duplicate info associated with a read. Duplicate marking is only done for
    /// reads with a valid barcode.
    ///
    /// In the default implementation, this function is called exactly once for every read that
    /// confidently maps to transcriptome and has a valid barcode.
    /// See `walk_read_annotations` for the traversal logic.
    fn visit_dup_info(&mut self, _annotation: &ReadAnnotations, _dup_info: &DupInfo) {}
}

/// The default traversal from the read node in the barcode annotation graph that visits various
/// annotation property nodes of the read.
pub fn walk_read_annotation<V: AnnotatedReadVisitor>(
    visitor: &mut V,
    annotation: &ReadAnnotations,
) {
    visitor.visit_every_read(&annotation.read);

    if annotation.is_mapped() {
        visitor.visit_mapped_read(annotation);
        visitor.visit_mapped_read_genomes(annotation, annotation.mapped_genomes());
        visitor.visit_mapped_read_genes(annotation, annotation.mapped_genes());
        visitor.visit_mapped_read_regions(annotation, annotation.mapped_regions());

        if let Some(genome) = annotation.conf_mapped_genome() {
            visitor.visit_conf_mapped_read_genome(annotation, genome);
        }

        if let Some((genome, gene)) = annotation.conf_mapped_gene() {
            visitor.visit_conf_mapped_read_gene(annotation, genome, gene);
        }

        if let Some((genome, region)) = annotation.conf_mapped_region() {
            visitor.visit_conf_mapped_read_region(annotation, genome, region);
        }

        if annotation.is_conf_mapped_antisense() {
            visitor.visit_antisense_conf_mapped_read(annotation);
        }
    } else {
        visitor.visit_unmapped_read(annotation);
    }
    // Dup info is only computed for reads with a valid barcode that confidently
    // maps to a transcriptome
    if let Some(ref dup_info) = annotation.dup_info {
        visitor.visit_dup_info(annotation, dup_info);
    };
    // Only accumulate feature info for feature-barcode reads
    if annotation.is_feature_read() {
        visitor.visit_feature_read(annotation);
    }
}
