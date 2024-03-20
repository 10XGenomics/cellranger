use crate::reference::reference_info::ReferenceInfo;
use crate::types::GenomeName;
use anyhow::Result;
use std::path::Path;

/// Compute the genome for each chromosome. The output is a vector of genome names in the same
/// order as the names of chromosomes in `star/chrName.txt`. If we have more than 1 genome in the
/// reference, by convention, we store the chromosome names as `{genome}_{chromosome name}`.
///
/// NOTE: We could store the list of genomes in a vector and an index to the vector instead of
/// explicitly storing the genome for each chromosome. However, for convenience, and due to the
/// small number of chromosomes (O(100)), we do not do that.
pub fn genome_of_chrom(reference_path: &Path) -> Result<Vec<GenomeName>> {
    let genomes = ReferenceInfo::from_reference_path(reference_path)?.genomes;
    assert!(!genomes.is_empty());
    let chrom_names: Vec<String> =
        crate::utils::load_txt(&reference_path.join("star/chrName.txt"))?;

    // If we have just one genome in the reference, it's trivial
    if genomes.len() == 1 {
        let genome = genomes[0].clone();
        return Ok(vec![genome; chrom_names.len()]);
    }

    // We have more than 1 genome in the reference. Use the prefix to identify the genomes
    let mut result = Vec::new();
    for name in chrom_names {
        // Find all the genomes that are prefixes of the chromosome name
        let mut genome_of_chrom = Vec::new();
        for genome in &genomes {
            let prefix = format!("{genome}_");
            if name.starts_with(&prefix) {
                genome_of_chrom.push(genome);
            }
        }
        // There should only be one genome.
        assert!(
            genome_of_chrom.len() == 1,
            "Could not find the unique genome of the chromosome {} among the genomes {:?}. Found {} possible genomes: {:?}.",
            name,
            genomes,
            genome_of_chrom.len(),
            genome_of_chrom
        );
        result.push(genome_of_chrom[0].clone());
    }
    Ok(result)
}
