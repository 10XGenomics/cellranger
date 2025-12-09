#![deny(missing_docs)]
use cr_types::rna_read::RnaRead;
use fastq_set::WhichEnd::{FivePrime, ThreePrime};
use fastq_set::adapter_trimmer::{Adapter, AdapterLoc, ReadAdapterCatalog};
use fastq_set::adapters::{
    ILLUMINA_P5_RC, ILLUMINA_P7_RC, ILLUMINA_R1_RC, ILLUMINA_R2_RC, POLY_A, POLY_T, RT_PRIMER,
    RT_PRIMER_RC, SPACER, SPACER_RC,
};
use fastq_set::read_pair::WhichRead;
use metric::{PercentMetric, TxHashMap};

pub(super) fn get_vdj_adapters() -> TxHashMap<WhichRead, Vec<Adapter>> {
    use AdapterLoc::{Anywhere, NonInternal};
    // The molecules which are sequenced looks like
    // P5-R1-BC-UMI-SPACER-INSERT-R2-SI-P7

    let mut adapters = TxHashMap::default();

    adapters.insert(
        WhichRead::R1,
        vec![
            // For the SCVDJ chemistry, we trim the first 41 bases of read1
            // But it is safe to check for the spacer adapter in case a custom
            // chemistry is specified
            Adapter::new("spacer", FivePrime, NonInternal, SPACER),
            // Adapters expected at the 3' end of Read1
            Adapter::new("R2_rc", ThreePrime, Anywhere, ILLUMINA_R2_RC),
            Adapter::new("P7_rc", ThreePrime, Anywhere, ILLUMINA_P7_RC),
            Adapter::new("polyA", ThreePrime, Anywhere, POLY_A),
            Adapter::new("rt_primer_rc", ThreePrime, Anywhere, RT_PRIMER_RC),
        ],
    );

    adapters.insert(
        WhichRead::R2,
        vec![
            Adapter::new("spacer_rc", ThreePrime, Anywhere, SPACER_RC),
            Adapter::new("R1_rc", ThreePrime, Anywhere, ILLUMINA_R1_RC),
            Adapter::new("P5_rc", ThreePrime, Anywhere, ILLUMINA_P5_RC),
            Adapter::new("polyT", FivePrime, Anywhere, POLY_T),
            Adapter::new("rt_primer", FivePrime, NonInternal, RT_PRIMER),
        ],
    );

    adapters
}

pub(super) struct VdjTrimmer<'a> {
    adapter_catalog: ReadAdapterCatalog<'a>,
    metrics: TxHashMap<WhichRead, TxHashMap<String, PercentMetric>>,
}

impl<'a> VdjTrimmer<'a> {
    pub(super) fn new(adapter_map: &'a TxHashMap<WhichRead, Vec<Adapter>>) -> Self {
        let mut metrics = TxHashMap::default();
        for (&which_read, adapters) in adapter_map {
            let mut adapter_frac = TxHashMap::default();
            for adapter in adapters {
                adapter_frac.insert(adapter.name.clone(), PercentMetric::default());
            }
            metrics.insert(which_read, adapter_frac);
        }

        VdjTrimmer {
            adapter_catalog: ReadAdapterCatalog::from(adapter_map),
            metrics,
        }
    }

    pub(super) fn trim(&mut self, rna_read: &mut RnaRead) {
        let adapter_positions = rna_read.trim_adapters(&mut self.adapter_catalog);
        for metrics in self.metrics.values_mut() {
            for (name, metric) in metrics {
                metric.increment(adapter_positions.contains_key(name));
            }
        }
    }
}
