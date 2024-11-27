//! Martian stage WRITE_MATRIX_MARKET.

use crate::env;
use crate::types::FeatureReferenceFormat;
use crate::utils::estimate_mem::{barcode_mem_gib, get_total_barcodes_detected};
use anyhow::Result;
use barcode::Barcode;
use cr_types::barcode_index::BarcodeIndex;
use cr_types::reference::feature_reference::{FeatureReference, FeatureType};
use cr_types::types::{BarcodeIndexFormat, FeatureBarcodeCount};
use cr_types::{BarcodeThenFeatureOrder, CountShardFile, MetricsFile};
use martian::{MartianFileType, MartianRover, MartianStage, MartianVoid, Resource, StageDef};
use martian_derive::{make_mro, MartianStruct};
use martian_filetypes::gzip_file::Gzip;
use martian_filetypes::tabular_file::TsvFileNoHeader;
use martian_filetypes::{FileTypeRead, FileTypeWrite};
use serde::{Deserialize, Serialize};
use shardio::ShardReader;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Serialize, Deserialize)]
struct FeatureTsvRow {
    // `ENSG00000243485` for example
    id: String,
    // `MIR1302-2HG` for example
    name: String,
    // `Gene Expression` for example
    feature_type: FeatureType,
}

struct MtxWriter {
    folder: PathBuf,
}

impl MtxWriter {
    fn new(folder: &Path) -> Result<Self> {
        std::fs::create_dir(folder)?;
        Ok(MtxWriter {
            folder: folder.to_path_buf(),
        })
    }
    fn write_files(
        &self,
        chunks: &[CountShardFile],
        feature_ref: &FeatureReference,
        barcode_index: &BarcodeIndex,
        software_version: &str,
    ) -> Result<()> {
        self.write_matrix_mtx(chunks, feature_ref, barcode_index, software_version)?;
        self.write_barcodes_tsv(barcode_index.sorted_barcodes())?;
        self.write_features_tsv(feature_ref)
    }

    fn write_barcodes_tsv(&self, barcodes: &[Barcode]) -> Result<()> {
        let barcodes_path = self.folder.join("barcodes.tsv.gz");
        let mut barcodes_writer = BufWriter::new(flate2::write::GzEncoder::new(
            std::fs::File::create(barcodes_path)?,
            flate2::Compression::fast(),
        ));
        for bc in barcodes {
            writeln!(barcodes_writer, "{bc}")?;
        }
        barcodes_writer.flush()?;
        Ok(())
    }

    fn write_features_tsv(&self, feature_ref: &FeatureReference) -> Result<()> {
        let rows = feature_ref
            .feature_defs
            .iter()
            .map(|fdef| FeatureTsvRow {
                id: fdef.id.clone(),
                name: fdef.name.clone(),
                feature_type: fdef.feature_type,
            })
            .collect();
        Gzip::<TsvFileNoHeader<_>>::new(&self.folder, "features").write(&rows)
    }

    fn write_matrix_mtx(
        &self,
        chunks: &[CountShardFile],
        feature_ref: &FeatureReference,
        barcode_index: &BarcodeIndex,
        software_version: &str,
    ) -> Result<()> {
        let mtx_path = self.folder.join("matrix.mtx.gz");
        let mut writer = BufWriter::new(flate2::write::GzEncoder::new(
            std::fs::File::create(mtx_path)?,
            flate2::Compression::fast(),
        ));
        let reader: ShardReader<FeatureBarcodeCount, BarcodeThenFeatureOrder> =
            ShardReader::open_set(chunks)?;

        // Print the header, matrix dimensions, and number of non-zero entries.
        writeln!(
            writer,
            r#"%%MatrixMarket matrix coordinate integer general
%metadata_json: {{"software_version": "{} {}", "format_version": 2}}
{} {} {}"#,
            env::get_tenx_product_name(),
            software_version,
            feature_ref.feature_defs.len(),
            barcode_index.len(),
            reader.len()
        )?;

        for count in reader.iter()? {
            let count = count?;
            // indices are 1-based
            writeln!(
                writer,
                "{} {} {}",
                1 + count.feature_idx,
                1 + barcode_index.get_index(&count.barcode),
                count.umi_count
            )?;
        }

        writer.flush()?;
        Ok(())
    }
}

/// Output the feature-barcode Matrix Market file `raw_feature_bc_matrix/matrix.mtx.gz`.
pub struct WriteMatrixMarket;

#[derive(Deserialize, Clone, MartianStruct)]
pub struct WriteMatrixMarketStageInputs {
    pub counts: Vec<CountShardFile>,
    pub feature_reference: FeatureReferenceFormat,
    pub barcode_index: BarcodeIndexFormat,
    pub barcode_correction_summary: MetricsFile,
}

#[derive(Serialize, Deserialize, Clone, MartianStruct)]
pub struct WriteMatrixMarketStageOutputs {
    pub feature_bc_matrix: PathBuf,
}

#[make_mro(volatile = strict)]
impl MartianStage for WriteMatrixMarket {
    type StageInputs = WriteMatrixMarketStageInputs;
    type StageOutputs = WriteMatrixMarketStageOutputs;
    type ChunkInputs = MartianVoid;
    type ChunkOutputs = MartianVoid;

    fn split(
        &self,
        args: Self::StageInputs,
        _rover: MartianRover,
    ) -> Result<StageDef<Self::ChunkInputs>> {
        let barcodes_count = get_total_barcodes_detected(&args.barcode_correction_summary.read()?);
        // bytes_per_barcode and offset_gib are empirically determined.
        let mem_gib = barcode_mem_gib(barcodes_count, 220, 2);
        println!("barcode_count={barcodes_count},mem_gib={mem_gib}");
        Ok(StageDef::with_join_resource(Resource::with_mem_gb(mem_gib)))
    }

    fn main(
        &self,
        _args: Self::StageInputs,
        _chunk_args: MartianVoid,
        _rover: MartianRover,
    ) -> Result<Self::ChunkOutputs> {
        unreachable!()
    }

    fn join(
        &self,
        args: Self::StageInputs,
        _chunk_defs: Vec<MartianVoid>,
        _chunk_outs: Vec<MartianVoid>,
        rover: MartianRover,
    ) -> Result<Self::StageOutputs> {
        let barcode_index = args.barcode_index.read()?;
        let feature_bc_matrix: PathBuf = rover.make_path("raw_feature_bc_matrix");
        let mtx_writer = MtxWriter::new(&feature_bc_matrix)?;
        mtx_writer.write_files(
            &args.counts,
            &args.feature_reference.read()?,
            &barcode_index,
            &rover.pipelines_version(),
        )?;
        Ok(Self::StageOutputs { feature_bc_matrix })
    }
}
