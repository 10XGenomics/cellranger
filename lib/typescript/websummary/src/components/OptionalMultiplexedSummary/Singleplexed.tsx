import { mergeWith } from "lodash";
import React, { useMemo } from "react";
import { Tab } from "react-bootstrap";
import Tabs from "react-bootstrap/Tabs";

import ClusteringSelector from "../../shared/websummary/components/CellRanger/ClusteringSelector";
import {
  TitledChangeablePlot,
  TitledPlot,
  TitledPlotTable,
} from "../../shared/websummary/components/CellRanger/hoc";
import { TabContents } from "./AlertElement";
import { SinglePlexConfig as config } from "./config/singleplex";
import { MetricTableConfig } from "./config/types";
import { TabTitle } from "./consts";
import ExperimentalDesign from "./ExperimentalDesign";
import Metrics from "./Metrics";
import SuccessCallout from "./SuccessCallout";
import { LibrarySummaryHeader } from "./SummaryHeader";
import { optionalTab, TabContentsProps } from "./TabbedSection";
import {
  DiffExpTableCardAndTitle,
  VegaLitePlotWithCardAndTitle,
} from "./TitledComponents";
import { LibraryType, MultiWebSummary } from "./types";
import { getMetricsForLibraryType, hasAlerts } from "./utilities";
import { getDataTest } from "./utils";

const GEX = ({ content, metrics, experimentalDesign }: TabContentsProps) => {
  const libraryMetrics = useMemo(
    () => getMetricsForLibraryType("Gene Expression", metrics),
    [metrics],
  );

  return (
    <>
      <Metrics
        config={
          experimentalDesign.is_barnyard
            ? config.GEX.cellCallingBarnyard
            : config.GEX.cellCalling
        }
        experimentalDesign={experimentalDesign}
        hero={true}
        metrics={libraryMetrics}
      />

      {experimentalDesign.is_barnyard && (
        <Metrics
          config={config.GEX.barnyardPerGenome}
          experimentalDesign={experimentalDesign}
          metrics={libraryMetrics}
        />
      )}
      <div className="half-plot">
        <TitledChangeablePlot
          dataTest={getDataTest("barcode_rank_plot")}
          {...content.barcode_rank_plot}
        />
        {experimentalDesign.is_barnyard && (
          <TitledChangeablePlot
            dataTest={getDataTest("barnyard_biplot")}
            {...content.barnyard_biplot}
          />
        )}
      </div>
      <Metrics
        config={config.GEX.mapping}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />
      <Metrics
        config={config.GEX.libraryQuality}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />
      <Metrics
        config={config.GEX.sequencing}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />
      <div className="half-plot">
        {content?.sequencing_saturation_plot && (
          <TitledPlot
            dataTest={getDataTest("sequencing_saturation_plot")}
            {...content.sequencing_saturation_plot}
          />
        )}
        {content?.median_genes_per_cell_plot && (
          <TitledPlot
            dataTest={getDataTest("median_genes_per_cell_plot")}
            {...content.median_genes_per_cell_plot}
          />
        )}
      </div>
      {content?.clustering_and_diffexp_plots && (
        <ClusteringSelector
          dataTest={getDataTest("clustering_and_diffexp_plots")}
          {...content.clustering_and_diffexp_plots}
        />
      )}
      <Metrics
        config={config.GEX.gdna_metrics}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />
    </>
  );
};

const VDJ_T = (tabProps: TabContentsProps) => {
  return VDJ(tabProps, "VDJ T", config.VDJ_T);
};

const VDJ_T_GD = (tabProps: TabContentsProps) => {
  return VDJ(tabProps, "VDJ T GD", config.VDJ_T_GD);
};

const VDJ_B = (tabProps: TabContentsProps) => {
  return VDJ(tabProps, "VDJ B", config.VDJ_B);
};

interface VdjTableConfigs {
  cellCalling: MetricTableConfig;
  enrichmentSensitivity: MetricTableConfig;
  libraryQuality: MetricTableConfig;
  sequencingQuality: MetricTableConfig;
  clonotype: MetricTableConfig;
}

function VDJ(
  tabProps: TabContentsProps,
  libraryType: LibraryType,
  vdjConfig: VdjTableConfigs,
) {
  const { metrics, content: data, experimentalDesign } = tabProps;
  const libraryMetrics = useMemo(
    () => getMetricsForLibraryType(libraryType, metrics),
    [libraryType, metrics],
  );
  return (
    <>
      <Metrics
        config={vdjConfig.cellCalling}
        experimentalDesign={experimentalDesign}
        hero={true}
        metrics={libraryMetrics}
      />

      {data.barcode_rank_plot && (
        <div className="half-plot">
          {data?.barcode_rank_plot && (
            <TitledChangeablePlot
              dataTest={getDataTest("barcode_rank_plot")}
              {...data.barcode_rank_plot}
            />
          )}
        </div>
      )}

      <Metrics
        config={vdjConfig.enrichmentSensitivity}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />
      <Metrics
        config={vdjConfig.libraryQuality}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />
      <Metrics
        config={vdjConfig.sequencingQuality}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />

      <Metrics
        config={vdjConfig.clonotype}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />

      {data.clonotype_info && (
        <TitledPlotTable
          dataTest={getDataTest("clonotype_info")}
          {...data.clonotype_info}
        />
      )}
    </>
  );
}

const AB = ({ content, metrics, experimentalDesign }: TabContentsProps) => {
  const libraryMetrics = useMemo(
    () => getMetricsForLibraryType("Antibody Capture", metrics),
    [metrics],
  );
  return (
    <>
      <Metrics
        config={config.AB.hero}
        experimentalDesign={experimentalDesign}
        hero={true}
        metrics={libraryMetrics}
      />

      {content.feature_histogram && (
        <TitledPlot
          dataTest={getDataTest("feature_histogram")}
          {...content.feature_histogram}
        />
      )}
      <Metrics
        config={config.AB.libraryQuality}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />
      <Metrics
        config={config.AB.sequencingQuality}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />

      {content.antibody_treemap && (
        <TitledPlot
          dataTest={getDataTest("antibody_treemap")}
          {...content.antibody_treemap}
        />
      )}

      <div className="half-plot">
        {content.barcode_rank_plot && (
          <TitledPlot
            dataTest={getDataTest("barcode_rank_plot")}
            {...content.barcode_rank_plot}
          />
        )}
        {content?.sequencing_saturation_plot && (
          <TitledPlot
            dataTest={getDataTest("sequencing_saturation_plot")}
            {...content.sequencing_saturation_plot}
          />
        )}
      </div>

      {content.clustering_and_diffexp_plots && (
        <ClusteringSelector
          dataTest={getDataTest("clustering_and_diffexp_plots")}
          {...content.clustering_and_diffexp_plots}
        />
      )}

      <div className="half-plot">
        {content.projection_plot && (
          <TitledPlot
            dataTest={getDataTest("projection_plot")}
            {...content.projection_plot}
          />
        )}
      </div>
    </>
  );
};

const AG = ({ content, metrics, experimentalDesign }: TabContentsProps) => {
  const libraryMetrics = useMemo(
    () => getMetricsForLibraryType("Antigen Capture", metrics),
    [metrics],
  );
  return (
    <>
      <Metrics
        config={config.AG.hero}
        experimentalDesign={experimentalDesign}
        hero={true}
        metrics={libraryMetrics}
      />

      <Metrics
        config={config.AG.expression}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />

      {content.clonotype_clustermap && (
        <TitledPlot
          dataTest={getDataTest("clonotype_clustermap")}
          {...content.clonotype_clustermap}
        />
      )}

      <Metrics
        config={config.AG.libraryQuality}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />

      <Metrics
        config={config.AG.sequencingQuality}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />

      {content.antigen_treemap && (
        <TitledPlot
          dataTest={getDataTest("antigen_treemap")}
          {...content.antigen_treemap}
        />
      )}

      {content?.feature_histogram && (
        <TitledPlot
          dataTest={getDataTest("feature_histogram")}
          {...content.feature_histogram}
        />
      )}

      <div className="half-plot">
        {content?.barcode_rank_plot && (
          <TitledChangeablePlot
            dataTest={getDataTest("barcode_rank_plot")}
            {...content.barcode_rank_plot}
          />
        )}
        {content?.sequencing_saturation_plot && (
          <TitledPlot
            dataTest={getDataTest("sequencing_saturation_plot")}
            {...content.sequencing_saturation_plot}
          />
        )}
      </div>
    </>
  );
};

const CRISPR = ({ content, metrics, experimentalDesign }: TabContentsProps) => {
  const libraryMetrics = useMemo(
    () => getMetricsForLibraryType("CRISPR Guide Capture", metrics),
    [metrics],
  );
  return (
    <>
      <Metrics
        config={config.CRISPR.hero}
        experimentalDesign={experimentalDesign}
        hero={true}
        metrics={libraryMetrics}
      />
      <Metrics
        config={config.CRISPR.libraryQuality}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />
      <Metrics
        config={config.CRISPR.sequencingQuality}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />

      <div className="half-plot">
        {content.barcode_rank_plot && (
          <TitledPlot
            dataTest={getDataTest("barcode_rank_plot")}
            {...content.barcode_rank_plot}
          />
        )}
        {content.sequencing_saturation_plot && (
          <TitledPlot
            dataTest={getDataTest("sequencing_saturation_plot")}
            {...content.sequencing_saturation_plot}
          />
        )}
      </div>

      <div className="half-plot">
        {content.projection_plot && (
          <TitledPlot
            dataTest={getDataTest("projection_plot")}
            {...content.projection_plot}
          />
        )}
      </div>
    </>
  );
};

const Custom = ({ content, metrics, experimentalDesign }: TabContentsProps) => {
  const libraryMetrics = useMemo(
    () => getMetricsForLibraryType("Custom Feature", metrics),
    [metrics],
  );
  return (
    <>
      <Metrics
        config={config.Custom.hero}
        experimentalDesign={experimentalDesign}
        hero={true}
        metrics={libraryMetrics}
      />

      <Metrics
        config={config.Custom.libraryQuality}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />

      <Metrics
        config={config.Custom.sequencingQuality}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />

      <div className="half-plot">
        {content.barcode_rank_plot && (
          <TitledPlot
            dataTest={getDataTest("barcode_rank_plot")}
            {...content.barcode_rank_plot}
          />
        )}
        {content.sequencing_saturation_plot && (
          <TitledPlot
            dataTest={getDataTest("sequencing_saturation_plot")}
            {...content.sequencing_saturation_plot}
          />
        )}
      </div>

      <div className="half-plot">
        {content.projection_plot && (
          <TitledPlot
            dataTest={getDataTest("projection_plot")}
            {...content.projection_plot}
          />
        )}
      </div>
    </>
  );
};

const Annotate = ({ content }: TabContentsProps) => {
  return (
    <>
      {content.cell_annotation_cell_types_chart && (
        <VegaLitePlotWithCardAndTitle
          dataTest={getDataTest("cell_annotation_cell_types_chart")}
          {...content.cell_annotation_cell_types_chart}
        />
      )}
      <div
        style={{
          display: "flex",
          flexDirection: "row",
          flexWrap: "wrap",
          gap: "2rem",
          justifyContent: "space-between",
        }}
      >
        {content.cell_annotation_violin_plot_chart && (
          <VegaLitePlotWithCardAndTitle
            dataTest={getDataTest("cell_annotation_violin_plot_chart")}
            {...content.cell_annotation_violin_plot_chart}
          />
        )}
        {content.cell_annotation_umap_plot_chart && (
          <VegaLitePlotWithCardAndTitle
            dataTest={getDataTest("cell_annotation_umap_plot_chart")}
            {...content.cell_annotation_umap_plot_chart}
          />
        )}
      </div>
      {content.cell_annotation_diffexp_table && (
        <DiffExpTableCardAndTitle
          dataTest={getDataTest("cell_annotation_diffexp_table")}
          {...content.cell_annotation_diffexp_table}
        />
      )}
    </>
  );
};

// lodash merge but also ignore nulls
function mergeIgnoreNull<T>(a: T, b: T): T {
  return mergeWith({}, a, b, (x, y) => (y === null ? x : undefined));
}

// Top-level component for a non-multiplexed analysis run through the multi pipeline.
const Singleplexed = (props: MultiWebSummary) => {
  const {
    library,
    per_sample: perSample,
    experimental_design: experimentalDesign,
  } = props;

  const showSuccessCallout =
    !hasAlerts(library.data) && !hasAlerts(perSample[0].data);

  const combinedMetrics = [...perSample[0].metrics, ...library.metrics];

  const gexTabData = mergeIgnoreNull(
    perSample[0].data.gex_tab,
    library.data.gex_tab,
  );

  const vdjBTabData = mergeIgnoreNull(
    perSample[0].data.vdj_b_tab,
    library.data.vdj_b_tab,
  );

  const vdjTTabData = mergeIgnoreNull(
    perSample[0].data.vdj_t_tab,
    library.data.vdj_t_tab,
  );
  const vdjTGdTabData = mergeIgnoreNull(
    perSample[0].data.vdj_t_gd_tab,
    library.data.vdj_t_gd_tab,
  );
  const antibodyTabData = mergeIgnoreNull(
    perSample[0].data.antibody_tab,
    library.data.antibody_tab,
  );
  const antigenTabData = mergeIgnoreNull(
    perSample[0].data.antigen_tab,
    library.data.antigen_tab,
  );
  const crisprTabData = mergeIgnoreNull(
    perSample[0].data.crispr_tab,
    library.data.crispr_tab,
  );
  const customFeatureTabData = mergeIgnoreNull(
    perSample[0].data.custom_feature_tab,
    library.data.custom_feature_tab,
  );

  return (
    <article
      className="summary"
      data-test={getDataTest("library")}
      key="Library"
    >
      <LibrarySummaryHeader
        description={library.data.description}
        id={library.data.id}
        pipelineVersion={props.pipeline_version}
        sampleId={perSample[0].data.id}
      />

      {showSuccessCallout && <SuccessCallout />}

      <Tabs mountOnEnter={true}>
        {optionalTab(
          TabTitle.GEX,
          gexTabData,
          <GEX
            content={gexTabData?.content}
            experimentalDesign={experimentalDesign}
            metrics={combinedMetrics}
          />,
        )}
        {optionalTab(
          TabTitle.VDJ_B,
          vdjBTabData,
          <VDJ_B
            content={vdjBTabData?.content}
            experimentalDesign={experimentalDesign}
            metrics={combinedMetrics}
          />,
        )}
        {optionalTab(
          TabTitle.VDJ_T,
          vdjTTabData,
          <VDJ_T
            content={vdjTTabData?.content}
            experimentalDesign={experimentalDesign}
            metrics={combinedMetrics}
          />,
        )}
        {optionalTab(
          TabTitle.VDJ_T_GD,
          vdjTGdTabData,
          <VDJ_T_GD
            content={vdjTGdTabData?.content}
            experimentalDesign={experimentalDesign}
            metrics={combinedMetrics}
          />,
        )}
        {optionalTab(
          TabTitle.AB,
          antibodyTabData,
          <AB
            content={antibodyTabData?.content}
            experimentalDesign={experimentalDesign}
            metrics={combinedMetrics}
          />,
        )}
        {optionalTab(
          TabTitle.AG,
          antigenTabData,
          <AG
            content={antigenTabData?.content}
            experimentalDesign={experimentalDesign}
            metrics={combinedMetrics}
          />,
        )}
        {optionalTab(
          TabTitle.CRISPR,
          crisprTabData,
          <CRISPR
            content={crisprTabData?.content}
            experimentalDesign={experimentalDesign}
            metrics={combinedMetrics}
          />,
        )}
        {optionalTab(
          TabTitle.Custom,
          customFeatureTabData,
          <Custom
            content={customFeatureTabData?.content}
            experimentalDesign={experimentalDesign}
            metrics={combinedMetrics}
          />,
        )}
        {optionalTab(
          TabTitle.Annotate,
          perSample[0].data.cell_annotation_tab,
          <Annotate
            content={perSample[0].data.cell_annotation_tab?.content}
            experimentalDesign={experimentalDesign}
            metrics={combinedMetrics}
          />,
        )}
        <Tab
          eventKey={"experimental_design"}
          key={"experimental_design"}
          title={<TabContents title={"Experimental Design"} />}
        >
          <ExperimentalDesign {...experimentalDesign} />
        </Tab>
      </Tabs>
    </article>
  );
};

export default Singleplexed;
