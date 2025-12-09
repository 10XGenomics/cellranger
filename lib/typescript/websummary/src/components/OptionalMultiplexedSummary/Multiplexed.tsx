import React, { useMemo } from "react";
import { Tab } from "react-bootstrap";
import Tabs from "react-bootstrap/Tabs";

import {
  TitledChangeablePlot,
  TitledPlot,
  TitledSelectableBiplotsWithThresholding,
} from "../../shared/websummary/components/CellRanger/hoc";
import { TabContents } from "./AlertElement";
import config from "./config/multiplex-library";
import { MetricTableConfig } from "./config/types";
import { TabTitle } from "./consts";
import ExperimentalDesign from "./ExperimentalDesign";
import Metrics from "./Metrics";
import Sample from "./Sample";
import SampleRollup from "./SampleRollup";
import SuccessCallout from "./SuccessCallout";
import { LibrarySummaryHeader } from "./SummaryHeader";
import { optionalTab, TabContentsProps } from "./TabbedSection";
import {
  ExperimentalDesign as ExperimentalDesignType,
  LibraryType,
  MultiWebSummary,
  MultiWebSummarySampleData,
  SampleWebSummary,
} from "./types";
import { getMetricsForLibraryType, hasAlerts } from "./utilities";
import {
  getAlertSeverity,
  getAlertSeverityFromAlerts,
  getDataTest,
} from "./utils";

const GEX = ({ content, metrics, experimentalDesign }: TabContentsProps) => {
  const libraryMetrics = useMemo(
    () => getMetricsForLibraryType("Gene Expression", metrics),
    [metrics],
  );
  return (
    <>
      <Metrics
        config={config.GEX.hero}
        experimentalDesign={experimentalDesign}
        hero={true}
        metrics={libraryMetrics}
      />

      <div className="half-plot">
        {content?.barcode_rank_plot && (
          <TitledChangeablePlot
            dataTest={getDataTest("barcode_rank_plot")}
            {...content.barcode_rank_plot}
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

      {experimentalDesign.multiplexing_method == "RTL" && (
        <Metrics
          config={config.GEX.perRtlProbeBarcode}
          experimentalDesign={experimentalDesign}
          metrics={libraryMetrics}
        />
      )}

      {experimentalDesign.multiplexing_method == "OH" && (
        <Metrics
          config={config.GEX.perOcmBarcode}
          experimentalDesign={experimentalDesign}
          metrics={libraryMetrics}
        />
      )}

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

      <Metrics
        config={config.GEX.gdna}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />
    </>
  );
};

interface VdjTableConfigs {
  cellCalling: MetricTableConfig;
  enrichment: MetricTableConfig;
  libraryQuality: MetricTableConfig;
  sequencing: MetricTableConfig;
  perOcmBarcode: MetricTableConfig;
  perHashtag: MetricTableConfig;
}

const VDJ_T = (tabProps: TabContentsProps) => {
  return VDJ(tabProps, "VDJ T", config.VDJ_T);
};

const VDJ_T_GD = (tabProps: TabContentsProps) => {
  return VDJ(tabProps, "VDJ T GD", config.VDJ_T_GD);
};

const VDJ_B = (tabProps: TabContentsProps) => {
  return VDJ(tabProps, "VDJ B", config.VDJ_B);
};

const VDJ = (
  tabProps: TabContentsProps,
  libraryType: LibraryType,
  vdjConfig: VdjTableConfigs,
) => {
  const { metrics, content, experimentalDesign } = tabProps;
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

      <Metrics
        config={vdjConfig.enrichment}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />

      <Metrics
        config={vdjConfig.libraryQuality}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />

      <Metrics
        config={vdjConfig.sequencing}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />

      {experimentalDesign.multiplexing_method == "OH" && (
        <Metrics
          config={vdjConfig.perOcmBarcode}
          experimentalDesign={experimentalDesign}
          metrics={libraryMetrics}
        />
      )}

      {experimentalDesign.multiplexing_method == "Hashtag" && (
        <Metrics
          config={vdjConfig.perHashtag}
          experimentalDesign={experimentalDesign}
          metrics={libraryMetrics}
        />
      )}

      <div className="half-plot">
        {content.barcode_rank_plot && (
          <TitledChangeablePlot
            dataTest={getDataTest("barcode_rank_plot")}
            {...content.barcode_rank_plot}
          />
        )}
      </div>
    </>
  );
};

const AB = ({ content, metrics, experimentalDesign }: TabContentsProps) => {
  const libraryMetrics = useMemo(
    () => getMetricsForLibraryType("Antibody Capture", metrics),
    [metrics],
  );
  return (
    <>
      {content?.feature_histogram && (
        <TitledPlot
          dataTest={getDataTest("feature_histogram")}
          {...content.feature_histogram}
        />
      )}

      <Metrics
        config={config.AB.hero}
        experimentalDesign={experimentalDesign}
        hero={true}
        metrics={libraryMetrics}
      />

      <Metrics
        config={config.AB.libraryQuality}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />

      <Metrics
        config={config.AB.sequencing}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />

      {experimentalDesign.multiplexing_method == "RTL" && (
        <Metrics
          config={config.AB.rtl_probe_barcode_metrics}
          experimentalDesign={experimentalDesign}
          metrics={libraryMetrics}
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

const AG = ({ content, metrics, experimentalDesign }: TabContentsProps) => {
  const libraryMetrics = useMemo(
    () => getMetricsForLibraryType("Antigen Capture", metrics),
    [metrics],
  );
  return (
    <>
      <Metrics
        config={config.AG.antigen_reads_table}
        experimentalDesign={experimentalDesign}
        hero={true}
        metrics={libraryMetrics}
      />

      {content?.feature_histogram && (
        <TitledPlot
          dataTest={getDataTest("feature_histogram")}
          {...content.feature_histogram}
        />
      )}

      <Metrics
        config={config.AG.libraryQuality}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />

      <Metrics
        config={config.AG.sequencing}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />

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

const CMO = ({ content, metrics, experimentalDesign }: TabContentsProps) => {
  const libraryMetrics = useMemo(
    () => getMetricsForLibraryType("Multiplexing Capture", metrics),
    [metrics],
  );
  return (
    <>
      <Metrics
        config={config.CMO.hero}
        experimentalDesign={experimentalDesign}
        hero={true}
        metrics={libraryMetrics}
      />

      {content?.jibes_histogram && (
        <TitledPlot
          dataTest={getDataTest("jibes_histogram")}
          {...content.jibes_histogram}
        />
      )}

      <Metrics
        config={config.CMO.perCmoBarcode}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />

      <Metrics
        config={config.CMO.libraryQuality}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />

      <Metrics
        config={config.CMO.sequencing}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />

      {content?.jibes_biplot && (
        <TitledSelectableBiplotsWithThresholding
          dataTest={getDataTest("jibes_biplot")}
          {...content.jibes_biplot}
        />
      )}

      <div className="half-plot">
        {content.cmo_umi_projection_plot && (
          <TitledPlot
            dataTest={getDataTest("cmo_umi_projection_plot")}
            {...content.cmo_umi_projection_plot}
          />
        )}
        {content.cmo_tags_projection_plot && (
          <TitledPlot
            dataTest={getDataTest("cmo_tags_projection_plot")}
            {...content.cmo_tags_projection_plot}
          />
        )}
      </div>
      <div className="half-plot">
        {content?.barcode_rank_plot && (
          <TitledChangeablePlot
            dataTest={getDataTest("barcode_rank_plot")}
            {...content.barcode_rank_plot}
          />
        )}
      </div>
    </>
  );
};

const Hashtag = ({
  content,
  metrics,
  experimentalDesign,
}: TabContentsProps) => {
  const libraryMetrics = useMemo(
    () => getMetricsForLibraryType("Hashtag", metrics),
    [metrics],
  );
  return (
    <>
      <Metrics
        config={config.Hashtag.hero}
        experimentalDesign={experimentalDesign}
        hero={true}
        metrics={libraryMetrics}
      />

      {content?.jibes_histogram && (
        <TitledPlot
          dataTest={getDataTest("jibes_histogram")}
          {...content.jibes_histogram}
        />
      )}

      <Metrics
        config={config.Hashtag.perHashtag}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />

      {content?.jibes_biplot && (
        <TitledSelectableBiplotsWithThresholding
          dataTest={getDataTest("jibes_biplot")}
          {...content.jibes_biplot}
        />
      )}

      <div className="half-plot">
        {content.hashtag_umi_projection_plot && (
          <TitledPlot
            dataTest={getDataTest("hashtag_umi_projection_plot")}
            {...content.hashtag_umi_projection_plot}
          />
        )}
        {content.hashtag_tags_projection_plot && (
          <TitledPlot
            dataTest={getDataTest("hashtag_tags_projection_plot")}
            {...content.hashtag_tags_projection_plot}
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
        config={config.CRISPR.sequencing}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />

      <Metrics
        config={config.CRISPR.rtl_probe_barcode_metrics}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />
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
        config={config.Custom.sequencing}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />

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

const Multiplexed = (props: MultiWebSummary) => {
  const {
    library,
    per_sample: perSample,
    experimental_design: experimentalDesign,
  } = props;

  const anySampleHasAlert = perSample
    .map((sampleData) => hasAlerts(sampleData.data))
    .some((hasAlert) => hasAlert);

  const showSuccessCallout = !hasAlerts(library.data) && !anySampleHasAlert;

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
        sampleId={perSample.length === 1 ? perSample[0].data.id : ""}
      />

      {showSuccessCallout && <SuccessCallout />}

      <Tabs mountOnEnter={true}>
        {optionalTab(
          `${TabTitle.GEX} Library`,
          library.data.gex_tab,
          <GEX
            content={library.data.gex_tab?.content}
            experimentalDesign={experimentalDesign}
            metrics={library.metrics}
          />,
        )}
        {optionalTab(
          `${TabTitle.VDJ_T} Library`,
          library.data.vdj_t_tab,
          <VDJ_T
            content={library.data.vdj_t_tab?.content}
            experimentalDesign={experimentalDesign}
            metrics={library.metrics}
          />,
        )}
        {optionalTab(
          `${TabTitle.VDJ_B} Library`,
          library.data.vdj_b_tab,
          <VDJ_B
            content={library.data.vdj_b_tab?.content}
            experimentalDesign={experimentalDesign}
            metrics={library.metrics}
          />,
        )}
        {optionalTab(
          `${TabTitle.VDJ_T_GD} Library`,
          library.data.vdj_t_gd_tab,
          <VDJ_T_GD
            content={library.data.vdj_t_gd_tab?.content}
            experimentalDesign={experimentalDesign}
            metrics={library.metrics}
          />,
        )}
        {optionalTab(
          `${TabTitle.AB} Library`,
          library.data.antibody_tab,
          <AB
            content={library.data.antibody_tab?.content}
            experimentalDesign={experimentalDesign}
            metrics={library.metrics}
          />,
        )}
        {optionalTab(
          `${TabTitle.AG} Library`,
          library.data.antigen_tab,
          <AG
            content={library.data.antigen_tab?.content}
            experimentalDesign={experimentalDesign}
            metrics={library.metrics}
          />,
        )}
        {optionalTab(
          `${TabTitle.CRISPR} Library`,
          library.data.crispr_tab,
          <CRISPR
            content={library.data.crispr_tab?.content}
            experimentalDesign={experimentalDesign}
            metrics={library.metrics}
          />,
        )}
        {optionalTab(
          `${TabTitle.CMO} Library`,
          library.data.cmo_tab,
          <CMO
            content={library.data.cmo_tab?.content}
            experimentalDesign={experimentalDesign}
            metrics={library.metrics}
          />,
        )}
        {optionalTab(
          TabTitle.Hashtag,
          library.data.hashtag_tab,
          <Hashtag
            content={library.data.hashtag_tab?.content}
            experimentalDesign={experimentalDesign}
            metrics={library.metrics}
          />,
        )}
        {optionalTab(
          `${TabTitle.Custom} Library`,
          library.data.custom_feature_tab,
          <Custom
            content={library.data.custom_feature_tab?.content}
            experimentalDesign={experimentalDesign}
            metrics={library.metrics}
          />,
        )}
        {perSampleSection(perSample, experimentalDesign)}
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

// Dispatch helper for what we put in the per-sample section.
const perSampleSection = (
  perSample: MultiWebSummarySampleData[],
  experimentalDesign: ExperimentalDesignType,
) => {
  if (perSample.length === 0) {
    // Technically this should be unreachable, we should always have at least
    // one sample's worth of data.
    return null;
  }

  const mostSevereSampleAlertLevel = getAlertSeverity(perSample, (sample) => {
    return getAlertSeverityFromAlerts(getAllSampleAlerts(sample.data));
  });

  if (perSample.length === 1) {
    // Standard case, rendering websummary for a single sample.
    const sample = perSample[0];
    return (
      <Tab
        eventKey={"per_sample"}
        key={"per_sample"}
        title={
          <TabContents
            alertType={mostSevereSampleAlertLevel}
            title={"Sample"}
          />
        }
      >
        <Sample
          data={sample}
          experimentalDesign={experimentalDesign}
          key={sample.data.id}
        />
      </Tab>
    );
  }
  // Multiple sample case - render a rolled-up per-sample summary.
  return (
    <Tab
      eventKey={"sample_rollup"}
      key={"sample_rollup"}
      title={
        <TabContents
          alertType={mostSevereSampleAlertLevel}
          title={"Samples Summary"}
        />
      }
    >
      <SampleRollup
        data={perSample}
        experimentalDesign={experimentalDesign}
        key={"sample_rollup"}
      />
    </Tab>
  );
};

/**
 * Collect alerts from all sample tabs.
 */
const getAllSampleAlerts = (ws: SampleWebSummary) =>
  [
    ws.gex_tab?.alerts || [],
    ws.vdj_t_tab?.alerts || [],
    ws.vdj_t_gd_tab?.alerts || [],
    ws.vdj_b_tab?.alerts || [],
    ws.antibody_tab?.alerts || [],
    ws.antigen_tab?.alerts || [],
    ws.crispr_tab?.alerts || [],
    ws.custom_feature_tab?.alerts || [],
    ws.cell_annotation_tab?.alerts || [],
  ].flat();

export default Multiplexed;
