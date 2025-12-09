import React from "react";
import Tab from "react-bootstrap/Tab";
import Tabs from "react-bootstrap/Tabs";

import { MetricTableConfig } from "../../config/types";
import { MultiWebSummarySampleData } from "../../types";
import { useMetricsPerSamplePerLibraryType } from "../hooks";
import type { MultiSampleTabContentsProps } from "../types";
import { BeeswarmMetrics } from "./BeeswarmMetrics";
import { CollatedPlot } from "./CollatedPlot";
import { PlateAlerts } from "./PlateAlerts";
import { PlateMetrics } from "./PlateMetrics";
import { isPlateBased } from "./plateUtils";
import { groupMetricsWithBarcodes } from "./utils";

interface SectionWithVisualizations extends MultiSampleTabContentsProps {
  // TODO
  // We can be more specific on the types for metric keys.
  beeswarmMetricKeys: string[];
  configs: MetricTableConfig[];
  sampleData: MultiWebSummarySampleData[];
}

function useProps({
  experimentalDesign,
  libraryType,
  sampleData,
}: Pick<
  SectionWithVisualizations,
  "experimentalDesign" | "libraryType" | "sampleData"
>) {
  const sampleIds = sampleData.map((d) => d.data.id);
  const barcodesPerSample = sampleData.map(
    (d) => d.data.multiplexing_barcode_ids,
  );

  const metricsPerSampleThisLib = useMetricsPerSamplePerLibraryType(
    sampleData.map((d) => d.metrics),
    libraryType,
  );
  const metricsBarcodes = groupMetricsWithBarcodes(
    sampleIds,
    barcodesPerSample,
    metricsPerSampleThisLib,
  );

  const showPlateUI =
    barcodesPerSample.length > 0 &&
    isPlateBased(experimentalDesign, barcodesPerSample[0]);

  return { metricsBarcodes, showPlateUI };
}

type ComputedProps = ReturnType<typeof useProps>;

function MetricsPlots({
  configs,
  experimentalDesign,
  beeswarmMetricKeys,
  metricsBarcodes,
  showPlateUI,
}: Pick<
  SectionWithVisualizations,
  "configs" | "experimentalDesign" | "beeswarmMetricKeys"
> &
  ComputedProps) {
  return (
    <div className="summary_border">
      <h2>Metrics</h2>
      {showPlateUI ? (
        <Tabs className="nav-tabs" mountOnEnter={true}>
          <Tab
            data-test="metrics-beeswarm-plot-tab"
            eventKey="metrics-beeswarm-plot"
            title="Beeswarm plot"
          >
            <BeeswarmMetrics
              beeswarmMetricKeys={beeswarmMetricKeys}
              configs={configs}
              experimentalDesign={experimentalDesign}
              sampleMetricsBarcodes={metricsBarcodes}
            />
          </Tab>
          <Tab
            data-test="metrics-heat-map-tab"
            eventKey="metrics-heat-map"
            title="Heat map"
          >
            <PlateMetrics
              configs={configs}
              experimentalDesign={experimentalDesign}
              metricKeys={beeswarmMetricKeys}
              sampleMetricsBarcodes={metricsBarcodes}
            />
          </Tab>
        </Tabs>
      ) : (
        <BeeswarmMetrics
          beeswarmMetricKeys={beeswarmMetricKeys}
          configs={configs}
          experimentalDesign={experimentalDesign}
          sampleMetricsBarcodes={metricsBarcodes}
        />
      )}
    </div>
  );
}

function Antibody(props: SectionWithVisualizations) {
  const {
    configs,
    experimentalDesign,
    libraryType,
    beeswarmMetricKeys,
    sampleData,
  } = props;
  const { showPlateUI, metricsBarcodes } = useProps(props);

  return (
    <>
      {showPlateUI && <PlateAlerts sampleMetricsBarcodes={metricsBarcodes} />}
      <CollatedPlot
        batchDataPerSample={true}
        libraryType={libraryType}
        plotKey="barcode_rank_plot"
        plotTitle="Barcode Rank Plot"
        sampleData={sampleData}
      />
      <MetricsPlots
        beeswarmMetricKeys={beeswarmMetricKeys}
        configs={configs}
        experimentalDesign={experimentalDesign}
        metricsBarcodes={metricsBarcodes}
        showPlateUI={showPlateUI}
      />
    </>
  );
}

function Antigen(props: SectionWithVisualizations) {
  const { configs, experimentalDesign, beeswarmMetricKeys } = props;
  const { showPlateUI, metricsBarcodes } = useProps(props);

  return (
    <>
      {showPlateUI && <PlateAlerts sampleMetricsBarcodes={metricsBarcodes} />}
      {beeswarmMetricKeys.length && (
        <MetricsPlots
          beeswarmMetricKeys={beeswarmMetricKeys}
          configs={configs}
          experimentalDesign={experimentalDesign}
          metricsBarcodes={metricsBarcodes}
          showPlateUI={showPlateUI}
        />
      )}
    </>
  );
}

function CRISPR(props: SectionWithVisualizations) {
  const {
    beeswarmMetricKeys,
    configs,
    experimentalDesign,
    libraryType,
    sampleData,
  } = props;
  const { showPlateUI, metricsBarcodes } = useProps(props);

  return (
    <>
      {showPlateUI && <PlateAlerts sampleMetricsBarcodes={metricsBarcodes} />}
      <CollatedPlot
        batchDataPerSample={true}
        libraryType={libraryType}
        plotKey="barcode_rank_plot"
        plotTitle="Barcode Rank Plot"
        sampleData={sampleData}
      />
      <MetricsPlots
        beeswarmMetricKeys={beeswarmMetricKeys}
        configs={configs}
        experimentalDesign={experimentalDesign}
        metricsBarcodes={metricsBarcodes}
        showPlateUI={showPlateUI}
      />
    </>
  );
}

function Custom(props: SectionWithVisualizations) {
  // We hide plate-based UI this because user-provided CMO names could inadvertently
  // match the regex we use for extracting plate/well locations from sample barcodes.
  const {
    configs,
    experimentalDesign,
    libraryType,
    beeswarmMetricKeys,
    sampleData,
  } = props;
  const { metricsBarcodes } = useProps(props);

  return (
    <>
      <CollatedPlot
        batchDataPerSample={true}
        libraryType={libraryType}
        plotKey="barcode_rank_plot"
        plotTitle="Barcode Rank Plot"
        sampleData={sampleData}
      />
      <MetricsPlots
        beeswarmMetricKeys={beeswarmMetricKeys}
        configs={configs}
        experimentalDesign={experimentalDesign}
        metricsBarcodes={metricsBarcodes}
        showPlateUI={false}
      />
      <CollatedPlot
        batchDataPerSample={false}
        // The legend would be the same value repeated for each sample, so there's no value in
        // showing it.
        legendConfig="none"
        libraryType={libraryType}
        plotKey="median_genes_per_cell_plot"
        plotTitle="Median Genes per Cell"
        sampleData={sampleData}
      />
    </>
  );
}

function GeneExpression(props: SectionWithVisualizations) {
  const {
    configs,
    experimentalDesign,
    libraryType,
    beeswarmMetricKeys,
    sampleData,
  } = props;
  const { showPlateUI, metricsBarcodes } = useProps(props);

  return (
    <>
      {showPlateUI && <PlateAlerts sampleMetricsBarcodes={metricsBarcodes} />}
      <CollatedPlot
        batchDataPerSample={true}
        libraryType={libraryType}
        plotKey="barcode_rank_plot"
        plotTitle="Barcode Rank Plot"
        sampleData={sampleData}
      />
      <MetricsPlots
        beeswarmMetricKeys={beeswarmMetricKeys}
        configs={configs}
        experimentalDesign={experimentalDesign}
        metricsBarcodes={metricsBarcodes}
        showPlateUI={showPlateUI}
      />
      <CollatedPlot
        batchDataPerSample={false}
        // The legend would be the same value repeated for each sample, so there's no value in
        // showing it.
        legendConfig="none"
        libraryType={libraryType}
        plotKey="median_genes_per_cell_plot"
        plotTitle="Median Genes per Cell"
        sampleData={sampleData}
      />
    </>
  );
}

function VDJ(props: SectionWithVisualizations) {
  const {
    configs,
    experimentalDesign,
    libraryType,
    beeswarmMetricKeys,
    sampleData,
  } = props;
  const { showPlateUI, metricsBarcodes } = useProps(props);

  return (
    <>
      {showPlateUI && <PlateAlerts sampleMetricsBarcodes={metricsBarcodes} />}
      <CollatedPlot
        batchDataPerSample={true}
        libraryType={libraryType}
        plotKey="barcode_rank_plot"
        plotTitle="Barcode Rank Plot"
        sampleData={sampleData}
      />
      <MetricsPlots
        beeswarmMetricKeys={beeswarmMetricKeys}
        configs={configs}
        experimentalDesign={experimentalDesign}
        metricsBarcodes={metricsBarcodes}
        showPlateUI={showPlateUI}
      />
    </>
  );
}

/**
 * Render a section for the provided library type that contains the specified tables.
 */
export const SectionWithVisualizations = (props: SectionWithVisualizations) => {
  const { libraryType } = props;

  return (
    <div
      className="section-with-visualizations"
      data-test="section-with-visualizations"
    >
      {libraryType === "Antibody Capture" && <Antibody {...props} />}
      {libraryType === "Antigen Capture" && <Antigen {...props} />}
      {libraryType === "CRISPR Guide Capture" && <CRISPR {...props} />}
      {libraryType === "Custom Feature" && <Custom {...props} />}
      {libraryType === "Gene Expression" && <GeneExpression {...props} />}
      {libraryType === "VDJ B" && <VDJ {...props} />}
      {libraryType === "VDJ T" && <VDJ {...props} />}
      {libraryType === "VDJ T GD" && <VDJ {...props} />}
    </div>
  );
};
