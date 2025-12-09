import React from "react";

import { MetricTableConfig } from "../../config/types";
import { ExperimentalDesign } from "../../types";
import { AlertsThresholdLegend } from "./AlertsThresholdLegend";
import { ColorScale } from "./ColorScale";
import { SECTION_STYLE } from "./common";
import { PlateGroup } from "./PlateGroup";
import { usePlateMetricsData } from "./usePlateMetricsData";
import { type SampleMetricsBarcodes } from "./utils";
import { VisualizationWrapper } from "./VisualizationWrapper";

export function PlateMetrics({
  configs,
  experimentalDesign,
  sampleMetricsBarcodes,
  metricKeys,
}: {
  configs: MetricTableConfig[];
  experimentalDesign: ExperimentalDesign;
  // List of metric keys for which to display these plots for. Anything
  // not in this list will be skipped when we process the data.
  metricKeys: string[];
  sampleMetricsBarcodes: SampleMetricsBarcodes;
}) {
  const plateMetrics = usePlateMetricsData(
    sampleMetricsBarcodes,
    configs,
    experimentalDesign,
    metricKeys,
  );

  return (
    <>
      <AlertsThresholdLegend />
      <section style={{ marginTop: 16, ...SECTION_STYLE }}>
        {plateMetrics.map(
          ({
            bounds,
            colorScale,
            helpText,
            key,
            plateDataMap,
            title,
            typedThreshold,
          }) => (
            <VisualizationWrapper
              helpText={helpText}
              key={key}
              style={{ marginBottom: "56px" }}
              title={title}
              uniqueId={key}
            >
              <ColorScale
                bounds={bounds}
                colorScale={colorScale}
                typedTheshold={typedThreshold ?? undefined}
              />
              <PlateGroup plateDataMap={plateDataMap} />
            </VisualizationWrapper>
          ),
        )}
      </section>
    </>
  );
}
