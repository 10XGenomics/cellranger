import React from "react";

import { GroupedBeeswarmPlot } from "../../../../shared/websummary/components/BeeswarmPlot";
import { MetricTableConfig } from "../../config/types";
import { ExperimentalDesign } from "../../types";
import { AlertsThresholdLegend } from "./AlertsThresholdLegend";
import { SECTION_STYLE } from "./common";
import { useBeeswarmMetricsData } from "./useBeeswarmMetricsData";
import { type SampleMetricsBarcodes } from "./utils";
import { VisualizationWrapper } from "./VisualizationWrapper";

export function BeeswarmMetrics({
  configs,
  experimentalDesign,
  sampleMetricsBarcodes,
  beeswarmMetricKeys,
}: {
  configs: MetricTableConfig[];
  experimentalDesign: ExperimentalDesign;
  // List of metric keys for which to display these plots for. Anything
  // not in this list will be skipped when we process the data.
  beeswarmMetricKeys: string[];
  sampleMetricsBarcodes: SampleMetricsBarcodes;
}) {
  const beeswarmMetrics = useBeeswarmMetricsData(
    sampleMetricsBarcodes,
    configs,
    experimentalDesign,
    beeswarmMetricKeys,
  );

  if (beeswarmMetricKeys.length === 0) {
    return null;
  }

  // Layout:
  // 3 beeswarm plots per row.
  // Spacing of 7% between them.
  // (100% - 7% * 2) = 86%
  // So 86% / 3 is percentage width we have for the 3 plots.
  const plotsPerRow = 3;
  const gapPercent = 7;
  const widthPercent = (100 - gapPercent * (plotsPerRow - 1)) / plotsPerRow;
  return (
    <>
      <AlertsThresholdLegend />
      <section
        style={{
          display: "flex",
          flexFlow: "row wrap",
          gap: `0 ${gapPercent}%`,
          marginTop: 16,
          ...SECTION_STYLE,
        }}
      >
        {beeswarmMetrics.map(
          ({ groups, helpText, key, title, typedThreshold }) => (
            <VisualizationWrapper
              helpText={helpText}
              key={key}
              style={{ width: `${widthPercent}%` }}
              title={title}
              uniqueId={key}
            >
              <GroupedBeeswarmPlot
                data={{
                  groups,
                  typedThreshold: typedThreshold ?? undefined,
                }}
                style={{
                  width: "100%",
                }}
                title=""
              />
            </VisualizationWrapper>
          ),
        )}
      </section>
    </>
  );
}
