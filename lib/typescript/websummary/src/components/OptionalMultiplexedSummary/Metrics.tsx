import React, { useMemo } from "react";

import { MetricTableConfig } from "./config/types";
import { MetricsTable } from "./MetricsTable";
import MetricsTableSingleRow from "./MetricsTableSingleRow";
import type { ExperimentalDesign, JsonMetricSummary } from "./types";
import { getRowMetrics } from "./utilities";

interface Props {
  metrics: JsonMetricSummary[];
  config: MetricTableConfig;
  experimentalDesign: ExperimentalDesign;
  hero?: boolean;
}

const Metrics = ({ metrics, config, experimentalDesign, hero }: Props) => {
  const { entries } = config;
  const shouldRender = useMemo(
    () => getRowMetrics(entries, metrics).length > 0,
    [entries, metrics],
  );
  // If we have more than one instance of the first metric config entry,
  // and the metrics have grouping keys, assume this is multi-row data.
  const useTable =
    metrics.filter(
      (metric) => entries[0].key === metric.key && metric.grouping_key,
    ).length > 1;

  if (!shouldRender) {
    return null;
  }

  return (
    <>
      {useTable ? (
        <MetricsTable
          config={config}
          experimentalDesign={experimentalDesign}
          hero={!!hero}
          metrics={metrics}
        />
      ) : (
        <MetricsTableSingleRow
          config={config}
          experimentalDesign={experimentalDesign}
          hero={!!hero}
          metrics={metrics}
        />
      )}
    </>
  );
};

export default Metrics;
