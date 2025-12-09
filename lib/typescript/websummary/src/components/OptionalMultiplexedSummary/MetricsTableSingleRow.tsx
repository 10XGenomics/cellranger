import React, { useMemo } from "react";

import { MetricTableConfig } from "./config/types";
import { MetricsTableInner } from "./MetricsTable";
import type { ExperimentalDesign, JsonMetricSummary } from "./types";
import { getRowMetrics } from "./utilities";

interface Props {
  metrics: JsonMetricSummary[];
  config: MetricTableConfig;
  experimentalDesign: ExperimentalDesign;
  hero: boolean;
}

/**
 * Render a table of metrics that only has one row of data.
 */
const MetricsTableSingleRow = ({
  metrics,
  config,
  experimentalDesign,
  hero,
}: Props) => {
  const rowMetrics = useMemo(
    () => getRowMetrics(config.entries, metrics),
    [config.entries, metrics],
  );

  return (
    <MetricsTableInner
      config={config}
      experimentalDesign={experimentalDesign}
      hero={hero}
      initialSortOrder={[]}
      showEmptyMessaging={false}
      tableRows={[rowMetrics]}
    />
  );
};

export default MetricsTableSingleRow;
