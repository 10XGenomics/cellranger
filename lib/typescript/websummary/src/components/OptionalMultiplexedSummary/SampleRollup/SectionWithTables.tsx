import React from "react";

import { MetricsTableInner } from "../MetricsTable";
import { cleanTableStructure } from "../MetricsTable/utils";
import type { JsonMetricSummary, LibraryType } from "../types";
import { getRowMetrics } from "../utilities";
import { useMetricsPerSamplePerLibraryType } from "./hooks";
import type { MultiSampleTabContentsProps } from "./types";

const sampleIdMetric = (
  id: string,
  libraryType: LibraryType,
): JsonMetricSummary => ({
  alerts: [],
  category: "Cells",
  config: {
    alerts: [],
    header: "Sample ID",
    help: "The ID given to this sample.",
    type: "String",
  },
  key: "sample_id",
  library_type: libraryType,
  value: id,
});

/**
 * Render a section for the provided library type that contains the specified tables.
 */
export const SectionWithTables = ({
  configs,
  experimentalDesign,
  libraryType,
  sampleData,
}: MultiSampleTabContentsProps) => {
  const metricsPerSampleThisLib = useMetricsPerSamplePerLibraryType(
    sampleData.map((d) => d.metrics),
    libraryType,
  );

  // For each config, we will get a table, where each table has a header,
  // rows, etc.
  return configs.map((tableConfig) => {
    let rows = metricsPerSampleThisLib.map((metrics) =>
      getRowMetrics(tableConfig.entries, metrics),
    );

    // Normalize the table structure.
    rows = cleanTableStructure(tableConfig.entries, rows);

    // Prepare each row by first putting the sample id at the front, and
    // then determining whether the row matches the filtering criteria.
    rows.forEach((row, i) => {
      // Push the sample ID into the zeroth column.
      const sampleId = sampleData[i].data.id;
      row.unshift(sampleIdMetric(sampleId, libraryType));
    });

    const initialSortOrder = ["sample_id", ...(tableConfig.orderBy ?? [])];

    return (
      <div data-test="section-with-tables" key={tableConfig.title}>
        <MetricsTableInner
          config={tableConfig}
          experimentalDesign={experimentalDesign}
          hero={false}
          initialSortOrder={initialSortOrder}
          showEmptyMessaging={true}
          tableRows={rows}
        />
      </div>
    );
  });
};
