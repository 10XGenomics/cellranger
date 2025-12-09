import React, { useMemo } from "react";

import { TitledPlot } from "../../../shared/websummary/components/CellRanger/hoc";
import config from "../config/multiplex-sample";
import Metrics from "../Metrics";
import { TabContentsProps } from "../TabbedSection";
import { getMetricsForLibraryType } from "../utilities";
import { getDataTest } from "../utils";

export const Antigen = ({
  content,
  metrics,
  experimentalDesign,
}: TabContentsProps) => {
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

      {content.antigen_treemap && (
        <TitledPlot
          dataTest={getDataTest("antigen_treemap")}
          {...content.antigen_treemap}
        />
      )}

      {content.clonotype_clustermap && (
        <TitledPlot
          dataTest={getDataTest("clonotype_clustermap")}
          {...content.clonotype_clustermap}
        />
      )}
    </>
  );
};
