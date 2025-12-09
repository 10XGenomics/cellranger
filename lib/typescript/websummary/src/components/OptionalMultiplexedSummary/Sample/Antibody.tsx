import React, { useMemo } from "react";

import ClusteringSelector from "../../../shared/websummary/components/CellRanger/ClusteringSelector";
import { TitledPlot } from "../../../shared/websummary/components/CellRanger/hoc";
import config from "../config/multiplex-sample";
import Metrics from "../Metrics";
import { TabContentsProps } from "../TabbedSection";
import { getMetricsForLibraryType } from "../utilities";
import { getDataTest } from "../utils";

export const Antibody = ({
  content,
  metrics,
  experimentalDesign,
}: TabContentsProps) => {
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

      {content.antibody_treemap && (
        <TitledPlot
          dataTest={getDataTest("antibody_treemap")}
          {...content.antibody_treemap}
        />
      )}

      {content.clustering_and_diffexp_plots && (
        <ClusteringSelector
          dataTest={getDataTest("clustering_and_diffexp_plots")}
          {...content.clustering_and_diffexp_plots}
        />
      )}

      <div className="half-plot">
        {content.barcode_rank_plot && (
          <TitledPlot
            dataTest={getDataTest("barcode_rank_plot")}
            {...content.barcode_rank_plot}
          />
        )}
      </div>
    </>
  );
};
