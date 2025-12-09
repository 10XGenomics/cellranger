import React, { useMemo } from "react";

import ClusteringSelector from "../../../shared/websummary/components/CellRanger/ClusteringSelector";
import { TitledPlot } from "../../../shared/websummary/components/CellRanger/hoc";
import config from "../config/multiplex-sample";
import Metrics from "../Metrics";
import { TabContentsProps } from "../TabbedSection";
import { getMetricsForLibraryType } from "../utilities";
import { getDataTest } from "../utils";

export const GEX = ({
  content,
  metrics,
  experimentalDesign,
}: TabContentsProps) => {
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
        {content.barcode_rank_plot && (
          <TitledPlot
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

      <div className="half-plot">
        {content.median_genes_per_cell_plot && (
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
        config={config.GEX.gdna}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />
    </>
  );
};
