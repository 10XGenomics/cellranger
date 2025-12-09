import React, { useMemo } from "react";

import { TitledPlot } from "../../../shared/websummary/components/CellRanger/hoc";
import config from "../config/multiplex-sample";
import Metrics from "../Metrics";
import { TabContentsProps } from "../TabbedSection";
import { getMetricsForLibraryType } from "../utilities";
import { getDataTest } from "../utils";

export const CRISPR = ({
  content,
  metrics,
  experimentalDesign,
}: TabContentsProps) => {
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

      <div className="half-plot">
        {content.projection_plot && (
          <TitledPlot
            dataTest={getDataTest("projection_plot")}
            {...content.projection_plot}
          />
        )}
      </div>

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
