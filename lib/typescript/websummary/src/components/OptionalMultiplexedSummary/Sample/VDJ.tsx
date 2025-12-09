import React, { useMemo } from "react";

import { TitledPlotTable } from "../../../shared/websummary/components/CellRanger/hoc";
import { TitledChangeablePlot } from "../../../shared/websummary/components/CellRanger/hoc";
import config from "../config/multiplex-sample";
import { MetricTableConfig } from "../config/types";
import Metrics from "../Metrics";
import { TabContentsProps } from "../TabbedSection";
import { LibraryType } from "../types";
import { getMetricsForLibraryType } from "../utilities";
import { getDataTest } from "../utils";

export interface VdjSampleTableConfigs {
  hero: MetricTableConfig;
  sensitivity: MetricTableConfig;
  sensitivityEnrichment: MetricTableConfig;
  pairedClonotypeDiversity: MetricTableConfig;
}

const VDJ = (
  tabProps: TabContentsProps,
  libraryType: LibraryType,
  vdjConfig: VdjSampleTableConfigs,
) => {
  const { metrics, content, experimentalDesign } = tabProps;
  const libraryMetrics = useMemo(
    () => getMetricsForLibraryType(libraryType, metrics),
    [libraryType, metrics],
  );
  return (
    <>
      <Metrics
        config={vdjConfig.hero}
        experimentalDesign={experimentalDesign}
        hero={true}
        metrics={libraryMetrics}
      />

      {experimentalDesign.multiplexing_method == "Hashtag" && (
        <Metrics
          config={vdjConfig.sensitivity}
          experimentalDesign={experimentalDesign}
          metrics={libraryMetrics}
        />
      )}

      {experimentalDesign.multiplexing_method == "OH" && (
        <Metrics
          config={vdjConfig.sensitivityEnrichment}
          experimentalDesign={experimentalDesign}
          metrics={libraryMetrics}
        />
      )}

      <Metrics
        config={vdjConfig.pairedClonotypeDiversity}
        experimentalDesign={experimentalDesign}
        metrics={libraryMetrics}
      />

      {content.clonotype_info && (
        <TitledPlotTable
          dataTest={getDataTest("clonotype_info")}
          {...content.clonotype_info}
        />
      )}

      <div className="half-plot">
        {content?.barcode_rank_plot && (
          <TitledChangeablePlot
            dataTest={getDataTest("barcode_rank_plot")}
            {...content.barcode_rank_plot}
          />
        )}
      </div>
    </>
  );
};

export const VDJ_T = (tabProps: TabContentsProps) => {
  return VDJ(tabProps, "VDJ T", config.VDJ_T);
};

export const VDJ_T_GD = (tabProps: TabContentsProps) => {
  return VDJ(tabProps, "VDJ T GD", config.VDJ_T_GD);
};

export const VDJ_B = (tabProps: TabContentsProps) => {
  return VDJ(tabProps, "VDJ B", config.VDJ_B);
};
