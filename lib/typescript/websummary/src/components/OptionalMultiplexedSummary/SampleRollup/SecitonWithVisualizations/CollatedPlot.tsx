import kebabCase from "lodash/kebabCase";
import React from "react";

import { TitledChangeablePlot } from "../../../../shared/websummary/components/CellRanger/hoc";
import { LibraryType, MultiWebSummarySampleData } from "../../types";
import { SECTION_STYLE } from "./common";
import { collatePlots } from "./plotlyUtils";
import { contentForLibraryType } from "./utils";
import { VisualizationWrapper } from "./VisualizationWrapper";

type PlotKey = "barcode_rank_plot" | "median_genes_per_cell_plot";

export function CollatedPlot({
  batchDataPerSample = false,
  libraryType,
  plotKey,
  plotTitle,
  sampleData,
  legendConfig = "first-plot",
}: {
  /**
   * If provided, we will give a value to the plot by which it chunks up the
   * data to be rendered incrementally.
   */
  batchDataPerSample: boolean;
  libraryType: LibraryType;
  plotKey: PlotKey;
  plotTitle: string;
  sampleData: MultiWebSummarySampleData[];
  legendConfig?: "first-plot" | "none";
}) {
  const plotsToCombine = sampleData
    .map((d) => {
      const plot = contentForLibraryType(libraryType, d)?.[plotKey]?.plot;
      return plot ? { plot, sampleId: d.data.id } : null;
    })
    .filter((x) => !!x);

  const { batchSizes, plot: samplePlot } = collatePlots(
    plotsToCombine,
    legendConfig,
  );
  const plot = {
    ...samplePlot,
    batchSizes: batchDataPerSample ? batchSizes : undefined,
  };

  if (plot.data.length === 0) {
    return null;
  }

  return (
    <div className="summary_border">
      <section style={{ marginTop: 16, ...SECTION_STYLE }}>
        <h2>{plotTitle}</h2>
        <VisualizationWrapper uniqueId={`collated-${kebabCase(plotKey)}`}>
          <TitledChangeablePlot dataTest={kebabCase(plotKey)} plot={plot} />
        </VisualizationWrapper>
      </section>
    </div>
  );
}
