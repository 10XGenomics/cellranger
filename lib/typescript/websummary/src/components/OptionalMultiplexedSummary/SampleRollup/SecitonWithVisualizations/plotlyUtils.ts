import { isEqual } from "lodash";
import cloneDeep from "lodash/cloneDeep";
import type { Config, Layout, PlotData } from "plotly.js";

type Plot = {
  config: Partial<Config>;
  data: ReadonlyArray<Partial<PlotData>>;
  layout: Partial<Layout>;
};

type SamplePlot = {
  plot: Plot;
  sampleId: string;
};

function getCommonAxes(titledPlots: ReadonlyArray<SamplePlot>) {
  // Initialize trackers for the extreme ranges
  let minX = Infinity;
  let maxX = -Infinity;
  let minY = Infinity;
  let maxY = -Infinity;
  titledPlots.forEach(({ plot }) => {
    const xAxisRange = plot.layout?.xaxis?.range;
    if (Array.isArray(xAxisRange) && xAxisRange.length === 2) {
      minX = Math.min(minX, xAxisRange[0]);
      maxX = Math.max(maxX, xAxisRange[1]);
    }

    // Check for explicit yaxis range and update extremes
    const yAxisRange = plot.layout?.yaxis?.range;
    if (Array.isArray(yAxisRange) && yAxisRange.length === 2) {
      minY = Math.min(minY, yAxisRange[0]);
      maxY = Math.max(maxY, yAxisRange[1]);
    }
  });

  return {
    xAxisRange: isFinite(minX) && isFinite(maxX) ? [minX, maxX] : null,
    yAxisRange: isFinite(minY) && isFinite(maxY) ? [minY, maxY] : null,
  };
}
/**
 * Takes a collection of plots which would normally be plotted discretely
 * across different samples and collates them into a single plotlyJS object.
 *
 * For layout and config, the first plot will be used as a template.
 * For this reason, there should be some kind of relationship between the
 * plots, rather than trying to combine arbitrary plots.
 *
 * Also, multi-axis traces are not supported (e.g. xaxis,xaxis2,yaxis,yaxis2,
 * etc.) only single axis traces are (e.g. xaxis, yaxis).
 *
 * @param titledPlots - The plots to collate, each with a sampleId to identify it.
 * @param legendConfig - How to handle legends:
 *  - 'first-plot' means only the first plot's legend is kept / shown.
 *  - 'none' means no legends are kept / shown.
 * @returns A single PlotlyJS object combining all the input plots.:
 */
export function collatePlots(
  titledPlots: ReadonlyArray<SamplePlot>,
  legendConfig: "first-plot" | "none" = "first-plot",
) {
  if (titledPlots.length === 0) {
    return {
      batchSizes: [],
      plot: {
        config: {},
        data: [],
        layout: {},
      },
    };
  }

  const layoutsEqual = titledPlots.every((titledPlot) =>
    isEqual(titledPlot.plot.layout, titledPlots[0].plot.layout),
  );
  if (!layoutsEqual) {
    console.warn(
      "Not all plots have the same layout. This could mean you're trying to collate unrelated plots.",
    );
  }

  // Use the layout and config from the first plot as a template.
  const [first] = titledPlots;
  const finalData: Partial<PlotData>[] = [];
  const finalLayout = cloneDeep(first.plot.layout);
  const batchSizes: number[] = [];

  titledPlots.forEach(({ plot, sampleId }, plotIdx) => {
    batchSizes.push(plot.data.length);

    plot.data.forEach((trace) => {
      // We're going to be changing values on the initially-passed in traces,
      // so we should deep copy them because there may be other references to
      // them.
      const newTrace = cloneDeep(trace);

      // We only keep the legend from the first plot's trace since the legend
      // should be common across plots.
      if (legendConfig === "first-plot") {
        newTrace.showlegend = trace.showlegend && plotIdx === 0;
      } else if (legendConfig === "none") {
        newTrace.showlegend = false;
      }

      // We want to include the sampleId in the hover info no matter where
      // the user hovers over the trace, so we must update it for each chunk
      // in the plot.
      newTrace.hovertemplate = `Sample Id: ${sampleId}${newTrace.text ? "<br><br>%{text}<extra></extra>" : ""}`;
      finalData.push(newTrace);
    });
  });

  const { xAxisRange, yAxisRange } = getCommonAxes(titledPlots);

  if (xAxisRange) {
    finalLayout.xaxis = {
      ...(finalLayout.xaxis ?? {}),
      range: xAxisRange,
    };
  }

  if (yAxisRange) {
    finalLayout.yaxis = {
      ...(finalLayout.yaxis ?? {}),
      range: yAxisRange,
    };
  }

  return {
    batchSizes,
    plot: {
      config: first.plot.config,
      data: finalData,
      layout: finalLayout,
    },
  };
}
