import {
  afterEach,
  beforeEach,
  describe,
  expect,
  jest,
  test,
} from "@jest/globals";

import { collatePlots } from "./plotlyUtils";

describe("OptionalMultiplexedSummary", () => {
  describe("plotlyUtils", () => {
    describe("collatePlots", () => {
      let warnSpy: jest.SpiedFunction<typeof console.warn>;

      beforeEach(() => {
        warnSpy = jest.spyOn(console, "warn").mockImplementation(jest.fn());
      });

      afterEach(() => {
        warnSpy.mockRestore();
      });

      test("returns empty object for empty input", () => {
        const result = collatePlots([]);
        expect(result).toEqual({
          batchSizes: [],
          plot: {
            config: {},
            data: [],
            layout: {},
          },
        });
      });

      test("collates similar plots", () => {
        const samplePlots = [
          {
            plot: {
              config: { responsive: true },
              data: [
                {
                  name: "trace1",
                  showlegend: true,
                  text: ["a", "b"],
                  x: [1, 2],
                  y: [3, 4],
                },
              ],
              layout: { title: "Plot 1" },
            },
            sampleId: "Plot 1",
          },
          {
            plot: {
              config: { responsive: false },
              data: [
                {
                  name: "trace2",
                  showlegend: true,
                  text: ["a", "b"],
                  x: [5, 6],
                  y: [7, 8],
                },
              ],
              layout: {
                legend: { orientation: "h" as const },
                title: "Plot 2",
              },
            },
            sampleId: "Plot 2",
          },
        ];

        const result = collatePlots(samplePlots);

        // Takes config, layout from plot 1, which is why there's no:
        // * legend or
        // * responsive: false
        expect(result).toEqual({
          batchSizes: [1, 1],
          plot: {
            config: { responsive: true },
            data: [
              {
                hovertemplate:
                  "Sample Id: Plot 1<br><br>%{text}<extra></extra>",
                name: "trace1",
                showlegend: true,
                text: ["a", "b"],
                x: [1, 2],
                y: [3, 4],
              },
              {
                hovertemplate:
                  "Sample Id: Plot 2<br><br>%{text}<extra></extra>",
                name: "trace2",
                // We only show it for traces on the first plot since it
                // should be shared across all traces.
                showlegend: false,
                text: ["a", "b"],
                x: [5, 6],
                y: [7, 8],
              },
            ],
            layout: { title: "Plot 1" }, // Layout from the first plot
          },
        });

        // Warns because of the dissimilar layouts, which could mean mismatched
        // plot types.
        expect(warnSpy).toHaveBeenCalledWith(
          "Not all plots have the same layout. This could mean you're trying to collate unrelated plots.",
        );
      });

      test("collates plots with different axes ranges", () => {
        const samplePlots = [
          {
            plot: {
              config: {},
              data: [{ name: "trace1", text: ["a"], x: [1], y: [10] }],
              layout: {
                title: "First Plot",
                xaxis: { range: [0, 10] },
                yaxis: { range: [0, 20] },
              },
            },
            sampleId: "Plot A",
          },
          {
            plot: {
              config: {},
              data: [{ name: "trace2", text: ["a"], x: [19], y: [28] }],
              layout: {
                title: "Second Plot",
                // Overlaps and extends xaxis
                xaxis: { range: [5, 20] },
                // Overlaps and extends yaxis
                yaxis: { range: [-10, 30] },
              },
            },
            sampleId: "Plot B",
          },
        ];

        const result = collatePlots(samplePlots);

        // The final result should have:
        // 1. The config and base layout from the FIRST plot.
        // 2. Data traces from ALL plots, with modified names.
        // 3. Axis ranges that are the MIN/MAX of ALL plot ranges.
        expect(result).toEqual({
          batchSizes: [1, 1],
          plot: {
            config: {},
            data: [
              {
                hovertemplate:
                  "Sample Id: Plot A<br><br>%{text}<extra></extra>",
                name: "trace1",
                text: ["a"],
                x: [1],
                y: [10],
              },
              {
                hovertemplate:
                  "Sample Id: Plot B<br><br>%{text}<extra></extra>",
                name: "trace2",
                text: ["a"],
                x: [19],
                y: [28],
              },
            ],
            layout: {
              title: "First Plot", // From the first plot
              xaxis: {
                // min of [0, 5] is 0; max of [10, 20] is 20
                range: [0, 20],
              },
              yaxis: {
                // min of [0, -10] is -10; max of [20, 30] is 30
                range: [-10, 30],
              },
            },
          },
        });
      });
    });
  });
});
