/** @jest-environment jsdom */

import { describe, expect, jest, test } from "@jest/globals";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import React from "react";

import {
  makeExperimentalDesign,
  makeMetric,
  makeMultiWebsummaryData,
} from "../test-utils";
import SampleRollup from "./index";

// We don't want anything trying to render Plotly.
// We get it imported via the Beeswarm plots or indirectly
// from the hoc file.
jest.mock("@websummary/components/Plotly.js", () => ({
  react: jest.fn(),
  register: jest.fn(),
}));
jest.mock("@websummary/components/CellRanger/hoc.js", () => ({
  TitledChangeablePlot: jest.fn(),
}));

/**
 * Rendering the component assumes we have some metrics for the library
 * type setup. That is, it will try and find these keys for GEX to group
 * data for the sample.
 * If we provide no mock metrics we will get console warnings.
 * Ideally we would use real metrics data but it isn't important for this
 * test.
 */
function mockMetrics() {
  const keys = [
    "total_singlets",
    "confidently_mapped_reads_in_cells",
    "median_genes_per_singlet",
    "median_umi_per_singlet",
    "reads_confidently_mapped_to_filtered_probe_set",
    "confidently_mapped_to_transcriptome",
  ];

  return keys.map((key) => makeMetric(key, 0, 0));
}

describe("OptionalMultiplexedSummary", () => {
  describe("SampleRollup", () => {
    describe("Default view type", () => {
      test("should default to 'visualization' view in all cases", () => {
        const data = [
          makeMultiWebsummaryData(
            {
              id: "sample 1",
            },
            mockMetrics(),
          ),
        ];
        const experimentalDesign = makeExperimentalDesign({
          is_rtl: true,
        });

        render(
          <SampleRollup data={data} experimentalDesign={experimentalDesign} />,
        );

        // View component should be rendered when plate-based
        const viewComponent = screen.queryByTestId("sample-view");
        expect(viewComponent).toBeTruthy();

        // Visualization should be selected by default
        const vizRadio =
          screen.getByTestId<HTMLInputElement>("view-visualization");
        expect(vizRadio.checked).toBe(true);

        const tableRadio = screen.getByTestId<HTMLInputElement>("view-table");
        expect(tableRadio.checked).toBe(false);

        // Visualization section should be rendered
        const vizSection = screen.queryByTestId("section-with-visualizations");
        expect(vizSection).toBeTruthy();
      });
    });

    test("should render table view when clicked", async () => {
      const user = userEvent.setup();

      const data = [
        makeMultiWebsummaryData(
          {
            id: "sample 1",
            // Non-plate barcode, still renders the visualization by default.
            multiplexing_barcode_ids: ["barcode-1"],
          },
          mockMetrics(),
        ),
      ];
      render(
        <SampleRollup
          data={data}
          experimentalDesign={makeExperimentalDesign({ is_rtl: true })}
        />,
      );

      // View component should not be rendered when not plate-based
      const viewComponent = screen.queryByTestId("sample-view");
      expect(viewComponent).toBeTruthy();

      const tableCheck = screen.getByTestId<HTMLInputElement>("view-table");

      await user.click(tableCheck);

      // We fetch multiple because there are is more than one per tab,
      // but only one is on screen at a time.
      const tableSections = await screen.findAllByTestId("section-with-tables");

      // We expect 3 because we're rendering three tables in this section
      // for GEX by default:
      // 1. 'Cell Calling Quality'
      // 2. 'Mapping Quality (Amongst Reads From Cells Assigned To Sample)'
      // 3. 'UMIs from Genomic DNA',
      expect(tableSections.length).toEqual(3);
    });
  });
});
