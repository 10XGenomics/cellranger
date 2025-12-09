/** @jest-environment jsdom */
import { describe, expect, test } from "@jest/globals";
import { render, screen } from "@testing-library/react";
import React from "react";

import { makeMetric } from "../test-utils";
import { ExperimentalDesign } from "../types";
import { MetricsTableInner } from "./index";

const baseConfig = {
  entries: [
    { help: "sample help", key: "sample_id" },
    { help: "singlets help", key: "total_singlets" },
  ],
  help: "table help",
  title: "test table",
};

const baseTableRows = [
  [makeMetric("sample_id", 0, "abc"), makeMetric("total_singlets", 0, "20")],
  [makeMetric("sample_id", 0, "def"), makeMetric("total_singlets", 0, "10")],
  [makeMetric("sample_id", 0, "ghi"), makeMetric("total_singlets", 0, "50")],
];

function getDefaultProps() {
  return {
    config: baseConfig,
    experimentalDesign: {} as ExperimentalDesign,
    hero: false,
    initialSortOrder: ["sample_id"],
    showEmptyMessaging: false,
    tableRows: baseTableRows,
  };
}

describe("OptionalMultiplexedSummary", () => {
  describe("MetricsTable", () => {
    describe("MetricsTableInner", () => {
      test("should render the main title and set up its tooltip", () => {
        render(<MetricsTableInner {...getDefaultProps()} />);

        const title = screen.getByTestId("metrics-table-title");
        expect(title).toHaveAttribute("data-tooltip-id", "test table-table");
        expect(title).toHaveAttribute(
          "data-tooltip-html",
          `<span>${baseConfig.help}</span>`,
        );
      });

      test("should render nothing if the title is omitted from config", () => {
        const props = {
          ...getDefaultProps(),
          config: {
            ...baseConfig,
            title: "",
          },
        };
        const { container } = render(<MetricsTableInner {...props} />);
        // The wrapper div exists, but it should be empty
        expect(
          container.querySelector(".header-wrapper"),
        ).toBeEmptyDOMElement();
      });

      test("should render correct column headers with tooltips", () => {
        render(<MetricsTableInner {...getDefaultProps()} />);

        const columns = screen.getAllByTestId("metrics-table-column-title");
        // Two metrics: sample_id, total_singlets.
        expect(columns).toHaveLength(baseConfig.entries.length);

        for (let i = 0; i < columns.length; i++) {
          const column = columns[i];
          const entry = baseConfig.entries[i];
          expect(column).toHaveAttribute("data-tooltip-id", "test table-table");
          expect(column).toHaveAttribute(
            "data-tooltip-html",
            `<span>${entry.help}</span>`,
          );
          expect(column.textContent).toEqual("Test Header");
        }
      });

      test("should not render sort buttons for tables with 2 or fewer rows", () => {
        const props = {
          ...getDefaultProps(),
          tableRows: baseTableRows.filter((_, i) => i < 2),
        };
        render(<MetricsTableInner {...props} />);

        expect(screen.queryAllByTestId("metrics-table-sort")).toHaveLength(0);
      });

      test("should render sortable rows", () => {
        render(<MetricsTableInner {...getDefaultProps()} />);

        // We expect one set of buttons for each of the 2 columns
        expect(screen.queryAllByTestId("metrics-table-sort")).toHaveLength(2);

        // 1 header row + 3 body rows
        const rows = [
          ...screen.getAllByTestId("metrics-table-header-row"),
          ...screen.getAllByTestId("metrics-table-body-row"),
        ];
        expect(rows).toHaveLength(4);

        // 3 rows * 2 cells/row = 6 cells
        const cells = screen.getAllByRole("cell");
        expect(cells).toHaveLength(6);

        const rowContent = [
          ["Test Header", "Test Header"],
          ...baseTableRows.map((r) => r.map((v) => v.value)),
        ];
        rows.forEach((r, i) => {
          expect(r.textContent).toEqual(rowContent[i].join(""));
        });
      });

      test("shows empty messaging", () => {
        render(
          <MetricsTableInner
            {...getDefaultProps()}
            showEmptyMessaging={true}
            tableRows={[]}
          />,
        );
        expect(screen.queryAllByTestId("metrics-table-empty")).toHaveLength(1);
      });

      test("hides empty messaging", () => {
        render(<MetricsTableInner {...getDefaultProps()} tableRows={[]} />);
        expect(screen.queryAllByTestId("metrics-table-empty")).toHaveLength(0);
      });
    });
  });
});
