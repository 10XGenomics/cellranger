/** @jest-environment jsdom */

import { describe, expect, test } from "@jest/globals";
import { act, renderHook } from "@testing-library/react";

import { makeMetric } from "../test-utils";
import { useRowSort } from "./hooks";

describe("OptionalMultiplexedSummary", () => {
  describe("MetricsTable", () => {
    describe("hooks", () => {
      describe("useRowSort", () => {
        const tableRows = [
          [
            makeMetric("sample_id", 0, "abc"),
            makeMetric("total_singlets", 0, "20"),
          ],
          [
            makeMetric("sample_id", 0, "def"),
            makeMetric("total_singlets", 0, "10"),
          ],
        ];

        test("should use explicit orderBy from config if provided", () => {
          const { result } = renderHook(() =>
            useRowSort(tableRows, ["total_singlets"]),
          );
          // The default ordering has the 20 total singlets before 10. If we order by default
          // ascending, then we would expect the row with 10 to come before the row of 20,
          // which is why our initial sort based on `total_singlets` is not the same as the
          // initial `tableRows`.
          expect(result.current.sortedTableRows).not.toEqual(tableRows);
          expect(result.current.sortedTableRows).toEqual([
            [
              makeMetric("sample_id", 0, "def"),
              makeMetric("total_singlets", 0, "10"),
            ],
            [
              makeMetric("sample_id", 0, "abc"),
              makeMetric("total_singlets", 0, "20"),
            ],
          ]);
        });

        test("should update sort order and direction when handleSort is called", () => {
          const { result } = renderHook(() =>
            useRowSort(tableRows, ["sample_id", "total_singlets"]),
          );
          expect(result.current.sortedTableRows).toEqual(tableRows);

          act(() => {
            // `sample_id` ascending means abc then def, if we switch to descending we expect
            // def then abc.
            result.current.handleSort("sample_id", "desc");
          });

          expect(result.current.sortedTableRows).not.toEqual(tableRows);
          expect(result.current.sortedTableRows).toEqual([
            [
              makeMetric("sample_id", 0, "def"),
              makeMetric("total_singlets", 0, "10"),
            ],
            [
              makeMetric("sample_id", 0, "abc"),
              makeMetric("total_singlets", 0, "20"),
            ],
          ]);
        });
      });
    });
  });
});
