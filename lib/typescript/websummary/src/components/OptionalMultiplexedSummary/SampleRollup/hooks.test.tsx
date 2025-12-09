/** @jest-environment jsdom */

import { describe, expect, test } from "@jest/globals";
import { renderHook } from "@testing-library/react";

import { makeMetric } from "../test-utils";
import { useMetricsPerSamplePerLibraryType } from "./hooks";

describe("OptionalMultiplexedSummary", () => {
  describe("SampleRollup", () => {
    describe("hooks", () => {
      describe("useMetricsPerSamplePerLibraryType", () => {
        test("filters down to the specific library type", () => {
          const metrics = [
            [
              {
                ...makeMetric("sampleId", 0, "0"),
                library_type: "Antibody Capture" as const,
              },
            ],
            [
              {
                ...makeMetric("sampleId", 0, "0"),
                library_type: "Antigen Capture" as const,
              },
            ],
            [
              makeMetric("sampleId", 0, "100%"),
              {
                ...makeMetric("sampleId", 0, "0"),
                library_type: "Antigen Capture" as const,
              },
            ],
          ];

          const { result } = renderHook(() =>
            useMetricsPerSamplePerLibraryType(metrics, "Gene Expression"),
          );

          expect(result.current).toEqual([
            // The two groups of metrics get filtered out in full because the library types
            // aren't 'Gene Expression'.
            [],
            [],
            [makeMetric("sampleId", 0, "100%")],
          ]);
        });
      });
    });
  });
});
