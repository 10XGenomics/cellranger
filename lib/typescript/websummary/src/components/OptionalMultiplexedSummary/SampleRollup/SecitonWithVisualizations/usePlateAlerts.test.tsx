/** @jest-environment jsdom */

import { describe, expect, test } from "@jest/globals";
import { renderHook } from "@testing-library/react";

import { makeAlert, makeMetric } from "../../test-utils";
import { defaultPlateData } from "./plateUtils";
import { usePlateAlerts } from "./usePlateAlerts";
import { groupMetricsWithBarcodes } from "./utils";

describe("OptionalMultiplexedSummary", () => {
  describe("SampleRollup", () => {
    describe("usePlateAlerts", () => {
      test("should map sampleIds with alerts to plate locations", () => {
        // Generated 4 plates: A, B, C, D
        // Where we'll have data for:
        // A: wells H1, A12, H12
        // B: wells C1
        // C: wells G5
        // D: wells E12
        const barcodesPerSample = [
          ["A-H001"],
          ["A-A12", "B-C01"],
          ["A-H12", "D-E12"],
          ["C-G5"],
        ];

        // Where there will be the following alerts for those wells:
        // H1: error
        // A12: warn
        // C1: warn
        // H12: error
        // E12: error
        // G5: info
        const metricsPerSample = [
          [
            {
              ...makeMetric("sampleId", 0, 0),
              alerts: [makeAlert("ERROR")],
            },
          ],
          [
            {
              ...makeMetric("sampleId", 0, 0),
              alerts: [makeAlert("WARN")],
            },
          ],
          [
            {
              ...makeMetric("sampleId", 0, 0),
              alerts: [makeAlert("ERROR")],
            },
          ],
          [
            {
              ...makeMetric("sampleId", 0, 0),
              alerts: [makeAlert("INFO")],
            },
          ],
        ];

        const sampleIds = ["sample-1", "sample-2", "sample-3", "sample-4"];

        const sampleMetricsBarcodes = groupMetricsWithBarcodes(
          sampleIds,
          barcodesPerSample,
          metricsPerSample,
        );

        const { result } = renderHook(() =>
          usePlateAlerts(sampleMetricsBarcodes),
        );

        expect(result.current).toEqual({
          A: {
            ...defaultPlateData(),
            A12: {
              color: "#F2994A",
              tooltip:
                "Sample barcode ID: A-A12<br>Sample ID: sample-2<br><br>Alert(s): Test Alert",
            },
            H1: {
              color: "#D71715",
              tooltip:
                "Sample barcode ID: A-H001<br>Sample ID: sample-1<br><br>Alert(s): Test Alert",
            },
            H12: {
              color: "#D71715",
              tooltip:
                "Sample barcode ID: A-H12<br>Sample ID: sample-3<br><br>Alert(s): Test Alert",
            },
          },
          B: {
            ...defaultPlateData(),
            C1: {
              color: "#F2994A",
              tooltip:
                "Sample barcode ID: B-C01<br>Sample ID: sample-2<br><br>Alert(s): Test Alert",
            },
          },
          C: {
            ...defaultPlateData(),
            // "INFO" alert in this case. We could still show the tooltip, but the
            // color is grey (rather than the default white) because the well is being
            // used but has no alerts, compared to wells with alerts or unused wells.
            G5: {
              color: "#C4CBD5",
              tooltip:
                "Sample barcode ID: C-G5<br>Sample ID: sample-4<br><br>Alert(s): Test Alert",
            },
          },
          D: {
            ...defaultPlateData(),
            E12: {
              color: "#D71715",
              tooltip:
                "Sample barcode ID: D-E12<br>Sample ID: sample-3<br><br>Alert(s): Test Alert",
            },
          },
        });
      });
    });
  });
});
