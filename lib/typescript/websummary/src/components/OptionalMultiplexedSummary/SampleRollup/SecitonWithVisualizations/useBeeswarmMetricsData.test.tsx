/** @jest-environment jsdom */

import { describe, expect, test } from "@jest/globals";
import { renderHook } from "@testing-library/react";

import MultiplexSampleConfig from "../../config/multiplex-sample";
import {
  makeAlertConfig,
  makeConfigMetric,
  makeMetric,
} from "../../test-utils";
import { useBeeswarmMetricsData } from "./useBeeswarmMetricsData";
import { groupMetricsWithBarcodes } from "./utils";

function expectedValue(sampleId: string, value: number) {
  const description = [
    "Sample ID:<br>",
    `${sampleId}<br><br>`,
    "Value:<br>",
    `${value}`,
  ].join("");

  return { description, value };
}

describe("OptionalMultiplexedSummary", () => {
  describe("SampleRollup", () => {
    describe("useBeeswarmMetricsData", () => {
      test("generates plate-based viz metrics data", () => {
        // Ensure we use hero, mapping, gdna so that we can look up the help text.
        const configs = [
          MultiplexSampleConfig.GEX.hero,
          MultiplexSampleConfig.GEX.mapping,
          MultiplexSampleConfig.GEX.gdna,
        ];

        // Required for helpText resolution.
        const experimentalDesign = { include_introns: false, is_rtl: true };
        const barcodesPerSample = [
          ["A-H001"],
          ["A-A12", "B-C01"],
          ["A-H12", "D-E12"],
          ["C-G5"],
        ];

        // Metrics from each GEX table config (hero, mapping, gdna).
        // Note that physical_library_id is going to get dropped because by
        // default it's a string.
        const metricsPerSample = [
          [
            makeConfigMetric(
              {
                alerts: [makeAlertConfig({ error_threshold: 2 })],
                type: "f64",
              },
              "total_singlets",
              0,
              1,
            ),
            makeConfigMetric({ type: "Percent" }, "mapped_to_genome", 0, 0.5),
            makeConfigMetric(
              { type: "CountAndPercent" },
              "estimated_gdna_content",
              0,
              {
                denominator: 10,
                numerator: 1,
              },
            ),
            makeMetric("physical_library_id", 0, "qwerty"),
          ],
          [
            makeConfigMetric({ type: "f64" }, "total_singlets", 0, 2),
            makeConfigMetric({ type: "Percent" }, "mapped_to_genome", 0, 0.1),
            makeConfigMetric(
              { type: "CountAndPercent" },
              "estimated_gdna_content",
              0,
              {
                denominator: 10,
                numerator: 2,
              },
            ),
            makeMetric("physical_library_id", 0, "abcd"),
          ],
          [
            makeConfigMetric({ type: "f64" }, "total_singlets", 0, 3),
            makeConfigMetric({ type: "Percent" }, "mapped_to_genome", 0, 1),
            makeConfigMetric(
              { type: "CountAndPercent" },
              "estimated_gdna_content",
              0,
              {
                denominator: 10,
                numerator: 3,
              },
            ),
            makeMetric("physical_library_id", 0, "asdf"),
          ],
          [],
        ];

        const sampleIds = ["sample-1", "sample-2", "sample-3", "sample-4"];

        const sampleMetricsBarcodes = groupMetricsWithBarcodes(
          sampleIds,
          barcodesPerSample,
          metricsPerSample,
        );

        const { result } = renderHook(() =>
          useBeeswarmMetricsData(
            sampleMetricsBarcodes,
            configs,
            experimentalDesign,
            ["total_singlets", "mapped_to_genome", "estimated_gdna_content"],
          ),
        );

        // We expect 3 metrics (by their key) to be setup:
        // 1. total_singlets
        // 2. mapped_to_genome
        // 3. estimated_gdna_content
        // physical_library_id is not included because it's a string.
        expect(result.current).toHaveLength(3);
        const [singlets, genome, gdna] = result.current;

        expect(singlets.key).toEqual("total_singlets");
        expect(genome.key).toEqual("mapped_to_genome");
        expect(gdna.key).toEqual("estimated_gdna_content");

        expect(singlets.title).toEqual("Test Header");
        expect(genome.title).toEqual("Test Header");
        expect(gdna.title).toEqual("Test Header");

        expect(singlets.helpText).toEqual(
          "Number of cells called from this sample.",
        );
        expect(genome.helpText).toEqual(
          "Fraction of reads that mapped to the genome.",
        );
        expect(gdna.helpText).toEqual(
          "The estimated fraction of filtered UMIs derived " +
            "from genomic DNA based on the discordance between probes targeting " +
            "exon-junction-spanning regions and non-exon-junction-spanning regions.",
        );

        // 3 plates were mapped to: A, B, and D. C wasn't because the only metric
        // that had barcodes mapping onto plate C was physical_library_id,
        expect(singlets.groups).toHaveLength(3);
        expect(genome.groups).toHaveLength(3);
        expect(gdna.groups).toHaveLength(3);

        // For 'total_singlets' we were able to map metrics values
        // from each sample basd on the barcodes prefixed with A-
        // For example, 3 samples, we have A-H001, A-A12, A-H12 and
        // the values are 1, 2, and 3 respectively.
        // We then see in samples 2, 3 we have B- and D- respectively
        // and those metric values are 2 and 3 respectively which is
        // why they mapepd to those plates.
        expect(singlets.groups).toEqual([
          {
            title: "Plate A",
            values: [
              expectedValue("sample-1", 1),
              expectedValue("sample-2", 2),
              expectedValue("sample-3", 3),
            ],
          },
          { title: "Plate B", values: [expectedValue("sample-2", 2)] },
          { title: "Plate D", values: [expectedValue("sample-3", 3)] },
        ]);
        expect(genome.groups).toEqual([
          {
            title: "Plate A",
            values: [
              expectedValue("sample-1", 0.5),
              expectedValue("sample-2", 0.1),
              expectedValue("sample-3", 1),
            ],
          },
          {
            title: "Plate B",
            values: [expectedValue("sample-2", 0.1)],
          },
          { title: "Plate D", values: [expectedValue("sample-3", 1)] },
        ]);
        expect(gdna.groups).toEqual([
          {
            title: "Plate A",
            values: [
              expectedValue("sample-1", 0.1),
              expectedValue("sample-2", 0.2),
              expectedValue("sample-3", 0.3),
            ],
          },
          {
            title: "Plate B",
            values: [expectedValue("sample-2", 0.2)],
          },
          {
            title: "Plate D",
            values: [expectedValue("sample-3", 0.3)],
          },
        ]);

        expect(singlets.typedThreshold).toEqual({
          threshold: 2,
          thresholdCondition: null,
          thresholdType: "error",
        });
        expect(genome.typedThreshold).toBeUndefined();
        expect(gdna.typedThreshold).toBeUndefined();
      });

      test("generates run/non-plate-based viz metrics data", () => {
        // Ensure we use hero, mapping, gdna so that we can look up the help text.
        const configs = [
          MultiplexSampleConfig.GEX.hero,
          MultiplexSampleConfig.GEX.mapping,
          MultiplexSampleConfig.GEX.gdna,
        ];

        // Required for helpText resolution.
        const experimentalDesign = { include_introns: false, is_rtl: true };
        const barcodesPerSample = [["bc-01"], ["bc-02"]];

        // Metrics from each GEX table config (hero, mapping).
        // Note that physical_library_id is going to get dropped because by
        // default it's a string.
        const metricsPerSample = [
          [
            makeConfigMetric(
              {
                alerts: [makeAlertConfig({ error_threshold: 2 })],
                type: "f64",
              },
              "total_singlets",
              0,
              1,
            ),
            makeConfigMetric({ type: "Percent" }, "mapped_to_genome", 0, 0.5),
          ],
          [
            makeConfigMetric({ type: "f64" }, "total_singlets", 0, 2),
            makeConfigMetric({ type: "Percent" }, "mapped_to_genome", 0, 0.1),
          ],
        ];

        const sampleIds = ["sample-1", "sample-2"];

        const sampleMetricsBarcodes = groupMetricsWithBarcodes(
          sampleIds,
          barcodesPerSample,
          metricsPerSample,
        );

        const { result } = renderHook(() =>
          useBeeswarmMetricsData(
            sampleMetricsBarcodes,
            configs,
            experimentalDesign,
            ["total_singlets", "mapped_to_genome"],
          ),
        );

        // We expect 2 metrics (by their key) to be setup:
        // 1. total_singlets
        // 2. mapped_to_genome
        expect(result.current).toHaveLength(2);
        const [singlets, genome] = result.current;

        expect(singlets.key).toEqual("total_singlets");
        expect(genome.key).toEqual("mapped_to_genome");

        expect(singlets.title).toEqual("Test Header");
        expect(genome.title).toEqual("Test Header");

        expect(singlets.helpText).toEqual(
          "Number of cells called from this sample.",
        );
        expect(genome.helpText).toEqual(
          "Fraction of reads that mapped to the genome.",
        );

        // All the barcodes were grouped under the single run.
        expect(singlets.groups).toHaveLength(1);
        expect(genome.groups).toHaveLength(1);

        expect(singlets.groups).toEqual([
          {
            title: "",
            values: [
              expectedValue("sample-1", 1),
              expectedValue("sample-2", 2),
            ],
          },
        ]);
        expect(genome.groups).toEqual([
          {
            title: "",
            values: [
              expectedValue("sample-1", 0.5),
              expectedValue("sample-2", 0.1),
            ],
          },
        ]);

        expect(singlets.typedThreshold).toEqual({
          threshold: 2,
          thresholdCondition: null,
          thresholdType: "error",
        });
        expect(genome.typedThreshold).toBeUndefined();
      });
    });
  });
});
