/** @jest-environment jsdom */

import { describe, expect, test } from "@jest/globals";
import { renderHook } from "@testing-library/react";

import MultiplexSampleConfig from "../../config/multiplex-sample";
import {
  makeAlertConfig,
  makeConfigMetric,
  makeMetric,
} from "../../test-utils";
import { usePlateMetricsData } from "./usePlateMetricsData";
import { groupMetricsWithBarcodes } from "./utils";

describe("OptionalMultiplexedSummary", () => {
  describe("SampleRollup", () => {
    describe("usePlateMetricsData", () => {
      test("generates viz metrics data", () => {
        type Datum = ReturnType<typeof usePlateMetricsData>[number];

        // Ensure we use hero, mapping, gdna so that we can look up the help text.
        const configs = [
          MultiplexSampleConfig.GEX.hero,
          MultiplexSampleConfig.GEX.mapping,
          MultiplexSampleConfig.GEX.gdna,
        ];

        // Required for helpText resolution.
        const experimentalDesign = { include_introns: false, is_rtl: false };

        // Plate and well are determined by regex from the sampleId.
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
                alerts: [
                  // The bounds for this will end up being [1, 3] so the threshold
                  // must be within that range to not be discarded.
                  makeAlertConfig({ error_threshold: 2 }),
                ],
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
          usePlateMetricsData(
            sampleMetricsBarcodes,
            configs,
            experimentalDesign,
            ["total_singlets", "mapped_to_genome", "estimated_gdna_content"],
          ),
        );

        // We expect three metrics (by their key) to be setup:
        // 1. total_singlets
        // 2. mapped_to_genome
        // 3. estimated_gdna_content
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

        function checkBounds(a: Pick<Datum, "bounds">, b: Datum["bounds"]) {
          expect(a.bounds).toEqual(b);
        }
        checkBounds(singlets, { max: 3, min: 1, type: "f64" });
        checkBounds(genome, { max: 1, min: 0, type: "Percent" });
        checkBounds(gdna, { max: 10, min: 1, type: "CountAndPercent" });

        function checkColorScale(
          a: Pick<Datum, "colorScale">,
          b: Datum["colorScale"],
        ) {
          // 32 subdivisions which we could verify but it's really the start, midpoint and end
          // with more significant information, like labels.
          expect(a.colorScale).toHaveLength(32);
          expect(a.colorScale[0]).toEqual(b[0]);
          expect(a.colorScale[16]).toEqual(b[1]);
          expect(a.colorScale[31]).toEqual(b[2]);

          // Constant across all pieces of data since we won't include a label at this position.
          expect(a.colorScale[1]).toEqual({
            color: "#030312",
            label: undefined,
          });
          expect(a.colorScale[15]).toEqual({
            color: "#ae347b",
            label: undefined,
          });
          expect(a.colorScale[30]).toEqual({
            color: "#fceeb0",
            label: undefined,
          });
        }
        checkColorScale(singlets, [
          { color: "#000004", label: "1.0" },
          { color: "#bd3977", label: "2.0" },
          { color: "#fcfdbf", label: "3.0" },
        ]);
        checkColorScale(genome, [
          { color: "#000004", label: "0%" },
          { color: "#bd3977", label: "50%" },
          { color: "#fcfdbf", label: "100%" },
        ]);
        checkColorScale(gdna, [
          { color: "#000004", label: "1.0" },
          { color: "#bd3977", label: "5.5" },
          { color: "#fcfdbf", label: "10.0" },
        ]);

        // Only the singlets metric has a config with alerts where the alert
        // config has a non-null threshold.
        expect(singlets.typedThreshold).toEqual({
          threshold: 2,
          thresholdCondition: null,
          thresholdType: "error",
        });
        expect(genome.typedThreshold).toBeUndefined();
        expect(gdna.typedThreshold).toBeUndefined();
      });
    });
  });
});
