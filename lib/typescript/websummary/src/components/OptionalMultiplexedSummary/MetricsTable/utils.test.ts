import { describe, expect, test } from "@jest/globals";

import { makeMetric } from "../test-utils";
import {
  cleanTableStructure,
  placeholderMetric,
  processMultiRowMetrics,
} from "./utils";

describe("OptionalMultiplexedSummary", () => {
  describe("MetricsTable", () => {
    describe("utils", () => {
      test("processMultiRowMetrics", () => {
        // two-row table, two grouping keys
        // same valye for key0, different for key1, should order rows with key0
        // first
        const metrics = [
          makeMetric("key1", 0, "0"),
          makeMetric("key1", 1, "1"),
          makeMetric("key0", 0, "0"),
          makeMetric("key0", 1, "0"),
        ];

        const config = {
          entries: [
            { help: "", key: "key0" },
            { help: "", key: "key1" },
          ],
          orderBy: ["key0", "key1"],
          title: "test",
        };

        expect(processMultiRowMetrics(metrics, config)).toEqual([
          [metrics[2], metrics[0]],
          [metrics[3], metrics[1]],
        ]);

        // should fill placeholder for missing metric
        const sparseMetrics = [metrics[0], metrics[2], metrics[3]];

        const expected = [
          [metrics[2], metrics[0]],
          [metrics[3], placeholderMetric(sparseMetrics[0])],
        ];

        expect(processMultiRowMetrics(sparseMetrics, config)).toEqual(expected);

        // test that this behavior works for pre-collected rows using the inner
        // cleaning function
        const sparseRows = [[metrics[2], metrics[0]], [metrics[3]]];
        expect(cleanTableStructure(config.entries, sparseRows)).toEqual(
          expected,
        );
      });
    });
  });
});
