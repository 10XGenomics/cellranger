import { describe, expect, test } from "@jest/globals";

import { makeAlertConfig } from "../test-utils";
import { alertLevelColor, getAlertThreshold } from "./alertUtils";

describe("OptionalMultiplexedSummary", () => {
  describe("SampleRollup", () => {
    describe("alertUtils", () => {
      describe("alertLevelColor", () => {
        test("gives the right color per alert level", () => {
          const pairs = [
            ["ERROR", "#D71715"],
            ["WARN", "#F2994A"],
            ["INFO", "#C4CBD5"],
            [undefined, "#C4CBD5"],
          ] as const;

          pairs.forEach(([level, color]) =>
            expect(alertLevelColor(level)).toEqual(color),
          );
        });
      });

      describe("getAlertThreshold", () => {
        test("returns undefined if no alerts", () => {
          const alertConfigs = [];
          const design = { include_introns: false, is_rtl: false };

          expect(getAlertThreshold(alertConfigs, design)).toBeUndefined();
        });

        test("returns the error threshold before others", () => {
          const alertConfigs = [
            makeAlertConfig({ error_threshold: 0, warn_threshold: 1 }),
          ];
          const design = { include_introns: false, is_rtl: false };

          expect(getAlertThreshold(alertConfigs, design)).toEqual({
            threshold: 0,
            thresholdCondition: null,
            thresholdType: "error",
          });
        });

        test("returns the warn threshold if no error threshold", () => {
          const alertConfigs = [makeAlertConfig({ warn_threshold: 5 })];
          const design = { include_introns: false, is_rtl: false };

          expect(getAlertThreshold(alertConfigs, design)).toEqual({
            threshold: 5,
            thresholdCondition: null,
            thresholdType: "warn",
          });
        });

        test("filters alerts based on experimental design", () => {
          const alertConfigs = [
            makeAlertConfig({
              conditions: { include_introns: false, is_rtl: true },
              error_threshold: 1,
            }),
            makeAlertConfig({
              conditions: { include_introns: false, is_rtl: false },
              warn_threshold: 2,
            }),
          ];
          // Both conditions are set, and they are false, so we want alerts which have
          // the same conditions matching the design's conditions, which is the warning alert.
          const design = { include_introns: false, is_rtl: false };

          expect(getAlertThreshold(alertConfigs, design)).toEqual({
            threshold: 2,
            thresholdCondition: null,
            thresholdType: "warn",
          });
        });
      });
    });
  });
});
