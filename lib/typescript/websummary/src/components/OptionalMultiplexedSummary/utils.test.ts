import { describe, expect, test } from "@jest/globals";

import { makeAlert } from "./test-utils";
import {
  getAlertSeverityFromAlerts,
  getMostSevereAlertTitlesFromAlerts,
} from "./utils";

describe("OptionalMultiplexedSummary", () => {
  describe("utils", () => {
    describe("getMostSevereAlertFromAlerts", () => {
      test("returns undefined for empty alerts", () => {
        const result = getMostSevereAlertTitlesFromAlerts([], undefined);
        expect(result).toBeUndefined();
      });

      test("returns the most severe alert", () => {
        const alerts = [
          { ...makeAlert("INFO"), message: "info alert" },
          { ...makeAlert("WARN"), message: "warn alert" },
          {
            ...makeAlert("ERROR"),
            message: "error alert",
            title: "the only error",
          },
        ];
        const severity = getAlertSeverityFromAlerts(alerts);
        const result = getMostSevereAlertTitlesFromAlerts(alerts, severity);
        expect(result).toEqual(["the only error"]);
      });

      test("returns the first alert if all are equally severe", () => {
        const alerts = [
          { ...makeAlert("WARN"), message: "warn alert 1", title: "1" },
          { ...makeAlert("WARN"), message: "warn alert 2", title: "2" },
        ];
        const severity = getAlertSeverityFromAlerts(alerts);
        const result = getMostSevereAlertTitlesFromAlerts(alerts, severity);
        expect(result).toEqual(["1", "2"]);
      });

      test("handles mixed alert levels", () => {
        const alerts = [
          { ...makeAlert("INFO"), message: "info alert" },
          { ...makeAlert("WARN"), message: "warn alert" },
        ];
        const severity = getAlertSeverityFromAlerts(alerts);
        const result = getMostSevereAlertTitlesFromAlerts(alerts, severity);
        expect(result).toEqual(["Test Alert"]);
      });
    });
  });
});
