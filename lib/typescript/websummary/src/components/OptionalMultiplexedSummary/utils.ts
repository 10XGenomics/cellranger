import kebabCase from "lodash/kebabCase";

import { ALERT_LEVEL } from "../../shared/websummary/components/consts";
import { AlertLevel, AlertSpec } from "./types";

/**
 * Find the most severe alert type from a collection of things that may contain
 * alerts.
 *
 * Provide the collection, as well as the function to unpack the alert severity
 * for each element in the collection.
 */
export function getAlertSeverity<T>(
  alertContainers: T[],
  getSeverity: (container: T) => AlertLevel | undefined,
) {
  let alertLevel: AlertLevel | undefined = undefined;
  for (const container of alertContainers) {
    alertLevel = getMostSevereAlertLevel(alertLevel, getSeverity(container));
    // If we know it's an error, we're done.
    if (alertLevel === ALERT_LEVEL.ERROR) {
      break;
    }
  }
  return alertLevel;
}

/**
 * Return the most severe alert type from a collection of alerts.
 *
 * This is a special-cased version of getAlertSeverity.
 */
export function getAlertSeverityFromAlerts(alerts: AlertSpec[]) {
  return getAlertSeverity(alerts, (alertSpec) => alertSpec.level);
}

export function getMostSevereAlertTitlesFromAlerts(
  alerts: AlertSpec[],
  severity: AlertLevel | undefined,
) {
  if (alerts.length === 0) {
    return undefined;
  }

  const alertTitles = alerts
    .map((alert) => {
      if (alert.level === severity) {
        return alert.title;
      }
      return null;
    })
    .filter((x): x is string => Boolean(x));

  return alertTitles;
}

const alertLevelIndex = (level: AlertLevel | undefined): number => {
  if (level === undefined) {
    return 0;
  } else if (level === ALERT_LEVEL.INFO) {
    return 1;
  } else if (level === ALERT_LEVEL.WARN) {
    return 2;
  } else if (level === ALERT_LEVEL.ERROR) {
    return 3;
  } else {
    // If we have an unknown alert level, treat it as "WARN"
    return 2;
  }
};

/**
 * Compare two AlertLevels and return the most severe of the two.
 */
const getMostSevereAlertLevel = (
  a: AlertLevel | undefined,
  b: AlertLevel | undefined,
) => {
  if (alertLevelIndex(a) > alertLevelIndex(b)) {
    return a;
  }
  return b;
};

export function getDataTest(key: string) {
  return kebabCase(key);
}
