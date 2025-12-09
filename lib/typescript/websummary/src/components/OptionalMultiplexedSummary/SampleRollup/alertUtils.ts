import {
  type AlertConfig,
  alertIsError,
  alertIsWarn,
  type AlertLevel,
  ErrorAlertConfig,
  type ExperimentalDesign,
  WarnAlertConfig,
} from "../types";

export function alertLevelColor(l: AlertLevel | undefined) {
  const feedbackErrorColor = "#D71715";
  const feedbackWarningColor = "#F2994A";
  const steelLightColor = "#C4CBD5";

  switch (l) {
    case "ERROR":
      return feedbackErrorColor;
    case "WARN":
      return feedbackWarningColor;
    default:
      return steelLightColor;
  }
}

/**
 * If the conditions of an alert has the boolean to filter, then compare it against
 * the design. Otherwise, keep the alert.
 */
function filterAlertConfigsToExperiment(
  alerts: AlertConfig[],
  design: Pick<ExperimentalDesign, "is_rtl" | "include_introns">,
) {
  return alerts.filter(
    ({ conditions }) =>
      (conditions.is_rtl !== null
        ? conditions.is_rtl === design.is_rtl
        : true) &&
      (conditions.include_introns !== null
        ? conditions.include_introns === design.include_introns
        : true),
  );
}

export function getAlertThreshold(
  alertConfigs: AlertConfig[],
  experimentDesign: Pick<ExperimentalDesign, "is_rtl" | "include_introns">,
  bounds?: { min: number; max: number },
) {
  const alerts = filterAlertConfigsToExperiment(alertConfigs, experimentDesign);
  if (alerts.length === 0) {
    return undefined;
  }

  // Prefer error_threshold if it exists, otherwise use warn_threshold.
  let alert: ErrorAlertConfig | WarnAlertConfig | undefined;
  for (const a of alerts) {
    if (alertIsError(a)) {
      alert = a;
      // We found an error threshold, so we can stop looking, because
      // either the next alert config will be a warning or it will have
      // a null error_threshold.
      break;
    }

    if (alertIsWarn(a)) {
      alert = a;
      // Keep looking, we may still encounter an error_threshold.
    }
  }

  if (!alert) {
    // No warning or error thresholds found.
    return undefined;
  }

  const threshold = alertIsError(alert)
    ? alert.error_threshold
    : alert.warn_threshold;

  if (bounds) {
    // If the threshold lies outside the bounds for the data we have
    // which we then use for plotting, there's no point in showing it
    // because it won't be visible.
    if (threshold < bounds.min || threshold > bounds.max) {
      return undefined;
    }
  }

  const thresholdType = alertIsError(alert)
    ? ("error" as const)
    : ("warn" as const);
  const thresholdCondition = alert.if_metric_is;

  return { threshold, thresholdCondition, thresholdType };
}
