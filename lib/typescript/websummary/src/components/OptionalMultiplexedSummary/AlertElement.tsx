import React from "react";

import { ALERT_LEVEL } from "../../shared/websummary/components/consts";
import ErrorShield from "../../shared/websummary/components/icons/ErrorShield";
import ExclamationTriangle from "../../shared/websummary/components/icons/ExclamationTriangle";
import InfoCircle from "../../shared/websummary/components/icons/InfoCircle";
import { AlertLevel } from "./types";

interface AlertElementProps {
  title: string;
  alertType?: AlertLevel;
}

const typeToIcon = {
  [ALERT_LEVEL.WARN]: <ExclamationTriangle key="warning-icon" />,
  [ALERT_LEVEL.ERROR]: <ErrorShield key="danger-icon" />,
  [ALERT_LEVEL.INFO]: <InfoCircle key="info-icon" />,
};

/**
 * Return the Icon for the given alert level.
 */
export const alertIcon = (level: AlertLevel) => {
  return typeToIcon[level];
};

const AlertElement = (className: string) => {
  const Component = ({ title, alertType }: AlertElementProps) => {
    if (!alertType) {
      return <>{title}</>;
    }

    return (
      <div className={className}>
        {title}
        {alertIcon(alertType)}
      </div>
    );
  };

  return Component;
};

export const TabContents = AlertElement("tab-warning");
