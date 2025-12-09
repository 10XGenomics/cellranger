import React, { CSSProperties, ReactNode } from "react";
import { Tooltip } from "react-tooltip";

import { getRenderedHelpText } from "../../MetricsTable/utils";

/**
 * Common container for a single beeswarm plot or plate metrics that
 * handles title and layout.
 */
export function VisualizationWrapper({
  children,
  helpText,
  style = {},
  title,
  uniqueId,
}: {
  children: ReactNode | ReactNode[];
  helpText?: string;
  style?: CSSProperties;
  title?: string;
  /** Uniquely identify this component for its tooltip. */
  uniqueId: string;
}) {
  const headerDataTooltipId = `${title}-${uniqueId}-tooltip`;

  return (
    <div
      style={{
        display: "flex",
        flexFlow: "column nowrap",
        gap: "16px 0",
        ...style,
      }}
    >
      {title && (
        <h3
          data-tooltip-html={
            helpText ? getRenderedHelpText(helpText) : undefined
          }
          data-tooltip-id={headerDataTooltipId}
          style={{
            // Needed so the element doesn't grow to take up the full width
            // so that tooltips are then centred relative to the element which
            // should be no wider than the text it contains.
            alignSelf: "flex-start",
            fontSize: "16px",
            fontWeight: 600,
            // It can be user supplied data that's too long for the parent.
            // Make it sure it doesn't become larger than the parent and that
            // it wraps when it does.
            maxWidth: "100%",
            overflowWrap: "break-word",
          }}
        >
          {title}
        </h3>
      )}
      {children}
      {helpText && (
        <Tooltip className={"tooltip table-tooltip"} id={headerDataTooltipId} />
      )}
    </div>
  );
}
