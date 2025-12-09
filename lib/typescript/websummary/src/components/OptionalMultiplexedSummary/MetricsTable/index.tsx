import React, { useMemo } from "react";
import { Tooltip } from "react-tooltip";

import { alertIcon } from "../AlertElement";
import { MetricTableConfig } from "../config/types";
import {
  type AlertLevel,
  type ExperimentalDesign,
  type JsonMetricSummary,
  metricIsFraction,
} from "../types";
import { getFormattedValue } from "../utilities";
import { getAlertSeverity, getAlertSeverityFromAlerts } from "../utils";
import { ArrowDown, ArrowUp } from "./Arrow";
import { useRowSort } from "./hooks";
import {
  getMetricHelpText,
  getRenderedHelpText,
  processMultiRowMetrics,
} from "./utils";

interface BaseProps {
  config: MetricTableConfig;
  experimentalDesign: ExperimentalDesign;
  hero: boolean;
}

namespace MetricsTable {
  export interface Props extends BaseProps {
    metrics: JsonMetricSummary[];
  }
}

namespace MetricsTableInner {
  export interface Props extends BaseProps {
    initialSortOrder: string[];
    /** Whether to show inline UI when the tableRows is empty ([]) */
    showEmptyMessaging: boolean;
    tableRows: JsonMetricSummary[][];
  }
}

function tableTooltipId(config: Pick<MetricTableConfig, "title">) {
  return `${config.title}-table`;
}

function alertLevelColor(l: AlertLevel | undefined) {
  const feedbackErrorLightColor = "#FEEDED";
  const feedbackWarningLightColor = "#FEF5ED";

  switch (l) {
    case "ERROR":
      return feedbackErrorLightColor;
    case "WARN":
      return feedbackWarningLightColor;
    default:
      return "transparent";
  }
}

function Header({
  config,
}: {
  config: Pick<MetricTableConfig, "help" | "title">;
}) {
  const { title, help } = config;
  return (
    <div className="header-wrapper">
      {title ? (
        <div>
          <h2
            data-test={"metrics-table-title"}
            data-tooltip-html={help ? getRenderedHelpText(help) : null}
            data-tooltip-id={tableTooltipId(config)}
            id={title}
          >
            {title}
          </h2>
        </div>
      ) : null}
    </div>
  );
}

const buttonStyle = {
  // TODO
  // This is just a default browser style reset. We should probably
  // get this "for free" from bootstap or the other .scss?
  background: "unset",
  border: 0,
  lineHeight: 0,
  padding: 0,
  // And this is to make the button eaiser to click because if it takes
  // the size of the internal SVG, it's just too small to easily click.
  // The scaling value was chosen arbitrarily.
  transform: "scale(1.25)",
};

const buttonsHeightPixels = 13;

function TableHeader({
  config,
  experimentalDesign,
  sortDirection,
  sortKey,
  tableRows,
  onSort,
}: {
  config: MetricTableConfig;
  experimentalDesign: ExperimentalDesign;
  sortDirection: "asc" | "desc";
  sortKey: string;
  tableRows: JsonMetricSummary[][];
  onSort?: (orderKey: string, direction: "asc" | "desc") => void;
}) {
  // Derive the header from the first metrics row because the first row's
  // header value should represent all subsequent rows' contents.
  const [tableHeaders] = tableRows;
  const alignRight = tableHeaders.length > 1 ? "table-right-aligned" : null;

  // We don't bother supporting sorting on tables with <= 2 rows because
  // it adds no benefit to the user.
  const sortable = onSort !== undefined && tableRows.length > 2;

  return (
    <thead
      style={{
        // These styles are necessary to allow the body of the table to
        // scroll while keeping the header in place.
        backgroundColor: "white",
        position: "sticky",
        top: 0,
        zIndex: 1,
      }}
    >
      <tr data-test="metrics-table-header-row">
        {tableHeaders.map((metric) => {
          const ascColor =
            metric.key === sortKey && sortDirection === "asc"
              ? "black"
              : "#8A97AB";
          const descColor =
            metric.key === sortKey && sortDirection === "desc"
              ? "black"
              : "#8A97AB";

          const helpText = getMetricHelpText(
            metric,
            config,
            experimentalDesign,
          );
          return (
            <th
              key={`${metric.key}-header`}
              style={{
                // TODO
                // Reuse the CSS targeting:
                // optional-multiplexed-summary .table td:not(:last-child) here as well
                borderRight: `1px solid #dee2e6`,
                position: "relative",
                verticalAlign: "top",
              }}
            >
              <div
                className={`table-header-container ${alignRight}`}
                style={{
                  // We've got an absolutely positioned child here,
                  // the sorting arrows, not contributing to the size
                  // of the element.
                  // As a result we have to add some padding to avoid
                  // colliding the elements.
                  // Since the height is set on it below at 13px, we
                  // have enough room with that and some extra..
                  paddingBottom: buttonsHeightPixels * 2,
                }}
              >
                <div
                  className={`table-header-content ${alignRight}`}
                  data-test={"metrics-table-column-title"}
                  data-tooltip-html={getRenderedHelpText(helpText)}
                  data-tooltip-id={tableTooltipId(config)}
                >
                  {metric.config.header}
                </div>
                {sortable && (
                  <button
                    data-test="metrics-table-sort"
                    onClick={() =>
                      onSort(
                        metric.key,
                        sortDirection === "asc" ? "desc" : "asc",
                      )
                    }
                    style={{
                      ...buttonStyle,
                      // Positioning matches the padding on the <th> from the style
                      // sheet.
                      // TODO
                      // Either import the value from the sheet or set them both here
                      // for easier maintenance.
                      bottom: "0.75rem",
                      display: "flex",
                      flexFlow: "column",
                      height: `${buttonsHeightPixels}px`,
                      justifyContent: "space-between",
                      position: "absolute",
                      right: "0.75rem",
                    }}
                  >
                    <ArrowUp color={ascColor} />
                    <ArrowDown color={descColor} />
                  </button>
                )}
              </div>
            </th>
          );
        })}
      </tr>
    </thead>
  );
}

function TableRow({
  rowMetrics,
  rowIndex,
}: {
  rowMetrics: JsonMetricSummary[];
  rowIndex: number;
}) {
  const rightAligned = rowMetrics.length > 1 ? "right-aligned" : undefined;
  const alertLevel = getAlertSeverity(
    rowMetrics.map((m) => m.alerts).flat(),
    (a) => a.level,
  );

  return (
    <tr
      data-test="metrics-table-body-row"
      style={{
        backgroundColor: alertLevelColor(alertLevel),
      }}
    >
      {rowMetrics.map((metric) => (
        <td className={rightAligned} key={`${metric.key}-value`}>
          <MetricValue metric={metric} rowIndex={rowIndex} />
        </td>
      ))}
    </tr>
  );
}

const MetricValue = ({
  metric,
  rowIndex,
}: {
  metric: JsonMetricSummary;
  rowIndex: number;
}) => {
  const alertLevel = getAlertSeverityFromAlerts(metric.alerts);
  const formattedValue = getFormattedValue(metric);
  if (!alertLevel) {
    return formattedValue;
  }
  const tooltipId = `${metric.key}-${rowIndex}-alert`;
  return (
    <>
      {/* manually set display to inline-flex to size the div to the content
          this ensures the tooltip is correctly positioned, and align-items
          to center to keep everything on one line if possible */}
      <div
        id={tooltipId}
        style={{ alignItems: "center", display: "inline-flex" }}
      >
        {formattedValue}
        {alertIcon(alertLevel)}
      </div>
      <Tooltip
        anchorSelect={`#${tooltipId}`}
        className={"tooltip table-tooltip"}
        clickable={true}
      >
        {metric.alerts.map((alertSpec, i) => (
          <div key={i}>
            {/* temporary hack to get consistent font size
                the current tooltip-related styles define font size as relative
                to the size of the table content, and this is different for hero
                metrics vs regular metrics
            */}
            <h2 style={{ fontSize: "13.68px", fontWeight: 400 }}>
              {alertIcon(alertLevel)}
              {alertSpec.title} ({alertSpec.formatted_value})
            </h2>
            <p style={{ fontSize: "13.68px", fontWeight: 400 }}>
              {alertSpec.message}
            </p>
          </div>
        ))}
      </Tooltip>
    </>
  );
};

function rowKey(index: number, metric: JsonMetricSummary) {
  if (metricIsFraction(metric)) {
    return `${index}-${metric.key}-${metric.value.numerator}/${metric.value.denominator}`;
  }

  return `${index}-${metric.key}-${metric.value}`;
}

/**
 * Render a table of metrics built from the provided rows.
 */
export const MetricsTableInner = ({
  config,
  hero,
  initialSortOrder,
  showEmptyMessaging,
  tableRows,
  ...rest
}: MetricsTableInner.Props) => {
  const { handleSort, sortDirection, sortKey, sortedTableRows } = useRowSort(
    tableRows,
    initialSortOrder,
  );
  const empty = tableRows.length === 0;

  return (
    <div className={`summary_border ${hero ? "hero-metrics" : null}`}>
      <Header config={config} />
      {empty && showEmptyMessaging && (
        <span data-test="metrics-table-empty">
          No data for current search criteria.
        </span>
      )}
      {!empty && (
        // We can't apply the max-height to the table itself because we don't
        // know the size until all the cells are laid out. Putting it on the
        // parent works though.
        <div style={{ maxHeight: "645px", overflowY: "scroll" }}>
          <table className="table table-hover">
            <TableHeader
              config={config}
              onSort={handleSort}
              sortDirection={sortDirection}
              sortKey={sortKey}
              tableRows={sortedTableRows}
              {...rest}
            />
            <tbody>
              {sortedTableRows.map((rowMetrics, rowIndex) => (
                <TableRow
                  key={rowKey(rowIndex, rowMetrics[0])}
                  rowIndex={rowIndex}
                  rowMetrics={rowMetrics}
                />
              ))}
            </tbody>
          </table>
        </div>
      )}
      <Tooltip
        className={"tooltip table-tooltip"}
        clickable={true}
        data-test={"metrics-table-tooltip"}
        id={tableTooltipId(config)}
      />
    </div>
  );
};

/**
 * Render a table of metrics that may have more than one row of data.
 *
 * The metrics will be grouped and sorted by the configured grouping.
 */
export const MetricsTable = ({
  metrics,
  config,
  ...rest
}: MetricsTable.Props) => {
  const groupedTableRows = useMemo(
    () => processMultiRowMetrics(metrics, config),
    [metrics, config],
  );
  const initialSortOrder =
    config.orderBy ?? config.entries.map(({ key }) => key);

  return (
    <MetricsTableInner
      config={config}
      initialSortOrder={initialSortOrder}
      showEmptyMessaging={false}
      tableRows={groupedTableRows}
      {...rest}
    />
  );
};
