import React from "react";

import { alertLevelColor } from "../alertUtils";
import { PlateGroup } from "./PlateGroup";
import { usePlateAlerts } from "./usePlateAlerts";
import { type SampleMetricsBarcodes } from "./utils";

function Swatch(props: { fill: string }) {
  return (
    <svg height="20" viewBox="0 0 20 20" width="20">
      <circle cx="10" cy="10" r="10" {...props}></circle>
    </svg>
  );
}

export function PlateAlerts({
  sampleMetricsBarcodes,
}: {
  sampleMetricsBarcodes: SampleMetricsBarcodes;
}) {
  const plateAlerts = usePlateAlerts(sampleMetricsBarcodes);
  return (
    <div className={`summary_border`}>
      <h2>Alert</h2>
      <ul
        style={{
          display: "flex",
          gap: "0 24px",
          // TODO
          // Needed to override style defaults. Maybe we can get it from bootstrap
          // CSS instead.
          listStyle: "none",
          padding: 0,
        }}
      >
        <li>
          <Swatch fill={alertLevelColor("ERROR")} /> Error
        </li>
        <li>
          <Swatch fill={alertLevelColor("WARN")} /> Warning
        </li>
        <li>
          <Swatch fill={alertLevelColor(undefined)} /> Normal
        </li>
      </ul>
      <PlateGroup plateDataMap={plateAlerts} />
    </div>
  );
}
