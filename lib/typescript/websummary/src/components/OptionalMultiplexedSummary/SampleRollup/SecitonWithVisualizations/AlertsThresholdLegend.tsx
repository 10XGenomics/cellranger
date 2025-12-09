import React from "react";

export function AlertsThresholdLegend() {
  return (
    <div style={{ marginBottom: 16, marginTop: 16 }}>
      <div>
        <svg height={4} width={90}>
          <line stroke="#D71715" strokeWidth={4} x1={0} x2={90} y1={0} y2={0} />
        </svg>
        <span style={{ marginLeft: 8 }}>Threshold for errors </span>
      </div>
      <div>
        <svg height={4} width={90}>
          <line stroke="#F2994A" strokeWidth={4} x1={0} x2={90} y1={0} y2={0} />
        </svg>
        <span style={{ marginLeft: 8 }}>Threshold for warnings </span>
      </div>
    </div>
  );
}
