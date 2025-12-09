import React from "react";
import { Tooltip } from "react-tooltip";

import {
  Plate,
  PLATE_WELL_TOOLTIP_ID,
} from "../../../../shared/websummary/components/CellRanger/Plate";
import type { PlateDataMap } from "../types";

export function PlateGroup({ plateDataMap }: { plateDataMap: PlateDataMap }) {
  const plateData = Object.entries(plateDataMap);
  plateData.sort(([plate0], [plate1]) => plate0.localeCompare(plate1));
  return (
    <div style={{ display: "flex", gap: "0 16px" }}>
      {plateData.map(([plate, wellData]) => (
        <div
          key={plate}
          style={{ flex: 1, maxWidth: "25%", textAlign: "center" }}
        >
          <Plate wellData={wellData} />
          <h3>{plate}</h3>
        </div>
      ))}
      <Tooltip className={"tooltip table-tooltip"} id={PLATE_WELL_TOOLTIP_ID} />
    </div>
  );
}
