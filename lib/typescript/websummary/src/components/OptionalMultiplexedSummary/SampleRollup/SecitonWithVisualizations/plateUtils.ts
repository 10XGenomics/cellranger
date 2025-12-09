import {
  WELL_COLUMNS,
  WELL_ROWS,
  WellColumn,
  WellData,
  WellRow,
} from "../../../../shared/websummary/components/CellRanger/Plate";
import { ExperimentalDesign, SampleWebSummary } from "../../types";

export function defaultPlateData() {
  const unusedWellColor = "#FFFFFF";

  const colors = {} as WellData;
  for (const row of WELL_ROWS) {
    for (const column of WELL_COLUMNS) {
      const wellId = `${row}${column}` as const;
      colors[wellId] = { color: unusedWellColor };
    }
  }
  return colors;
}

export function parsePlateBarcodeId(barcode: string) {
  const regex = /^([A-D])-([A-H])([0-9]+)$/;
  const match = barcode.match(regex);

  if (!match) {
    return null;
  }

  const [, plate, row, columnString] = match as [
    string,
    string,
    WellRow,
    string,
  ];
  const column = parseInt(columnString, 10) as WellColumn;
  return { column, plate, row };
}

/**
 * "Plate-based" is when is_rtl is true, but also when a candidate barcode out of
 * all of our barcodes has a plate prefix.
 */
export function isPlateBased(
  experimentalDesign: Pick<ExperimentalDesign, "is_rtl">,
  sampleBarcodes: SampleWebSummary["multiplexing_barcode_ids"],
) {
  const [bc] = sampleBarcodes;
  const chunks = parsePlateBarcodeId(bc);
  if (chunks === null) {
    return false;
  }
  return experimentalDesign.is_rtl;
}
