import { useMemo } from "react";

import { AlertLevel } from "../../types";
import {
  getAlertSeverityFromAlerts,
  getMostSevereAlertTitlesFromAlerts,
} from "../../utils";
import { alertLevelColor } from "../alertUtils";
import { PlateDataMap } from "../types";
import { defaultPlateData, parsePlateBarcodeId } from "./plateUtils";
import type { SampleMetricsBarcodes } from "./utils";

/**
 * Barcodes are specific to a sample. I.e. when we look at once sample discretely,
 * it should have barcodes that belong to no other sample. As a result, to get the
 * alerts for a plate, we need to go through each sample and figure out the alerts
 * for them, and then associate the most severe alert with the barcodes that are
 * specific to that sample.
 */
const useBarcodeAlerts = (sampleMetricsBarcodes: SampleMetricsBarcodes) =>
  useMemo(() => {
    type BarcodeSampleAlert = Record<
      string,
      {
        alertLevel: AlertLevel | undefined;
        alertTitles: string[] | undefined;
        sampleId: string;
      }
    >;
    const barcodeToSampleAlert: BarcodeSampleAlert = {};
    sampleMetricsBarcodes.forEach(
      ({ barcodes: sampleBarcodes, metrics: sampleMetrics, sampleId }) => {
        const alerts = sampleMetrics.map((m) => m.alerts).flat();
        const alertLevel = getAlertSeverityFromAlerts(alerts);
        const alertTitles = getMostSevereAlertTitlesFromAlerts(
          alerts,
          alertLevel,
        );

        sampleBarcodes.forEach((barcode) => {
          barcodeToSampleAlert[barcode] = { alertLevel, alertTitles, sampleId };
        });
      },
    );

    return barcodeToSampleAlert;
  }, [sampleMetricsBarcodes]);

/**
 * Using the barcodes across all samples, determine which plate it's associated with,
 * and then also figure out which well it references in the plate.
 * After that, determine whether we have any alerts for any smaples for that
 * well in that plate.
 *
 * This only applies to products where we expect barcodes in the form of:
 * [A-D]-[A-H][0-9]*
 */
export const usePlateAlerts = (
  sampleMetricsBarcodes: SampleMetricsBarcodes,
) => {
  const barcodeToSampleAlert = useBarcodeAlerts(sampleMetricsBarcodes);

  return useMemo(() => {
    const platesWellColors = Object.keys(
      barcodeToSampleAlert,
    ).reduce<PlateDataMap>((plateMap, barcodeId) => {
      // Barcode ids encode the plate, well row and well column in them, so
      // we need to parse them in order to figure out which plate to associate
      // the alert with.
      const chunks = parsePlateBarcodeId(barcodeId);
      if (!chunks) {
        return plateMap;
      }
      const { plate, column, row } = chunks;
      const well = `${row}${column}` as const;

      if (!plateMap[plate]) {
        plateMap[plate] = defaultPlateData();
      }

      const { alertLevel, alertTitles, sampleId } =
        barcodeToSampleAlert[barcodeId];
      const formattedAlertTitles = alertTitles?.join("<br>");
      const tooltip = [
        `Sample barcode ID: ${barcodeId}<br>`,
        `Sample ID: ${sampleId}<br>`,
        alertLevel ? `<br>Alert(s): ${formattedAlertTitles}` : null,
      ]
        .filter(Boolean)
        .join("");

      plateMap[plate][well] = { color: alertLevelColor(alertLevel), tooltip };
      return plateMap;
    }, {});

    return platesWellColors;
  }, [barcodeToSampleAlert]);
};
