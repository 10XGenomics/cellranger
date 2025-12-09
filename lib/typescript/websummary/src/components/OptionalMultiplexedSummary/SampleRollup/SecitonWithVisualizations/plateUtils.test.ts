import { describe, expect, test } from "@jest/globals";

import {
  defaultPlateData,
  isPlateBased,
  parsePlateBarcodeId,
} from "./plateUtils";

describe("OptionalMultiplexedSummary", () => {
  describe("SampleRollup", () => {
    describe("plateUtils", () => {
      describe("defaultPlateData", () => {
        test("sets up default plate colors", () => {
          const data = defaultPlateData();
          expect(data.A1.color).toEqual("#FFFFFF");
          expect(data.H12.color).toEqual("#FFFFFF");
        });
      });

      describe("isPlateBased", () => {
        test("returns true for plate-based barcodes", () => {
          const experimentalDesign = { is_rtl: true };
          const sampleBarcodes = ["A-A1"];
          expect(isPlateBased(experimentalDesign, sampleBarcodes)).toBe(true);
        });

        test("returns false for non plate-based barcodes", () => {
          const experimentalDesign = { is_rtl: true };
          const sampleBarcodes = ["A1", "B2"];
          expect(isPlateBased(experimentalDesign, sampleBarcodes)).toBe(false);
        });

        test("returns false for non-rtl", () => {
          const experimentalDesign = { is_rtl: false };
          const sampleBarcodes = ["A-A1", "B-B2"];
          expect(isPlateBased(experimentalDesign, sampleBarcodes)).toBe(false);
        });
      });

      describe("parsePlateBarcodeId", () => {
        test("parses a valid plate barcode", () => {
          const cases = [
            {
              barcode: "A-A1",
              result: {
                column: 1,
                plate: "A",
                row: "A",
              },
            },
            {
              barcode: "B-A12",
              result: {
                column: 12,
                plate: "B",
                row: "A",
              },
            },
            {
              barcode: "B-C011",
              result: {
                column: 11,
                plate: "B",
                row: "C",
              },
            },
            {
              barcode: "B-C000011",
              result: {
                column: 11,
                plate: "B",
                row: "C",
              },
            },
          ];

          cases.forEach(({ barcode, result }) => {
            expect(parsePlateBarcodeId(barcode)).toEqual(result);
          });
        });

        test("returns null for an invalid plate barcode", () => {
          const barcode = "X-X1";
          expect(parsePlateBarcodeId(barcode)).toBeNull();
        });

        test("returns null for a malformed plate barcode", () => {
          const barcode = "A-1";
          expect(parsePlateBarcodeId(barcode)).toBeNull();
        });

        test("returns null for a no plate barcode", () => {
          const barcode = "A1";
          expect(parsePlateBarcodeId(barcode)).toBeNull();
        });
      });
    });
  });
});
