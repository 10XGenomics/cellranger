import { describe, expect, test } from "@jest/globals";

import { contentForLibraryType } from "./utils";

describe("OptionalMultiplexedSummary", () => {
  describe("SampleRollup", () => {
    describe("utils", () => {
      describe("contentForLibraryType", () => {
        test("returns the correct content for each supported library type", () => {
          const dummySample = {
            data: {
              antibody_tab: { alerts: [], content: {} },
              antigen_tab: { alerts: [], content: {} },
              cmo_tab: { alerts: [], content: {} },
              crispr_tab: { alerts: [], content: {} },
              custom_feature_tab: { alerts: [], content: {} },
              description: "",
              gex_tab: { alerts: [], content: {} },
              id: "",
              multiplexing_barcode_ids: [],
              vdj_b_tab: { alerts: [], content: {} },
              vdj_t_gd_tab: { alerts: [], content: {} },
              vdj_t_tab: { alerts: [], content: {} },
            },
            metrics: [],
          };

          const supported = [
            "Antibody Capture",
            "CRISPR Guide Capture",
            "Custom Feature",
            "Gene Expression",
            "VDJ T",
            "VDJ T GD",
            "VDJ B",
          ] as const;

          const unsupported = [
            "Antigen Capture",
            "Multiplexing Capture",
            "Hashtag",
          ] as const;

          supported.forEach((libraryType) => {
            const content = contentForLibraryType(libraryType, dummySample);
            expect(content).toBeDefined();
          });

          unsupported.forEach((libraryType) => {
            const content = contentForLibraryType(libraryType, dummySample);
            expect(content).toBeUndefined();
          });
        });
      });
    });
  });
});
