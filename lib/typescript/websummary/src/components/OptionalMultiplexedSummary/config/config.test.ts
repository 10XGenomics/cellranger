import { describe, expect, test } from "@jest/globals";

import MultiplexLibraryConfig from "./multiplex-library";
import MultiplexSampleConfig from "./multiplex-sample";
import { SinglePlexConfig } from "./singleplex";
import { MetricTableConfig } from "./types";

describe("CellRanger", () => {
  describe("OptionalMultiplexedSummary", () => {
    describe("singleplex", () => {
      test("SinglePlexConfig", () => {
        expect(SinglePlexConfig).not.toBeNull();
        validateOrderBy(SinglePlexConfig);
      });
    });

    describe("multiplex-library", () => {
      test("MultiplexLibraryConfig", () => {
        expect(MultiplexLibraryConfig).not.toBeNull();
        validateOrderBy(MultiplexLibraryConfig);
      });
    });

    describe("multiplex-sample", () => {
      test("MultiplexSampleConfig", () => {
        expect(MultiplexSampleConfig).not.toBeNull();
        validateOrderBy(MultiplexSampleConfig);
      });
    });
  });
});

interface MetricTableConfigs {
  [key: string]: MetricTableConfig;
}

interface WsMetricTableConfigs {
  [key: string]: MetricTableConfigs;
}

const validateOrderBy = (cfg: WsMetricTableConfigs) => {
  Object.values(cfg)
    .flatMap((cfgs) => Object.values(cfgs))
    .filter((table) => table.orderBy !== undefined)
    .forEach((table) => {
      const keys = table.entries.map((cfg) => cfg.key);
      for (const orderBy of table.orderBy!) {
        expect(keys).toContain(orderBy);
      }
    });
};
