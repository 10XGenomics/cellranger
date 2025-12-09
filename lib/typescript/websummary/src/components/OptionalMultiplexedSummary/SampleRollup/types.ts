import type { WellData } from "../../../shared/websummary/components/CellRanger/Plate";
import { MetricTableConfig } from "../config/types";
import {
  ExperimentalDesign,
  LibraryType,
  MultiWebSummarySampleData,
} from "../types";

export interface MultiSampleTabContentsProps {
  configs: MetricTableConfig[];
  experimentalDesign: ExperimentalDesign;
  libraryType: LibraryType;
  sampleData: MultiWebSummarySampleData[];
}

export type PlateDataMap = Record<string, Partial<WellData>>;

export type ColorScale = { color: string; label?: string }[];
