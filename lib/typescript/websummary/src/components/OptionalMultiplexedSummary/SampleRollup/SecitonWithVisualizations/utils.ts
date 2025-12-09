import {
  JsonMetricSummary,
  LibraryType,
  MultiWebSummarySampleData,
  TabContent,
} from "../../types";

/**
 * Equivalent to a 'zip' over {barcodes|metrics}PerSample so we don't have to
 * rely on 'vectors of scalars' downstream.
 */
export function groupMetricsWithBarcodes(
  sampleIds: string[],
  barcodesPerSample: string[][],
  metricsPerSample: JsonMetricSummary[][],
) {
  if (
    barcodesPerSample.length !== metricsPerSample.length ||
    barcodesPerSample.length !== sampleIds.length
  ) {
    console.warn(
      `groupMetricsWithBarcodes: barcodesPerSample count ${barcodesPerSample.length}, metricsPerSample count ${metricsPerSample.length}, or sample count ${sampleIds.length} do not match`,
    );
  }
  const length = Math.min(barcodesPerSample.length, metricsPerSample.length);
  const sampleMetricsBarcodes: {
    barcodes: string[];
    metrics: JsonMetricSummary[];
    sampleId: string;
  }[] = [];
  for (let i = 0; i < length; i++) {
    sampleMetricsBarcodes.push({
      barcodes: barcodesPerSample[i],
      metrics: metricsPerSample[i],
      sampleId: sampleIds[i],
    });
  }
  return sampleMetricsBarcodes;
}

/**
 * Grouped metrics and barcodes for a sample.
 */
export type SampleMetricsBarcodes = ReturnType<typeof groupMetricsWithBarcodes>;

// Note: there will never be multiplexing data at the sample rollup level.
const TYPE_TO_KEY = {
  "Antibody Capture": "antibody_tab",
  "CRISPR Guide Capture": "crispr_tab",
  "Custom Feature": "custom_feature_tab",
  "Gene Expression": "gex_tab",
  "VDJ B": "vdj_b_tab",
  "VDJ T": "vdj_t_tab",
  "VDJ T GD": "vdj_t_gd_tab",
} as const;

function isTypeToKey(
  libraryType: string,
): libraryType is keyof typeof TYPE_TO_KEY {
  return libraryType in TYPE_TO_KEY;
}

export function contentForLibraryType(
  libraryType: LibraryType,
  sampleData: MultiWebSummarySampleData,
) {
  if (!isTypeToKey(libraryType)) {
    return;
  }

  const tabKey = TYPE_TO_KEY[libraryType];
  return sampleData.data[tabKey]?.content as TabContent;
}
