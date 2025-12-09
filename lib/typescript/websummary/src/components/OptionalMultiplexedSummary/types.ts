/* eslint-disable @typescript-eslint/no-explicit-any */
import type { Config, Layout, PlotData } from "plotly.js";

export interface TabData {
  // TODO
  // Try and narrow down content if possible/make it generic so that
  // it matches the specific TabData, e.g. TabDataGEX, TabDataVDJ, etc.
  content: TabContent & any;
  alerts: AlertSpec[];
}

type PlotWithHelp = {
  plot: {
    config: Partial<Config>;
    data: Array<Partial<PlotData>>;
    layout: Partial<Layout>;
  };
  help: {
    helpText: string;
    title: string;
  };
};

export interface TabContent {
  barcode_rank_plot?: PlotWithHelp;
  disclaimer?: string;
  median_genes_per_cell_plot?: PlotWithHelp;
  parameters_table?: string[][];
}

export type AlertLevel = "ERROR" | "WARN" | "INFO";

export type AlertSpec = {
  level: AlertLevel;
  title: string;
  formatted_value: string;
  message: string;
};

export type MultiWebSummary = {
  page_title: string;
  pipeline_version: string;
  library: MultiWebSummaryLibraryData;
  per_sample: MultiWebSummarySampleData[];
  experimental_design: ExperimentalDesign;
  diagnostics: any;
  sample_diagnostics: any;
  _resources: any;
};

export type MultiWebSummaryLibraryData = {
  data: LibraryWebSummary;
  types: string[];
  metrics: JsonMetricSummary[];
};

export type MultiWebSummarySampleData = {
  data: SampleWebSummary;
  metrics: JsonMetricSummary[];
};

export type MultiplexingMethod = "RTL" | "OH" | "CMO" | "Hashtag";

export type ExperimentalDesign = {
  csv: string;
  include_introns: boolean;
  is_rtl: boolean;
  is_barnyard: boolean;
  multiplexing_method: MultiplexingMethod | null;
};

export type LibraryWebSummary = {
  id: string;
  description: string;
  antibody_tab?: TabData;
  antigen_tab?: TabData;
  cmo_tab?: TabData;
  crispr_tab?: TabData;
  custom_feature_tab?: TabData;
  gex_tab?: TabData;
  hashtag_tab?: TabData;
  vdj_b_tab?: TabData;
  vdj_t_tab?: TabData;
  vdj_t_gd_tab?: TabData;
};

export type SampleWebSummary = {
  id: string;
  description: string;
  multiplexing_barcode_ids: string[];
  antibody_tab?: TabData;
  antigen_tab?: TabData;
  cell_annotation_tab?: TabData;
  crispr_tab?: TabData;
  custom_feature_tab?: TabData;
  gex_tab?: TabData;
  vdj_b_tab?: TabData;
  vdj_t_tab?: TabData;
  vdj_t_gd_tab?: TabData;
};

export type Fraction = { numerator: number; denominator: number };

export type JsonMetricSummary = {
  key: string;
  value: string | number | Fraction | null;
  category: string;
  library_type: LibraryType;
  config: MetricConfig;
  grouping_key?: any;
  alerts: AlertSpec[];
};

type JsonMetricSummaryValued<V extends JsonMetricSummary["value"]> = Omit<
  JsonMetricSummary,
  "value"
> & {
  value: V;
};

/** Narrows the JsonMetricSummary value to string */
export function metricIsString(
  jms: JsonMetricSummary,
): jms is JsonMetricSummaryValued<string> {
  return typeof jms.value === "string";
}

/** Narrows the JsonMetricSummary value to number */
export function metricIsNumber(
  jms: JsonMetricSummary,
): jms is JsonMetricSummaryValued<number> {
  return typeof jms.value === "number";
}

/** Narrows the JsonMetricSummary value to a fraction object */
export function metricIsFraction(
  jms: JsonMetricSummary,
): jms is JsonMetricSummaryValued<Fraction> {
  return (
    jms.value != null &&
    typeof jms.value === "object" &&
    "numerator" in jms.value &&
    "denominator" in jms.value
  );
}

/** Narrows the JsonMetricSummary value to null */
export function metricIsEmpty(
  jms: JsonMetricSummary,
): jms is JsonMetricSummaryValued<null> {
  return jms.value === null;
}

export type AlertConfig = {
  conditions: {
    include_introns: boolean | null;
    is_rtl: boolean | null;
  };
  detail: string;
  error_threshold: number | null;
  error_title: string | null;
  if_metric_is: "less_than_or_equal" | "greater_than_or_equal" | null;
  warn_threshold: number | null;
  warn_title: string | null;
};

export type ErrorAlertConfig = {
  error_threshold: number;
  error_title: string;
  warn_threshold: null;
  warn_title: null;
} & Pick<AlertConfig, "conditions" | "detail" | "if_metric_is">;

export type WarnAlertConfig = {
  error_threshold: null;
  error_title: null;
  warn_threshold: number;
  warn_title: string;
} & Pick<AlertConfig, "conditions" | "detail" | "if_metric_is">;

export function alertIsError(alert: AlertConfig): alert is ErrorAlertConfig {
  return alert.error_threshold !== null;
}

export function alertIsWarn(alert: AlertConfig): alert is WarnAlertConfig {
  return alert.warn_threshold !== null;
}

export type MetricConfig = {
  header: string;
  json_key?: string;
  type: MetricValueType;
  help?: string;
  alerts: AlertConfig[];
};

// NOTE: Hashtag isn't a true library type a la Cellranger, but the value is
// set for metrics related to the hashtag subset of antibody libraries.
export type LibraryType =
  | "Gene Expression"
  | "Antibody Capture"
  | "Antigen Capture"
  | "Multiplexing Capture"
  | "CRISPR Guide Capture"
  | "VDJ T"
  | "VDJ T GD"
  | "VDJ B"
  | "Custom Feature"
  | "Hashtag";

/**
 * For a metric with a value, what is the primitive type of that
 * value.
 */
export type MetricValueType =
  | "usize"
  | "FloatAsInt"
  | "f64"
  | "Percent"
  | "CountAndPercent"
  | "String"
  | undefined;
