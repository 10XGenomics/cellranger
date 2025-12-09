export type MetricTableConfig = {
  title: string;
  help?: string;
  entries: MetricEntryConfig[];
  orderBy?: string[];
};

export type MetricEntryConfig = {
  key: string;
  help: string | ConditionalMetricHelp[];
  category?: "Library" | "Cells";
};

/**
 * Define conditional help text for a metric.
 *
 * The help text will be used if all of the provided conditionals match the
 * experimental conditions.
 */
export type ConditionalMetricHelp = {
  msg: string;
  isRtl?: boolean;
};
