import camelCase from "lodash/camelCase";
import debounce from "lodash/debounce";
import React, { useMemo, useState } from "react";
import Tab from "react-bootstrap/Tab";
import Tabs from "react-bootstrap/Tabs";

import { TabContents } from "../AlertElement";
import config from "../config/multiplex-sample";
import { MetricTableConfig } from "../config/types";
import { TabTitle } from "../consts";
import { VdjSampleTableConfigs } from "../Sample/VDJ";
import { isEmptyTabData } from "../TabbedSection";
import {
  AlertLevel,
  ExperimentalDesign,
  LibraryType,
  MultiWebSummarySampleData,
  TabData,
} from "../types";
import { AlertDisplayControl } from "./AlertDisplayControl";
import {
  useFilteredPerSampleData,
  useLibraryAlertSeverity,
  useSampleHasAlerts,
} from "./hooks";
import { Search } from "./Search";
import { SectionWithVisualizations } from "./SecitonWithVisualizations";
import { SectionWithTables } from "./SectionWithTables";
import { View, ViewType } from "./View";

function generateVdjConfigs(
  vdjConfig: VdjSampleTableConfigs,
  experimentalDesign: ExperimentalDesign,
) {
  const configs = [vdjConfig.hero];
  if (experimentalDesign.multiplexing_method == "Hashtag") {
    configs.push(vdjConfig.sensitivity);
  } else if (experimentalDesign.multiplexing_method == "OH") {
    configs.push(vdjConfig.sensitivityEnrichment);
  }
  configs.push(vdjConfig.pairedClonotypeDiversity);
  return configs;
}

function generateGexConfigs(experimentalDesign: ExperimentalDesign) {
  const configs: MetricTableConfig[] = [config.GEX.hero, config.GEX.mapping];
  if (experimentalDesign.is_rtl) {
    configs.push(config.GEX.gdna);
  }
  return configs;
}

type GetTabData = {
  alertLevel?: AlertLevel;
  beeswarmMetricKeys: string[];
  configs: MetricTableConfig[];
  data?: TabData;
  libraryType: LibraryType;
  title: TabTitle;
};

function getTabData(
  perLibAlertSeverity: ReturnType<typeof useLibraryAlertSeverity>,
  baseSample: MultiWebSummarySampleData,
  experimentalDesign: ExperimentalDesign,
): GetTabData[] {
  return [
    {
      alertLevel: perLibAlertSeverity.gex,
      beeswarmMetricKeys: [
        "total_singlets",
        "confidently_mapped_reads_in_cells",
        "median_genes_per_singlet",
        "median_umi_per_singlet",
        "reads_confidently_mapped_to_filtered_probe_set",
        "confidently_mapped_to_transcriptome",
      ],
      configs: generateGexConfigs(experimentalDesign),
      data: baseSample.data.gex_tab,
      libraryType: "Gene Expression" as const,
      title: TabTitle.GEX,
    },
    {
      alertLevel: perLibAlertSeverity.vdjT,
      beeswarmMetricKeys: [
        "vdj_filtered_bcs",
        "multi_vdj_recombinome_mapped_reads_frac",
        "TRA_vdj_recombinome_mapped_reads_frac",
        "TRB_vdj_recombinome_mapped_reads_frac",
        "multi_vdj_assembly_contig_pair_productive_full_len_bc_frac",
        "vdj_filtered_bcs_cum_frac",
      ],
      configs: generateVdjConfigs(config.VDJ_T, experimentalDesign),
      data: baseSample.data.vdj_t_tab,
      libraryType: "VDJ T" as const,
      title: TabTitle.VDJ_T,
    },
    {
      alertLevel: perLibAlertSeverity.vdjB,
      beeswarmMetricKeys: [
        "vdj_filtered_bcs",
        "multi_vdj_recombinome_mapped_reads_frac",
        "IGH_vdj_recombinome_mapped_reads_frac",
        "IGK_vdj_recombinome_mapped_reads_frac",
        "IGL_vdj_recombinome_mapped_reads_frac",
        "multi_vdj_assembly_contig_pair_productive_full_len_bc_frac",
        "vdj_filtered_bcs_cum_frac",
      ],
      configs: generateVdjConfigs(config.VDJ_B, experimentalDesign),
      data: baseSample.data.vdj_b_tab,
      libraryType: "VDJ B" as const,
      title: TabTitle.VDJ_B,
    },
    {
      alertLevel: perLibAlertSeverity.vdjTGD,
      beeswarmMetricKeys: [
        "vdj_filtered_bcs",
        "multi_vdj_recombinome_mapped_reads_frac",
        "TRG_vdj_recombinome_mapped_reads_frac",
        "TRD_vdj_recombinome_mapped_reads_frac",
        "multi_vdj_assembly_contig_pair_productive_full_len_bc_frac",
        "vdj_filtered_bcs_cum_frac",
      ],
      configs: generateVdjConfigs(config.VDJ_T_GD, experimentalDesign),
      data: baseSample.data.vdj_t_gd_tab,
      libraryType: "VDJ T GD" as const,
      title: TabTitle.VDJ_T_GD,
    },
    {
      alertLevel: perLibAlertSeverity.antibody,
      beeswarmMetricKeys: [
        "total_singlets",
        "fraction_antibody_reads",
        "reads_in_cells",
        "median_umis_per_singlet",
        "antibody_reads_usable_per_cell",
        "fraction_reads_in_aggregate_barcodes",
      ],
      configs: [config.AB.hero],
      data: baseSample.data.antibody_tab,
      libraryType: "Antibody Capture" as const,
      title: TabTitle.AB,
    },
    {
      alertLevel: perLibAlertSeverity.antigen,
      beeswarmMetricKeys: [],
      configs: [config.AG.hero],
      data: baseSample.data.antigen_tab,
      libraryType: "Antigen Capture" as const,
      title: TabTitle.AG,
    },
    {
      alertLevel: perLibAlertSeverity.crispr,
      beeswarmMetricKeys: [
        "total_singlets",
        "reads_in_cells",
        "median_umis_per_singlet",
        "cells_with_one_or_more_protospacers_detected",
        "cells_with_two_or_more_protospacers_detected",
        "fraction_reads_with_putative_protospacer",
        "fraction_guide_reads",
        "guide_reads_usable_per_cell",
      ],
      configs: [config.CRISPR.hero],
      data: baseSample.data.crispr_tab,
      libraryType: "CRISPR Guide Capture" as const,
      title: TabTitle.CRISPR,
    },
    {
      alertLevel: perLibAlertSeverity.customFeature,
      beeswarmMetricKeys: [
        "total_singlets",
        "median_umis_per_singlet",
        "feature_reads_usable_per_cell",
      ],
      configs: [config.Custom.hero],
      data: baseSample.data.custom_feature_tab,
      libraryType: "Custom Feature" as const,
      title: TabTitle.Custom,
    },
  ];
}

const SampleRollup = ({
  data,
  experimentalDesign,
}: {
  data: MultiWebSummarySampleData[];
  experimentalDesign: ExperimentalDesign;
}) => {
  const sampleHasAlerts = useSampleHasAlerts(data);
  const perLibAlertSeverity = useLibraryAlertSeverity(data);
  const atLeastOneSampleHasAlert = sampleHasAlerts.some(Boolean);
  const [showOnlyWithAlerts, setShowOnlyWithAlerts] = useState(
    atLeastOneSampleHasAlert,
  );
  const [sampleIdFilter, setSampleIdFilter] = useState("");
  const sampleData = useFilteredPerSampleData(
    data,
    sampleIdFilter,
    sampleHasAlerts,
    showOnlyWithAlerts,
  );

  const debouncedSetSampleIdFilter = useMemo(
    () =>
      debounce((value: string) => {
        setSampleIdFilter(value);
        // 500ms delay because we don't want to fire this immediately
        // on each character change since it could impact plots where
        // a heavy amount of setup is required.
      }, 500),
    [],
  );

  const handleSearchChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    const cleaned = event.target.value.trim().toLowerCase();
    debouncedSetSampleIdFilter(cleaned);
  };

  // Derive the tabs from the 0th sample because we can assume
  // the first sample's tabs will be consistent across any subsequent
  // sample.
  const [baseSample] = data;
  const tabData = getTabData(
    perLibAlertSeverity,
    baseSample,
    experimentalDesign,
  );

  const [viewType, setViewType] = useState<ViewType>("visualization");

  return (
    <article
      className="summary per-sample"
      data-test="sample-rollup"
      key="sample-rollup"
      style={{
        display: "flex",
        flexFlow: "row wrap",
        gap: "0 32px",
      }}
    >
      <section
        style={{
          display: "flex",
          flex: 1,
          flexFlow: "column nowrap",
          minWidth: "280px",
        }}
      >
        <Search onChange={handleSearchChange} />
        <View onSelected={setViewType} viewType={viewType} />
        <AlertDisplayControl
          checked={showOnlyWithAlerts}
          disabled={!atLeastOneSampleHasAlert}
          onChange={setShowOnlyWithAlerts}
        />
      </section>
      <section style={{ flex: 4, minWidth: "75%" }}>
        <Tabs mountOnEnter={true}>
          {tabData.map((props) =>
            isEmptyTabData(props.data) ? null : (
              <Tab
                data-test={camelCase(props.title)}
                eventKey={props.title}
                key={props.title}
                title={
                  <TabContents
                    alertType={props.alertLevel}
                    title={props.title}
                  />
                }
              >
                {viewType === "table" && (
                  <SectionWithTables
                    configs={props.configs}
                    experimentalDesign={experimentalDesign}
                    libraryType={props.libraryType}
                    sampleData={sampleData}
                  />
                )}
                {viewType === "visualization" && (
                  <SectionWithVisualizations
                    beeswarmMetricKeys={props.beeswarmMetricKeys}
                    configs={props.configs}
                    experimentalDesign={experimentalDesign}
                    libraryType={props.libraryType}
                    sampleData={sampleData}
                  />
                )}
              </Tab>
            ),
          )}
        </Tabs>
      </section>
    </article>
  );
};

export default SampleRollup;
