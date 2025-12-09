import React from "react";
import Tabs from "react-bootstrap/Tabs";

import { TabTitle } from "../consts";
import { optionalTab } from "../TabbedSection";
import { ExperimentalDesign, MultiWebSummarySampleData } from "../types";
import { getDataTest } from "../utils";
import { Annotate } from "./Annotate";
import { Antibody } from "./Antibody";
import { Antigen } from "./Antigen";
import { CRISPR } from "./CRISPR";
import { Custom } from "./Custom";
import { GEX } from "./GEX";
import { VDJ_B, VDJ_T, VDJ_T_GD } from "./VDJ";

const Sample = ({
  data,
  experimentalDesign,
}: {
  data: MultiWebSummarySampleData;
  experimentalDesign: ExperimentalDesign;
}) => {
  return (
    <article
      className={"summary per-sample"}
      // TODO
      // 'Cells' used to be 'Sample', leaving for backwards compatibility with
      // tests for the moment.
      data-test={getDataTest("sample")}
      key={"Cells"}
    >
      <header
        className="sample-header per-sample-header"
        data-test="summary-header"
      >
        {data.data.description && (
          <div>
            <p
              className="subtitle"
              data-test="summary-header-sample-description"
            >
              {data.data.description}
            </p>
          </div>
        )}
      </header>
      <Tabs mountOnEnter={true}>
        {optionalTab(
          TabTitle.GEX,
          data.data.gex_tab,
          <GEX
            content={data.data.gex_tab?.content}
            experimentalDesign={experimentalDesign}
            metrics={data.metrics}
          />,
        )}
        {optionalTab(
          TabTitle.VDJ_T,
          data.data.vdj_t_tab,
          <VDJ_T
            content={data.data.vdj_t_tab?.content}
            experimentalDesign={experimentalDesign}
            metrics={data.metrics}
          />,
        )}
        {optionalTab(
          TabTitle.VDJ_B,
          data.data.vdj_b_tab,
          <VDJ_B
            content={data.data.vdj_b_tab?.content}
            experimentalDesign={experimentalDesign}
            metrics={data.metrics}
          />,
        )}
        {optionalTab(
          TabTitle.VDJ_T_GD,
          data.data.vdj_t_gd_tab,
          <VDJ_T_GD
            content={data.data.vdj_t_gd_tab?.content}
            experimentalDesign={experimentalDesign}
            metrics={data.metrics}
          />,
        )}
        {optionalTab(
          TabTitle.AB,
          data.data.antibody_tab,
          <Antibody
            content={data.data.antibody_tab?.content}
            experimentalDesign={experimentalDesign}
            metrics={data.metrics}
          />,
        )}
        {optionalTab(
          TabTitle.AG,
          data.data.antigen_tab,
          <Antigen
            content={data.data.antigen_tab?.content}
            experimentalDesign={experimentalDesign}
            metrics={data.metrics}
          />,
        )}
        {optionalTab(
          TabTitle.CRISPR,
          data.data.crispr_tab,
          <CRISPR
            content={data.data.crispr_tab?.content}
            experimentalDesign={experimentalDesign}
            metrics={data.metrics}
          />,
        )}
        {optionalTab(
          TabTitle.Custom,
          data.data.custom_feature_tab,
          <Custom
            content={data.data.custom_feature_tab?.content}
            experimentalDesign={experimentalDesign}
            metrics={data.metrics}
          />,
        )}
        {optionalTab(
          TabTitle.Annotate,
          data.data.cell_annotation_tab,

          <Annotate
            content={data.data.cell_annotation_tab?.content}
            experimentalDesign={experimentalDesign}
            metrics={data.metrics}
          />,
        )}
      </Tabs>
    </article>
  );
};

export default Sample;
