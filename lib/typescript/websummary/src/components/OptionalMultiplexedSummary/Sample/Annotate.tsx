import React from "react";

import { TabContentsProps } from "../TabbedSection";
import {
  DiffExpTableCardAndTitle,
  VegaLitePlotWithCardAndTitle,
} from "../TitledComponents";
import { getDataTest } from "../utils";

export const Annotate = ({ content }: TabContentsProps) => {
  return (
    <>
      {content.cell_annotation_cell_types_chart && (
        <VegaLitePlotWithCardAndTitle
          dataTest={getDataTest("cell_annotation_cell_types_chart")}
          {...content.cell_annotation_cell_types_chart}
        />
      )}
      <div
        style={{
          display: "flex",
          flexDirection: "row",
          flexWrap: "wrap",
          gap: "2rem",
          justifyContent: "space-between",
        }}
      >
        {content.cell_annotation_violin_plot_chart && (
          <VegaLitePlotWithCardAndTitle
            dataTest={getDataTest("cell_annotation_violin_plot_chart")}
            {...content.cell_annotation_violin_plot_chart}
          />
        )}
        {content.cell_annotation_umap_plot_chart && (
          <VegaLitePlotWithCardAndTitle
            dataTest={getDataTest("cell_annotation_umap_plot_chart")}
            {...content.cell_annotation_umap_plot_chart}
          />
        )}
      </div>
      {content.cell_annotation_diffexp_table && (
        <DiffExpTableCardAndTitle
          dataTest={getDataTest("cell_annotation_diffexp_table")}
          {...content.cell_annotation_diffexp_table}
        />
      )}
    </>
  );
};
