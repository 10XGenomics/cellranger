import React, { ComponentType } from "react";

import DifferentialExpressionTable from "../../shared/websummary/components/CellRanger/DifferentialExpressionTable";
import HeaderWithHelp from "../../shared/websummary/components/HeaderWithHelp";
import Plot from "../../shared/websummary/components/Plot";
import VegaLitePlot from "../../shared/websummary/components/VegaLitePlot";

interface HeaderWithHelpProps {
  title: string;
  helpText?: string;
}

/**
 * Used for the cell annotation web summary. Reads a title and inner component around things generated
 * using WithTitle in the main components
 */
const withCardAndTitle = <P extends object>(Component: ComponentType<P>) =>
  function addCardAndBorder(props: { title: HeaderWithHelpProps; inner: P }) {
    return (
      <div className="summary_border">
        <HeaderWithHelp {...props.title} />
        <Component {...props.inner} />
      </div>
    );
  };

export const VegaLitePlotWithCardAndTitle = withCardAndTitle(VegaLitePlot);
export const PlotWithCardAndTitle = withCardAndTitle(Plot);
export const DiffExpTableCardAndTitle = withCardAndTitle(
  DifferentialExpressionTable,
);
