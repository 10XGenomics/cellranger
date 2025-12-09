import React, { ReactNode } from "react";

import ExternalLink from "../../shared/websummary/components/icons/ExternalLink";

interface Props {
  children: ReactNode;
  headerText?: string;
}

const Page = ({ children }: Props) => {
  return (
    <div className={"page-layout"}>
      <header className={"page-header"}>
        <div className={"branding"}>
          <img
            alt={"10x Genomics"}
            src={
              "https://cdn.10xgenomics.com/image/upload/v1694905936/logos/10x-genomics/10x_Logo_Horizontal_Reverse.svg"
            }
          />
          <h3>Cell Ranger multi QC Report</h3>
        </div>
        <a
          className={"support-button"}
          href={
            "https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/cr-outputs-web-summary"
          }
        >
          <ExternalLink />
          <span>Support Center</span>
        </a>
      </header>
      <section className={`optional-multiplexed-summary content`}>
        {children}
      </section>
    </div>
  );
};

export default Page;
