import React from "react";

import type { ExperimentalDesign } from "./types";
import { getDataTest } from "./utils";

const CSV = ({ csv }: { csv: string }) => {
  return (
    <section className="csv summary_border">
      <h3>Input CSV</h3>
      <div className="code-wrapper">
        <code data-test={getDataTest("experimental-design-csv")}>{csv}</code>
      </div>
    </section>
  );
};

export default ({ csv }: ExperimentalDesign) => (
  <article
    className="experimental-design"
    data-test={getDataTest("experimental-design")}
    key="Experimental Design"
  >
    {csv && <CSV csv={csv} />}
  </article>
);
