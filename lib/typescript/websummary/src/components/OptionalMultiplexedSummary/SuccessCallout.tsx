import React from "react";

import SuccessCheck from "../../shared/websummary/components/icons/SuccessCheck";

export interface CalloutProps {
  title?: string;
  body?: string;
}

const Callout = (props: CalloutProps) => {
  return (
    <div className={"callout callout-success"}>
      <SuccessCheck />
      <div className={"callout-message"}>
        <h3>
          {props.title || "No critical issues were identified in this analysis"}
        </h3>
        <p>
          {props.body ||
            "However, we recommend that you examine key metrics below to double-check the results match your expectation."}
        </p>
      </div>
    </div>
  );
};

export default Callout;
