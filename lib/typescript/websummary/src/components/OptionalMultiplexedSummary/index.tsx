import React from "react";

import Multiplexed from "./Multiplexed";
import Page from "./Page";
import Singleplexed from "./Singleplexed";
import { MultiWebSummary } from "./types";
import { preprocess } from "./utilities";

const OptionalMultiplexedSummary = (props: MultiWebSummary) => {
  props.per_sample.forEach((sampleWebsummary) => {
    preprocess(sampleWebsummary.data);
  });

  const multiplex = props.experimental_design.multiplexing_method
    ? true
    : false;

  return (
    <Page>
      {multiplex ? <Multiplexed {...props} /> : <Singleplexed {...props} />}
    </Page>
  );
};

export default OptionalMultiplexedSummary;
