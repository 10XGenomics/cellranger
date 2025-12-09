import React from "react";
import { createRoot } from "react-dom/client";
import { MultiWebSummary } from "src/components/OptionalMultiplexedSummary/types";

import { OptionalMultiplexedSummary } from "./src/components/index";
import favicon from "./src/shared/websummary/assets/favicon.ico";
import expandSummary from "./src/shared/websummary/components/CellRanger/PreProcess";
import {
  expand_resources,
  RESOURCE_STRING,
} from "./src/shared/websummary/converters";

declare global {
  const data: MultiWebSummary;
}

// add favicon link tag
const faviconElement = document.createElement("link");
faviconElement.setAttribute("href", favicon);
faviconElement.setAttribute("rel", "icon");
faviconElement.setAttribute("type", "image/x-icon");
document.head.appendChild(faviconElement);

// FIXME: destroy this mechanism
if (Object.prototype.hasOwnProperty.call(data, RESOURCE_STRING)) {
  expand_resources(data);
}

// Do data expansion for ClusteringSelector components.  We feed only one copy
// of the x/y data, but then expand that to apply to every clustering to reduce
// data load.  Note that these expansions fail silently if the data keys are not
// present, so we can run them all without checks here.
// TODO: handle this expansion directly in the ClusteringSelector object
// instead of needing to special-case expansion locations in the data structure.
expandSummary(data);

document.title = data.page_title;

const root = createRoot(document.getElementById("root") as HTMLElement);
root.render(React.createElement(OptionalMultiplexedSummary, data));

// these styles get extracted into a different css file
import "./src/components/styles.scss";
