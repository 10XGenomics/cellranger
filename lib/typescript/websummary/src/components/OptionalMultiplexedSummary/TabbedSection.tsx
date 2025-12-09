import camelCase from "lodash/camelCase";
import React, { ReactNode } from "react";
import { Tab } from "react-bootstrap";

import Alert from "../../shared/websummary/components/Alert";
import InfoCircle from "../../shared/websummary/components/icons/InfoCircle";
import { TabContents } from "./AlertElement";
import TabHeader from "./TabHeader";
import {
  ExperimentalDesign,
  JsonMetricSummary,
  TabContent,
  TabData,
} from "./types";
import { getAlertSeverity } from "./utils";

export interface TabContentsProps {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  content: TabContent & any;
  metrics: JsonMetricSummary[];
  // Note: this data is injected into the tab contents by the websummary React app.
  experimentalDesign: ExperimentalDesign;
}

// NOTE: <Tab> MUST be a direct child of <Tabs>, so we cannot use wrapper components.
// Thus why these are pure typescript functions, NOT functional components.

const tabbedSection = (id: string, tabData: TabData, content: ReactNode) => (
  <Tab
    data-test={camelCase(id)}
    eventKey={id}
    key={id}
    title={
      <TabContents
        alertType={getAlertSeverity(
          tabData.alerts,
          (alertSpec) => alertSpec.level,
        )}
        title={id}
      />
    }
  >
    <div>
      {tabData.content?.disclaimer && (
        <div className="disclaimer">
          <tr>
            <td>
              <InfoCircle />
            </td>
            <td>
              <div
                dangerouslySetInnerHTML={{ __html: tabData.content.disclaimer }}
              />
            </td>
          </tr>
        </div>
      )}
      {tabData.content.parameters_table ? (
        <TabHeader parametersTable={tabData.content.parameters_table} />
      ) : null}
      {tabData.alerts?.length > 0 && (
        <div className="summary_border">
          <Alert alarms={tabData.alerts} />
        </div>
      )}
      {content}
    </div>
  </Tab>
);

export const isEmptyTabData = (data?: TabData) => {
  return !data || Object.keys(data).length === 0;
};

export const optionalTab = (
  id: string,
  tabData: TabData | undefined,
  content: ReactNode,
) => {
  const emptyData = isEmptyTabData(tabData);
  if (emptyData) {
    return null;
  }

  return tabbedSection(id, tabData!, content);
};
