import React from "react";

interface SummaryHeaderProps {
  id: string;
  description: string;
  sampleId: string;
  pipelineVersion: string;
}

export const LibrarySummaryHeader = (props: SummaryHeaderProps) => (
  <header className="sample-header" data-test="summary-header">
    <ul>
      {props.id && (
        <li className="sample-id">
          <p className="title">Run ID</p>
          <p className="subtitle" data-test="summary-header-sample-id">
            {props.id}
          </p>
        </li>
      )}

      {props.sampleId && (
        <li className="sample-id">
          <p className="title">Sample ID</p>
          <p className="subtitle" data-test="summary-header-sample-id">
            {props.sampleId}
          </p>
        </li>
      )}

      {props.description && (
        <li>
          <p className="title">Run Description</p>
          <p className="subtitle" data-test="summary-header-sample-description">
            {props.description}
          </p>
        </li>
      )}

      {props.pipelineVersion && (
        <li>
          <p className="title">Pipeline Version</p>
          <p className="subtitle" data-test="summary-header-pipeline-version">
            {props.pipelineVersion}
          </p>
        </li>
      )}
    </ul>
  </header>
);
