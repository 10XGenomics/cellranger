import React from "react";
import Form from "react-bootstrap/Form";

export type ViewType = "visualization" | "table";

/**
 * In a multisample context, determine whether to display data as
 * tables or as 'visualizations' (plates, charts, etc.)
 */
export function View({
  onSelected,
  viewType,
}: {
  onSelected: (viewType: ViewType) => void;
  viewType: ViewType;
}) {
  return (
    <div
      className="summary_border"
      data-test="sample-view"
      style={{
        overflow: "hidden",
        padding: "16px",
      }}
    >
      <h2 style={{ marginTop: 0 }}>View</h2>
      <Form.Group>
        <Form.Check
          checked={viewType === "visualization"}
          data-test="view-visualization"
          id="visualization-check"
          label="Visualization"
          onChange={() => onSelected("visualization")}
          type="radio"
        />
        <Form.Check
          checked={viewType === "table"}
          data-test="view-table"
          id="table-check"
          label="Table"
          onChange={() => onSelected("table")}
          type="radio"
        />
      </Form.Group>
    </div>
  );
}
