import React from "react";
import Form from "react-bootstrap/Form";
import InputGroup from "react-bootstrap/InputGroup";

/**
 * In a multisample context, search through a list of samples by some
 * predicate.
 * Other components are in charge of making sense of what to do with the
 * search results.
 */
export function Search({
  onChange,
}: {
  onChange: React.ChangeEventHandler<HTMLInputElement>;
}) {
  return (
    <div
      className="summary_border"
      data-test="sample-search"
      style={{
        overflow: "hidden",
        padding: "16px",
      }}
    >
      <h2 style={{ marginTop: 0 }}>Search</h2>
      <InputGroup className="search-box">
        <Form.Control
          aria-label="Sample ID"
          onChange={onChange}
          placeholder="Search for a Sample ID"
        />
      </InputGroup>
    </div>
  );
}
