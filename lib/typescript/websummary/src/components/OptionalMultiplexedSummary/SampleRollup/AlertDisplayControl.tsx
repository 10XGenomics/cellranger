import React from "react";
import Switch from "react-switch";

export function AlertDisplayControl({
  checked,
  disabled,
  onChange,
}: {
  checked: boolean;
  disabled: boolean;
  onChange: (checked: boolean) => void;
}) {
  return (
    <div
      className="summary_border"
      data-test="alert-display-control"
      style={{
        overflow: "hidden",
        padding: "16px",
      }}
    >
      <h2
        style={{
          marginTop: 0,
        }}
      >
        Alerts
      </h2>
      <label
        style={{
          alignItems: "center",
          cursor: disabled ? "not-allowed" : "default",
          display: "flex",
          marginBottom: "0",
        }}
      >
        <Switch
          checked={checked}
          disabled={disabled}
          handleDiameter={16}
          height={18}
          onChange={onChange}
          width={30}
        />
        <span
          style={{
            color: disabled ? "#888" : "#000",
            marginLeft: "8px",
          }}
        >
          Show samples with alerts only
        </span>
      </label>
      {disabled && (
        <div style={{ marginTop: "8px" }}>No samples have an alert</div>
      )}
    </div>
  );
}
