import React from "react";

type ParametersTable = {
  rows: [string, string][];
};

interface TabHeaderProps {
  parametersTable?: ParametersTable;
}

const TabHeader = ({ parametersTable }: TabHeaderProps) => {
  if (!parametersTable) {
    return null;
  }

  const { rows } = parametersTable;

  if (!rows) {
    return null;
  }

  return (
    <header className="tab-header">
      <ul className="parameters-table">
        {rows.map(([title, subtitle]) => (
          <li key={title}>
            <p className="title">{title}</p>
            <p className="subtitle">{subtitle}</p>
          </li>
        ))}
      </ul>
    </header>
  );
};

export default TabHeader;
