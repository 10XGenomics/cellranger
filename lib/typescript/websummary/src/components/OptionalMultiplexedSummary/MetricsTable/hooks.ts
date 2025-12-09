import { useMemo, useState } from "react";

import { JsonMetricSummary } from "../types";
import { sortMultiRowMetrics } from "./utils";

type SortDirection = "asc" | "desc";

/**
 * Manages MetricsTable sorting.
 */
export function useRowSort(
  tableRows: JsonMetricSummary[][],
  initialSortOrder: string[],
) {
  const [orderBy, setOrderBy] = useState(initialSortOrder);
  const [sortDirection, setSortDirection] = useState<SortDirection>("asc");

  const sortedTableRows = useMemo(
    () => sortMultiRowMetrics(tableRows, orderBy, sortDirection),
    [tableRows, orderBy, sortDirection],
  );

  const handleSort = (orderKey: string, direction: SortDirection) => {
    setSortDirection(direction);

    const nextOrderBy = [...orderBy];
    const idx = nextOrderBy.indexOf(orderKey);
    if (idx !== -1) {
      nextOrderBy.splice(idx, 1);
    }
    nextOrderBy.unshift(orderKey);
    setOrderBy(nextOrderBy);
  };

  return { handleSort, sortDirection, sortKey: orderBy[0], sortedTableRows };
}
