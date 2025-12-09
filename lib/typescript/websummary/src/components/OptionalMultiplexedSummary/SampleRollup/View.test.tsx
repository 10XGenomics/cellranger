/** @jest-environment jsdom */

import { beforeAll, describe, expect, jest, test } from "@jest/globals";
import { configure, render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import React from "react";

import { View } from "./View";

describe("OptionalMultiplexedSummary", () => {
  describe("SampleRollup", () => {
    describe("View", () => {
      beforeAll(() => {
        configure({ testIdAttribute: "data-test" });
      });

      test("should render the view selection", () => {
        render(<View onSelected={() => {}} viewType="table" />);

        const view = screen.queryByTestId("sample-view");
        expect(view).toBeTruthy();

        const vizCheck =
          screen.getByTestId<HTMLInputElement>("view-visualization");
        expect(vizCheck.checked).toBe(false);

        const tableCheck = screen.getByTestId<HTMLInputElement>("view-table");
        expect(tableCheck.checked).toBe(true);
      });

      test("should handle clicks to switch between table and visualization", async () => {
        const user = userEvent.setup();
        const handleSelected = jest.fn();

        render(<View onSelected={handleSelected} viewType="table" />);

        const tableCheck = screen.getByTestId<HTMLInputElement>("view-table");
        expect(tableCheck.checked).toBe(true);

        const vizCheck =
          screen.getByTestId<HTMLInputElement>("view-visualization");
        expect(vizCheck.checked).toBe(false);

        await user.click(vizCheck);
        expect(handleSelected).toHaveBeenCalledWith("visualization");
        expect(handleSelected).toHaveBeenCalledTimes(1);
      });
    });
  });
});
