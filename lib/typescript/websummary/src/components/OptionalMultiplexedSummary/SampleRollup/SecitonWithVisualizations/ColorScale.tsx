import React from "react";

import type { ColorScale } from "../types";

interface Props {
  bounds: { min: number; max: number };
  colorScale: ColorScale;
  typedTheshold?: {
    threshold: number;
    thresholdType: "error" | "warn";
  };
}
/**
 * We have three labels.
 * - first label: to be left aligned at the start
 * - middle label: center aligned on the center point
 * - last label: right aligned at the end.
 */
function getLabelPosition({
  scaleIndex,
  scaleLength,
  scaleWidth,
}: {
  scaleIndex: number;
  scaleLength: number;
  scaleWidth: number;
}) {
  if (scaleIndex === 0) {
    return { textAnchor: "start", x: 0 };
  } else if (scaleIndex === scaleLength - 1) {
    return { textAnchor: "end", x: scaleWidth };
  } else {
    // Or just scaleWidth * 0.5, but debatably more future proof.
    return { textAnchor: "middle", x: (scaleIndex / scaleLength) * scaleWidth };
  }
}

function ThresholdLine({
  bounds,
  height,
  typedTheshold,
  width,
}: Pick<Props, "bounds"> & {
  typedTheshold: NonNullable<Props["typedTheshold"]>;
  height: number;
  width: number;
}) {
  const { threshold, thresholdType } = typedTheshold;
  const thresholdX =
    ((threshold - bounds.min) / (bounds.max - bounds.min)) * width;
  let stroke = "transparent";
  if (thresholdType == "error") {
    stroke = "#D71715";
  } else if (thresholdType == "warn") {
    stroke = "#F2994A";
  }

  return (
    <g>
      <line
        stroke={stroke}
        strokeWidth={2}
        x1={thresholdX}
        x2={thresholdX}
        y1={0}
        y2={height}
      />
    </g>
  );
}

/**
 * Renders a color scale bar with discrete colors and labels.
 */
export function ColorScale({ bounds, colorScale, typedTheshold }: Props) {
  // Dimensions were determined from the design, which picked them arbitrarily.
  // All layout is done relative to those sizes in the SVG.
  const width = 330;
  const height = 47;
  const topMargin = 14;
  const bottomMargin = 16;

  // A rectangle in the color scale bar where the rectangle has the discrete
  // color.
  const rectHeight = height - topMargin - bottomMargin;
  const rectWidth = width / colorScale.length;

  return (
    <svg
      height={`${height}`}
      viewBox={`0 0 ${width} ${height}`}
      width={`${width}`}
    >
      <g transform={`translate(0, ${topMargin})`}>
        {colorScale.map((c, i) => (
          <rect
            fill={c.color}
            height={rectHeight}
            key={i}
            width={rectWidth}
            x={i * rectWidth}
            y={0}
          />
        ))}
      </g>
      {typedTheshold && (
        <ThresholdLine
          bounds={bounds}
          height={height}
          typedTheshold={typedTheshold}
          width={width}
        />
      )}
      <g transform={`translate(0, ${topMargin + rectHeight})`}>
        {colorScale.map((c, i) =>
          !c.label ? null : (
            <text
              alignmentBaseline="hanging"
              fill="black"
              key={c.label}
              // Text hangs downward, add an additional offset to give spacing from the colors.
              y={4}
              {...getLabelPosition({
                scaleIndex: i,
                scaleLength: colorScale.length,
                scaleWidth: width,
              })}
            >
              {c.label}
            </text>
          ),
        )}
      </g>
    </svg>
  );
}
