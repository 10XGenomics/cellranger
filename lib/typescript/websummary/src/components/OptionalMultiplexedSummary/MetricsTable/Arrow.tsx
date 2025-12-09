import React from "react";

export function ArrowDown({ color }: { color: string }) {
  return (
    <svg
      height="4"
      viewBox="0 0 8 4"
      width="8"
      xmlns="http://www.w3.org/2000/svg"
    >
      <path
        d="M0.542219 4.54759C0.0578604 4.54759 -0.184319 3.95559 0.165496 3.60578L3.60982 0.161453C3.82509 -0.0538176 4.17491 -0.0538176 4.39018 0.161453L7.8345 3.60578C8.18432 3.95559 7.94214 4.54759 7.45778 4.54759H0.542219Z"
        fill={color}
        transform="rotate(180 4 2)"
      />
    </svg>
  );
}

export function ArrowUp({ color }: { color: string }) {
  return (
    <svg
      height="4"
      viewBox="0 0 8 4"
      width="8"
      xmlns="http://www.w3.org/2000/svg"
    >
      <path
        d="M0.542219 4.54759C0.0578604 4.54759 -0.184319 3.95559 0.165496 3.60578L3.60982 0.161453C3.82509 -0.0538176 4.17491 -0.0538176 4.39018 0.161453L7.8345 3.60578C8.18432 3.95559 7.94214 4.54759 7.45778 4.54759H0.542219Z"
        fill={color}
      />
    </svg>
  );
}
