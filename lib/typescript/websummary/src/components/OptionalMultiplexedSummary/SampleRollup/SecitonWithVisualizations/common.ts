/**
 * Each seciton can contain a lot of content.
 * Cap them to a max height so that when the content gets
 * really long a user will have to scroll inside each section,
 * while also allowing them to see the overall contents of each
 * section.
 */
export const SECTION_STYLE = {
  maxHeight: 645,
  overflowY: "scroll" as const,
};
