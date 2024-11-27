/// Type alias for an owned text, i.e. ``Vec<u8>``.
pub type Text = Vec<u8>;
/// Type alias for a text slice, i.e. ``&[u8]``.
pub type TextSlice<'a> = &'a [u8];

/// Remove a trailing newline from the given string in place.
pub fn trim_newline(s: &mut String) {
    if s.ends_with('\n') {
        s.pop();
    }
}
