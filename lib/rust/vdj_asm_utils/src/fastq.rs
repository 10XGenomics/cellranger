#![deny(missing_docs)]
use std::collections::HashMap;

pub(super) struct CellrangerFastqHeader {
    pub(super) qname: String,
    pub(super) header: String,
    pub(super) tags: HashMap<String, String>,
}

impl CellrangerFastqHeader {
    pub(super) fn new(header: String) -> CellrangerFastqHeader {
        let parts: Vec<_> = header.split("|||").map(String::from).collect();
        let mut tags = HashMap::new();
        if parts.len() > 1 {
            for idx in 1..parts.len() {
                if idx % 2 == 0 {
                    continue;
                }
                tags.insert(parts[idx].to_string(), parts[idx + 1].to_string());
            }
        }

        CellrangerFastqHeader {
            qname: parts[0].to_string(),
            header,
            tags,
        }
    }

    #[cfg(test)]
    fn get_tag(&self, tag: &str) -> Option<&str> {
        self.tags.get(tag).map(String::as_str)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fastq_header() {
        let name = "D000684:775:H3YYTBCXY:1:1101:10541:57223|||BC|||TGGTAAAC|||UR|||AGAGCTGCCA|||UY|||IIIIIIIIII|||TR||||||CB|||TCAGATGCAGGCTCAC-1|||UB|||AGAGCTGCCA";
        let header = CellrangerFastqHeader::new(name.to_string());
        assert_eq!(header.get_tag("BC"), Some("TGGTAAAC"));
        assert_eq!(header.get_tag("TR"), Some(""));
        assert_eq!(header.get_tag("FOO"), None);
    }
}
