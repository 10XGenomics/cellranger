#![expect(missing_docs)]
use crate::read_pair::ReadHeader;
use crate::read_pair_iter::FastqSet;
use anyhow::{Result, anyhow};
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Clone, Debug)]
pub struct IlluminaHeaderInfo {
    pub instrument: String,
    pub run_number: u32,
    pub flowcell: String,
    pub lane: u32,
}

impl Default for IlluminaHeaderInfo {
    fn default() -> IlluminaHeaderInfo {
        IlluminaHeaderInfo {
            instrument: "unknown_instrument".to_string(),
            run_number: 0,
            flowcell: "unknow_flowcell".to_string(),
            lane: 0,
        }
    }
}

impl FastqSet {
    pub fn get_header_info(&self) -> Result<Option<IlluminaHeaderInfo>> {
        let mut iter = self.clone().read_iter(ReadHeader::KeepId, None)?;

        let read1 = iter
            .next()
            .transpose()?
            .ok_or_else(|| anyhow!("Empty fastq file: {self:?}"))?;

        let header_prefix = std::str::from_utf8(read1.id().expect("missing FASTQ ID"))?;
        let header_parts: Vec<&str> = header_prefix.split(':').collect();

        if header_parts.len() < 4 {
            Ok(None)
        } else {
            let instrument = header_parts[0].to_string();

            // Just bail out with None if we can parse integer fields in header
            let run_number = header_parts[1].parse();
            if run_number.is_err() {
                return Ok(None);
            }
            let run_number = run_number?;

            let flowcell = header_parts[2].to_string();

            let lane = header_parts[3].parse();
            if lane.is_err() {
                return Ok(None);
            }
            let lane = lane?;

            let res = IlluminaHeaderInfo {
                instrument,
                run_number,
                flowcell,
                lane,
            };

            Ok(Some(res))
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::filenames::{LaneSpec, bcl2fastq};
    use crate::read_pair_iter::FastqSet;
    use std::path::Path;

    #[test]
    fn test_parse_fastq_info() -> Result<()> {
        let path = "tests/filenames/bcl2fastq";

        let fqs = bcl2fastq::find_fastqs(Path::new(path), &"Infected".into(), &LaneSpec::Any)?;

        let info = fqs[0].get_header_info()?;

        let correct = Some(IlluminaHeaderInfo {
            instrument: "A00419".to_string(),
            run_number: 42,
            flowcell: "H7CL3DRXX".to_string(),
            lane: 1,
        });

        assert_eq!(info, correct);

        Ok(())
    }

    #[test]
    fn weird_header_info() -> Result<()> {
        let fq = FastqSet {
            r1: "tests/read_pair_iter/weird-header-R1.fastq".to_string(),
            r2: Some("tests/read_pair_iter/weird-header-R2.fastq".to_string()),
            i1: None,
            i2: None,
            r1_interleaved: false,
        };

        let info = fq.get_header_info()?.unwrap();
        println!("info: {info:?}");

        // This is an example of a wierd customer FASTQ adapter from MGI
        // this checks that we cna parse the 4th field correctly if
        // it's the last field before a space or /
        assert_eq!(info.instrument, "3");
        assert_eq!(info.lane, 1000);

        Ok(())
    }

    #[test]
    fn weird_header2_info() -> Result<()> {
        let fq = FastqSet {
            r1: "tests/read_pair_iter/weird-header2-R1.fastq".to_string(),
            r2: Some("tests/read_pair_iter/weird-header2-R2.fastq".to_string()),
            i1: None,
            i2: None,
            r1_interleaved: false,
        };

        let info = fq.get_header_info()?;
        println!("info: {info:?}");
        assert_eq!(info, None);

        Ok(())
    }

    #[test]
    fn weird_header_csi_1376() -> Result<()> {
        let fq = FastqSet {
            r1: "tests/read_pair_iter/csi-1376-R1.fastq".to_string(),
            r2: Some("tests/read_pair_iter/csi-1376-R2.fastq".to_string()),
            i1: None,
            i2: None,
            r1_interleaved: false,
        };

        let info = fq.get_header_info()?;
        println!("info: {info:?}");
        assert_eq!(info, None);

        Ok(())
    }
}
