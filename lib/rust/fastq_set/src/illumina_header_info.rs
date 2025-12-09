#![expect(missing_docs)]
use crate::read_pair::{ReadPart, WhichRead};
use crate::read_pair_iter::{InputFastqs, ReadPairIter};
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

impl InputFastqs {
    pub fn get_header_info(&self) -> Result<Option<IlluminaHeaderInfo>> {
        let mut iter = ReadPairIter::from_fastq_files(self)?;

        let read1 = iter
            .next()
            .transpose()?
            .ok_or_else(|| anyhow!("Empty fastq file: {self:?}"))?;

        let header = read1
            .get(WhichRead::R1, ReadPart::Header)
            .ok_or_else(|| anyhow!("No Read1 in FASTQ data"))?;

        let header = std::str::from_utf8(header)?;
        let header_prefix = header.split([' ', '/']).next();
        if header_prefix.is_none() {
            return Ok(None);
        }
        let header_prefix = header_prefix.unwrap();

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
    use crate::InputFastqs;
    use crate::filenames::bcl2fastq::Bcl2FastqDef;
    use crate::filenames::{FindFastqs, LaneSpec};

    #[test]
    fn test_parse_fastq_info() -> Result<()> {
        let path = "tests/filenames/bcl2fastq";

        let query = Bcl2FastqDef {
            fastq_path: path.to_string(),
            sample_name_spec: "Infected".into(),
            lane_spec: LaneSpec::Any,
        };

        let mut fqs = query.find_fastqs()?;
        fqs.sort();

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
        let fq = InputFastqs {
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
        let fq = InputFastqs {
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
        let fq = InputFastqs {
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
