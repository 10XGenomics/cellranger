
use rust_htslib::bam;
use regex::Regex;
use std;
use std::num::ParseIntError;
use std::str::FromStr;
use std::cmp::max;
use std::error::Error;
use std::fmt;
use rust_htslib::bam::Read;
use collections::FxHashMap;
use geometry::Point2D;

const NUM_READS_TO_ESTIMATE_TILE_EXTENTS: u64 = 2000;

#[derive(Debug, PartialEq, Eq)]
pub struct ReadLoc {
    flowcell: String,
    lane: String,
    surface: u32,
    swath: u32,
    tile: u32,
    x: u32,
    y: u32,
}

#[derive(Debug)]
pub enum ReadLocParseError {
    InvalidQname,
    ParseError(ParseIntError),
}

impl fmt::Display for ReadLocParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Invalid read Qname!")
    }
}

impl Error for ReadLocParseError {
    fn description(&self) -> &str {
        match *self {
            ReadLocParseError::InvalidQname => {
                "Read name does not match the expected regex pattern"
            }
            ReadLocParseError::ParseError(ref err) => err.description(),
        }
    }

    fn cause(&self) -> Option<&Error> {
        match *self {
            ReadLocParseError::ParseError(ref err) => Some(err as &Error),
            _ => None,
        }
    }
}

impl From<ParseIntError> for ReadLocParseError {
    fn from(err: ParseIntError) -> Self {
        ReadLocParseError::ParseError(err)
    }
}

impl FromStr for ReadLoc {
    type Err = ReadLocParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {

        let read_name_re = Regex::new("^[a-zA-Z0-9-]+:[a-z0-9-]+:([a-zA-Z0-9-]+):([0-9]):([0-9])([0-9])([0-9]+):([0-9]+):([0-9]+)$").unwrap();

        if let Some(captures) = read_name_re.captures(s) {
            let vals: Vec<_> = captures
                .iter()
                .skip(1)
                .map(|c| c.unwrap().as_str())
                .collect();

            let flowcell = vals[0].to_owned();
            let lane = vals[1].to_owned();
            let surface = vals[2].parse::<u32>()?;
            let swath = vals[3].parse::<u32>()?;
            let tile = vals[4].parse::<u32>()?;
            let x = vals[5].parse::<u32>()?;
            let y = vals[6].parse::<u32>()?;

            Ok(ReadLoc {
                   flowcell,
                   lane,
                   surface,
                   swath,
                   tile,
                   x,
                   y,
               })

        } else {
            Err(ReadLocParseError::InvalidQname)
        }
    }
}


impl ReadLoc {
    pub fn from_bam_record(record: &bam::record::Record) -> Option<Self> {
        let read_name = String::from_utf8(record.qname().to_vec())
            .expect("Failed to parse the record qname as a String");
        ReadLoc::from_str(&read_name).ok()
    }

    pub fn to_flowcell_lane_key(&self) -> String {
        format!("{}_{}", self.flowcell, self.lane)
    }
}


#[derive(Debug, Default)]
pub struct LaneLayout {
    tile_width: u32,
    tile_height: u32,
    num_swaths: u32,
    tiles_per_swath: u32,
    num_reads_observed: u64,
}

impl LaneLayout {
    pub fn new() -> Self {
        LaneLayout::default() // all zeros
    }

    pub fn estimate_extents(&mut self, read_loc: &ReadLoc) {
        self.tile_width = max(self.tile_width, read_loc.x);
        self.tile_height = max(self.tile_height, read_loc.y);
        self.num_swaths = max(self.num_swaths, read_loc.swath);
        self.tiles_per_swath = max(self.tiles_per_swath, read_loc.tile);
        self.num_reads_observed += 1;
    }

    pub fn has_diffusion_duplicates<T>(&self, diffusion_radius: T) -> bool
        where f64: From<T>
    {
        let fc_area = self.tile_width * self.tile_height * self.num_swaths * self.tiles_per_swath;
        let diffusion_area = std::f64::consts::PI * f64::from(diffusion_radius).powi(2);
        // TODO: Robust Divide
        diffusion_area / (fc_area as f64) < 0.025
    }
}


#[derive(Debug, Default)]
pub struct LaneCoordinateSystem {
    lanes: FxHashMap<String, LaneLayout>,
}

impl LaneCoordinateSystem {
    pub fn new() -> Self {
        LaneCoordinateSystem::default()
    }

    pub fn from_bam(bam: &mut bam::Reader) -> Self {
        let start_pos = bam.tell();
        let mut lane_coordinate_sytem = LaneCoordinateSystem::new();
        lane_coordinate_sytem.estimate_tile_extents(bam);
        bam.seek(start_pos).expect("Seeking back to start failed unexpectedly");
        lane_coordinate_sytem
    }

    pub fn estimate_tile_extents(&mut self, bam: &mut bam::Reader) {
        for r in bam.records() {
            let record = r.expect("Failed to retrieve a bam record");

            let read_loc = match ReadLoc::from_bam_record(&record) {
                Some(loc) => loc,
                None => continue, // Pathological case, but possible
            };

            let key = read_loc.to_flowcell_lane_key();
            self.lanes
                .entry(key)
                .or_insert_with(|| LaneLayout::new())
                .estimate_extents(&read_loc);

            if !self.lanes
                    .values()
                    .any(|ref x| x.num_reads_observed < NUM_READS_TO_ESTIMATE_TILE_EXTENTS) {
                break;
            }

        }
    }

    pub fn get_layout_for_read_loc(&self, read_loc: &ReadLoc) -> &LaneLayout {

        let key = read_loc.to_flowcell_lane_key();

        // If we encounter a lane we've never seen before, assume it's like
        // one we've already observed
        match self.lanes.get(&key) {
            Some(layout) => layout,
            None => {
                self.lanes
                    .values()
                    .nth(0)
                    .expect("Expected a non-empty lane layout")
            }
        }

    }

    pub fn convert_to_lane_coords(&self, read_loc: &ReadLoc) -> Option<Point2D<u32>> {

        if self.lanes.is_empty() {
            return None;
        }

        let layout = self.get_layout_for_read_loc(read_loc);

        let lane_x = (read_loc.swath - 1) * layout.tile_width + read_loc.x;
        let lane_y = (read_loc.tile - 1) * layout.tile_height + read_loc.y;
        Some(Point2D::new(lane_x, lane_y))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_read_loc() {
        let read_qname = "D000684:909:HKFV2BCXY:1:1216:5815:28716";
        let read_loc = ReadLoc::from_str(read_qname).unwrap();
        assert_eq!(read_loc,
                   ReadLoc {
                       flowcell: "HKFV2BCXY".into(),
                       lane: "1".into(),
                       surface: 1,
                       swath: 2,
                       tile: 16,
                       x: 5815,
                       y: 28716,
                   });
    }

    #[test]
    fn test_invalid_read_loc() {
        let read_qname = ":909:HKFV2BCXY:1:1216:5815:28716";
        let read_loc = ReadLoc::from_str(read_qname);
        assert!(read_loc.is_err());
    }
}
