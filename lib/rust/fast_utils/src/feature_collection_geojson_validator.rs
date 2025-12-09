use anyhow::{Result, bail};
use geojson::{Bbox, Feature, JsonValue};
use serde::de::{SeqAccess, Visitor};
use serde::{Deserialize, Deserializer, de};
use std::collections::HashSet;
use std::fmt;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

const FEATURE_COLLECTION_NAME: &str = "FeatureCollection";

#[derive(Deserialize, Default, Debug, PartialEq)]
/// Matches the FeatureCollection spec from https://docs.rs/geojson/latest/geojson/struct.FeatureCollection.html
/// The `features` field is replaced with `valid_features`
/// which is true if the features were all valid. This is accomplished
/// without loading all features into memory
/// The Feature collection spec also allows a `foreign_members`
/// which is an arbitrary JsonValue. We allow this to be present in the GoeJSON
/// without every deserialising it into memory
pub(crate) struct FeatureCollectionGeoJson {
    /// Type of the collection. Validate makes sure this is FEATURE_COLLECTION_NAME
    r#type: String,

    #[serde(deserialize_with = "deserialize_validate_features")]
    #[serde(rename(deserialize = "features"))]
    /// The field to deserialise features into without loading it into memory
    /// This is always true if the deserialisation was successful
    valid_features: bool,

    #[serde(skip_deserializing)]
    /// Arbitrary metadata allowed by GeoJSON FeatureCollection spec
    foreign_members: Option<JsonValue>,

    /// Optional bounding box allowed by the GeoJSON spec
    bbox: Option<Bbox>,
}

impl FeatureCollectionGeoJson {
    pub(crate) fn from_geojson_path(geojson_path: PathBuf) -> Result<Self> {
        let file = File::open(geojson_path)?;
        let reader = BufReader::new(file);
        let deserialiser = &mut serde_json::Deserializer::from_reader(reader);
        // We will build a set of paths to the unused elements.
        let mut unused_geojson_keys = HashSet::new();
        let ftr_collection: Self = serde_ignored::deserialize(deserialiser, |path| {
            unused_geojson_keys.insert(path.to_string());
        })?;

        let tolerated_extra_keys = HashSet::from(["foreign_members".to_string()]);
        if !unused_geojson_keys.is_subset(&tolerated_extra_keys) {
            bail!(
                "Got extra keys in the geoJson {:?}",
                unused_geojson_keys.difference(&tolerated_extra_keys)
            )
        } else if ftr_collection.r#type != *FEATURE_COLLECTION_NAME {
            bail!(
                "Got GeoJSON type {:?}. Only accept `{}`",
                ftr_collection.r#type,
                FEATURE_COLLECTION_NAME
            )
        }
        Ok(ftr_collection)
    }
}

/// Deserialize the vector of geojson::Features. The entire sequence
/// is not buffered into memory as it would be if we deserialize to Vec<Features>.
///
/// This function just checks that each feature is valid and moves on. Returns true
/// if no invalid feature was encountered.
fn deserialize_validate_features<'de, D>(deserializer: D) -> Result<bool, D::Error>
where
    D: Deserializer<'de>,
{
    struct ValidateVisitor;

    impl<'de> Visitor<'de> for ValidateVisitor {
        type Value = bool;

        fn expecting(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
            formatter.write_str("a sequence of features")
        }

        fn visit_seq<S>(self, mut seq: S) -> Result<bool, S::Error>
        where
            S: SeqAccess<'de>,
        {
            // Go through all features. Deserialise and break if something
            // is invalid.
            while let Some(_feature) = seq.next_element::<Feature>().map_err(|err| {
                de::Error::custom(format!("Came across invalid feature. Error: {err}"))
            })? {}

            Ok(true)
        }
    }

    // Create the visitor and ask the deserializer to drive it. The
    // deserializer will call visitor.visit_seq() if a seq is present in
    // the input data.
    let visitor = ValidateVisitor;
    deserializer.deserialize_seq(visitor)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_good_geojson() -> Result<()> {
        let test_geojson = r#"
        {
  "type": "FeatureCollection",
  "foreign_members": [1, 2],
  "features": [
    {
      "type": "Feature",
      "geometry": {
        "coordinates": [
          [
            [
              1806.711548,
              5072
            ],
            [
              1805.817383,
              5075.817383
            ],
            [
              1805.189941,
              5076.77417
            ],
            [
              1800.947144,
              5077.293213
            ],
            [
              1805.331909,
              5068.667969
            ],
            [
              1806.711548,
              5072
            ]
          ]
        ],
        "type": "Polygon"
      },
      "properties": {
        "cell_id": 1
      }
    },
    {
      "type": "Feature",
      "geometry": {
        "coordinates": [
          [
            [
              2240.851562,
              4936
            ],
            [
              2241.098633,
              4937.412109
            ],
            [
              2234.985596,
              4931.045654
            ],
            [
              2235.734619,
              4931.812012
            ],
            [
              2238.851318,
              4932.758545
            ],
            [
              2239.570068,
              4933.692871
            ],
            [
              2240.267334,
              4934.753418
            ],
            [
              2240.851562,
              4936
            ]
          ]
        ],
        "type": "Polygon"
      },
      "properties": {
        "cell_id": 2
      }
    },
    {
      "type": "Feature",
      "geometry": {
        "coordinates": [
          [
            [
              2304.005859,
              5048
            ],
            [
              2304.482178,
              5049.289307
            ],
            [
              2303.462158,
              5046.913574
            ],
            [
              2304.005859,
              5048
            ]
          ]
        ],
        "type": "Polygon"
      },
      "properties": {
        "cell_id": 3
      }
    }
  ]
}
    "#;
        let mut in_file = NamedTempFile::new()?;
        writeln!(in_file, "{test_geojson}")?;

        let out = FeatureCollectionGeoJson::from_geojson_path(in_file.path().to_path_buf())?;
        assert_eq!(
            out,
            FeatureCollectionGeoJson {
                r#type: "FeatureCollection".to_string(),
                valid_features: true,
                foreign_members: None,
                bbox: None
            }
        );
        Ok(())
    }

    #[test]
    #[should_panic]
    fn test_bad_geojson() {
        let test_geojson = r#"
        {
  "type": "FeatureCollection",
  "foreign_members": [1, 2],
  "features": [
    {
      "type": "Feature",
      "geometry": {
        "coordinates": [
          [
            [
              1806.711548,
              5072
            ],
            [
              "von Neumann",
              5075.817383
            ],
            [
              1805.189941,
              5076.77417
            ],
            [
              1800.947144,
              5077.293213
            ],
            [
              1805.331909,
              5068.667969
            ],
            [
              1806.711548,
              5072
            ]
          ]
        ],
        "type": "Polygon"
      },
      "properties": {
        "cell_id": 1
      }
    },
    {
      "type": "Feature",
      "geometry": {
        "coordinates": [
          [
            [
              2240.851562,
              4936
            ],
            [
              2241.098633,
              4937.412109
            ],
            [
              2234.985596,
              4931.045654
            ],
            [
              2235.734619,
              4931.812012
            ],
            [
              2238.851318,
              4932.758545
            ],
            [
              2239.570068,
              4933.692871
            ],
            [
              2240.267334,
              4934.753418
            ],
            [
              2240.851562,
              4936
            ]
          ]
        ],
        "type": "Polygon"
      },
      "properties": {
        "cell_id": 2
      }
    },
    {
      "type": "Feature",
      "geometry": {
        "coordinates": [
          [
            [
              2304.005859,
              5048
            ],
            [
              2304.482178,
              5049.289307
            ],
            [
              2303.462158,
              5046.913574
            ],
            [
              2304.005859,
              5048
            ]
          ]
        ],
        "type": "Polygon"
      },
      "properties": {
        "cell_id": 3
      }
    }
  ]
}
    "#;
        let mut in_file = NamedTempFile::new().unwrap();
        writeln!(in_file, "{test_geojson}").unwrap();

        let _out =
            FeatureCollectionGeoJson::from_geojson_path(in_file.path().to_path_buf()).unwrap();
    }
}
