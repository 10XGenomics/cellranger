use std::env;

/// Return the software product name, for example "Cell Ranger" or "Space Ranger"
/// or "Ranger" when the environment variable TENX_PRODUCT is not set.
pub fn get_tenx_product_name() -> String {
    match env::var("TENX_PRODUCT").as_deref() {
        Ok("cellranger") => "Cell Ranger",
        Ok("spaceranger") => "Space Ranger",
        Ok("cellranger-atac") => "Cell Ranger ATAC",
        Ok("cellranger-arc") => "Cell Ranger ARC",
        Ok(s) => panic!("Unknown product: {s}"),
        Err(env::VarError::NotPresent) => "Ranger",
        Err(env::VarError::NotUnicode(s)) => panic!("Invalid Unicode: {s:?}"),
    }
    .to_string()
}

/// Return the software product id, for example "cellranger" or "spaceranger"
/// or "ranger" when the environment variable TENX_PRODUCT is not set.
pub fn get_tenx_product_id() -> String {
    match env::var("TENX_PRODUCT") {
        Ok(p) => p,
        Err(env::VarError::NotPresent) => "ranger".to_string(),
        Err(env::VarError::NotUnicode(s)) => panic!("Invalid Unicode: {s:?}"),
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_get_tenx_product() {
        assert_eq!(env::var("TENX_PRODUCT"), Err(env::VarError::NotPresent));
        assert_eq!(get_tenx_product_name(), "Ranger");
        assert_eq!(get_tenx_product_id(), "ranger");
    }
}
