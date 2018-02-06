//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

// Custom logger

use log::LogLevelFilter;
use env_logger::LogBuilder;
use chrono::Local;

pub fn init_log() {
    let _ = LogBuilder::new()
        .format(|record| {
            format!("{} [{}] - {}",
                    Local::now().format("%Y-%m-%dT%H:%M:%S"),
                    record.level(),
                    record.args())
        })
        .filter(None, LogLevelFilter::Info)
        .init();
}
