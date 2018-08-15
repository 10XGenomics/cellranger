//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

// Custom logger

use log::LevelFilter;
use env_logger::Builder;
use chrono::Local;
use std::io::Write;

pub fn init_log() {
    let _ = Builder::new()
        .format(|buf, record| {
            write!(buf, "{} [{}] - {}",
                    Local::now().format("%Y-%m-%dT%H:%M:%S"),
                    record.level(),
                    record.args())
        })
        .filter(None, LevelFilter::Info)
        .init();
}
