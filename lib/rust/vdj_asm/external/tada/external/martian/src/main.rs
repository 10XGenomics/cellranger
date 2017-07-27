//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

#[allow(unused_variables)]

use martian::*;
use rustc_serialize::json::Json;
use std::collections::{BTreeMap, HashMap};
use std::thread;
use std::time;
use std::env::args;

#[macro_use]
extern crate log;
extern crate martian;
extern crate rustc_serialize;

pub struct TestStage;


fn call_func(v: f32) -> usize {
    info!("log a message -- call_func");
    panic!("failed in call_func");
}

impl MartianStage for TestStage {
    fn split(&self, args: JsonDict) -> JsonDict {
        info!("Running split!");
        let mut cc =  BTreeMap::new();
        cc.insert("chunks".to_string(), Json::F64(1.0));
        cc
    }

    fn main(&self, args: JsonDict, outs: JsonDict) -> JsonDict {

        thread::sleep(time::Duration::from_millis(120000));
        let mut cc =  BTreeMap::new();
        cc.insert("chunks".to_string(), Json::F64(1.0));
        cc
    }

    fn join(&self, _: JsonDict, _: JsonDict, _: Vec<JsonDict>, chunk_outs: Vec<JsonDict>) -> JsonDict {

        call_func(1.0);
        let mut cc =  BTreeMap::new();
        cc.insert("chunks".to_string(), Json::F64(1.0));
        cc
    }
}


fn main() {

    let mut stage_registry : HashMap<String, Box<MartianStage>> = HashMap::new();
    stage_registry.insert("test".to_string(), Box::new(TestStage));

    let args = args().skip(1).collect();

    // Run the built-in martian adapter
    martian::martian_main(args, stage_registry);
}
