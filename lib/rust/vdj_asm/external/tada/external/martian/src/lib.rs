//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

//! Martian adapter for Rust code
//!
//! WIP.
//! TODOs: error handling (trap panics?), heartbeat, memory usage monitor.


extern crate libc;
extern crate chrono;
extern crate rustc_serialize;
extern crate backtrace;

use std::{thread, time};
use std::cmp::min;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use backtrace::Backtrace;

#[macro_use]
extern crate log;
extern crate fern;

use std::fs::{File, OpenOptions, rename};
use std::io::{Read, Write};
use std::env;
use std::panic;
use std::collections::{BTreeMap, HashSet, HashMap};
use std::path::PathBuf;
use std::cmp::max;

use chrono::*;

use libc::{timeval, rusage, getrusage, getpid, rlimit, c_ulong};
use rustc_serialize::{Encodable, Decodable};
use rustc_serialize::json::{self, Json, ParserError, ToJson};


pub type JsonDict = BTreeMap<String, Json>;

/// Empty rustage struct
pub fn default_rusage() -> rusage {
    rusage {
        ru_utime: timeval {
            tv_sec: 0,
            tv_usec: 0,
        },
        ru_stime: timeval {
            tv_sec: 0,
            tv_usec: 0,
        },
        ru_maxrss: 0,
        ru_ixrss: 0,
        ru_idrss: 0,
        ru_isrss: 0,
        ru_minflt: 0,
        ru_majflt: 0,
        ru_nswap: 0,
        ru_inblock: 0,
        ru_oublock: 0,
        ru_msgsnd: 0,
        ru_msgrcv: 0,
        ru_nsignals: 0,
        ru_nvcsw: 0,
        ru_nivcsw: 0,
    }
}

/// Rusage for self
pub fn get_rusage_self() -> Json {
    let mut ru: rusage = default_rusage();
    unsafe {
        getrusage(0, &mut ru);
    }
    rusage_to_json(&ru)
}


/// Rusage for childen
pub fn get_rusage_child() -> Json {
    let mut ru: rusage = default_rusage();
    unsafe {
        getrusage(1, &mut ru);
    }
    rusage_to_json(&ru)
}


/// Convert Rusage to json dict
pub fn rusage_to_json(rusage: &rusage) -> Json {
    let mut d = BTreeMap::new();
    {
        let mut ins = |n: &str, v| d.insert(n.to_string(), Json::I64(v as i64));
        ins("ru_utime", rusage.ru_utime.tv_sec as i64);
        ins("ru_stime", rusage.ru_stime.tv_sec as i64);
        ins("ru_maxrss", rusage.ru_maxrss);
        ins("ru_ixrss", rusage.ru_ixrss);
        ins("ru_idrss", rusage.ru_idrss);
        ins("ru_isrss", rusage.ru_isrss);
        ins("ru_minflt", rusage.ru_minflt);
        ins("ru_majflt", rusage.ru_majflt);
        ins("ru_nswap", rusage.ru_nswap);
        ins("ru_inblock", rusage.ru_inblock);
        ins("ru_oublock", rusage.ru_oublock);
        ins("ru_msgsnd", rusage.ru_msgsnd);
        ins("ru_msgrcv", rusage.ru_msgrcv);
        ins("ru_nsignals", rusage.ru_nsignals);
        ins("ru_nvcsw", rusage.ru_nvcsw);
        ins("ru_nivcsw", rusage.ru_nivcsw);
    }
    Json::Object(d)
}

const METADATA_PREFIX: &'static str = "_";


/// Tracking the metadata for one Martian chunk invocation
#[derive(Debug, Clone)]
pub struct Metadata {
    stage_name: String,
    stage_type: String,
    metadata_path: String,
    files_path: String,
    run_file: String,
    jobinfo: JsonDict,
    cache: HashSet<String>,
    start_time: DateTime<Local>,
    monitor_memory: bool,
}

pub fn make_timestamp(datetime: DateTime<Local>) -> String {
    datetime.format("%Y-%m-%d %H:%M:%S").to_string()
}

pub fn make_timestamp_now() -> String {
    return make_timestamp(Local::now());
}

impl Metadata {
    pub fn new(args: Vec<String>) -> Metadata {

        // # Take options from command line.
        // shell_cmd, stagecode_path, metadata_path, files_path, run_file = argv
        let md = Metadata {
            stage_name: args[0].clone(),
            stage_type: args[1].clone(),
            metadata_path: args[2].clone(),
            files_path: args[3].clone(),
            run_file: args[4].clone(),
            cache: HashSet::new(),
            start_time: Local::now(),
            jobinfo: BTreeMap::new(),
            monitor_memory: false,
        };
        
        md
    }

    /// Path within chunk
    pub fn make_path(&self, name: &str) -> PathBuf {
        let mut pb = PathBuf::from(self.metadata_path.clone());
        pb.push(METADATA_PREFIX.to_string() + name);
        pb
    }

    /// Write to a file inside the chunk
    pub fn write_raw(&mut self, name: &str, text: String) {
        let f = File::create(self.make_path(name));
        match f {
            Ok(mut ff) => {
                ff.write(text.as_bytes()).expect("io error");
                self.update_journal(name);
            },
            Err(e) => println!("err: {:?}", e)
        }
    }

    /// Update the Martian journal -- so that Martian knows what we've updated
    fn update_journal_main(&mut self, name: &str, force: bool) {
        let journal_name = if self.stage_type != "main" {
            format!("{}_{}", self.stage_type, name)
        } else {
            name.to_string()
        };

        if !self.cache.contains(name) || force {
            let run_file = format!("{}.{}", self.run_file, journal_name);
            let tmp_run_file = run_file.clone() + ".tmp";

            {
                let mut f = File::create(&tmp_run_file).expect("couldn't open file");
                f.write(make_timestamp_now().as_bytes()).expect("io error");
            };
            rename(&tmp_run_file, &run_file).expect("couldn't move file");
            self.cache.insert(journal_name);

        }
    }

    fn update_journal(&mut self, name: &str) {
        self.update_journal_main(name, false)
    }

    /*
    fn write_json(&mut self, name: &str, object: &Json) {
        // Serialize using `json::encode`
        let encoded = json::encode(object).unwrap();
        self.write_raw(name, encoded);
    }
    */

    /// Write JSON to a chunk file
    fn write_json_obj(&mut self, name: &str, object: &JsonDict) {
        // Serialize using `json::encode`
        let obj = &Json::Object(object.clone());
        let encoded = json::as_pretty_json(&obj);
        self.write_raw(name, format!("{}", encoded));
    }

    /// Read JSON from a chunk file
    fn read_json(&self, name: &str) -> Result<Json, ParserError> {
        let mut f = try!(File::open(self.make_path(name)));
        let mut buf = String::new();
        try!(f.read_to_string(&mut buf));

        Json::from_str(&buf)
    }

    fn read_json_obj(&self, name: &str) -> JsonDict {
        let r = self.read_json(name).expect("bad json");
        r.as_object().unwrap().clone()
    }

    fn read_json_obj_array(&self, name: &str) -> Vec<JsonDict> {
        let json = self.read_json(name).unwrap();
        let arr = json.as_array().unwrap();
        let r : Vec<JsonDict> = arr.into_iter().map(|o| o.as_object().unwrap().clone()).collect();
        r
    }


    fn _append(&mut self, name: &str, message: &str) {
        let filename = self.make_path(name);
        let mut file = OpenOptions::new().create(true).append(true).open(filename).expect("couldn't open");
        file.write(message.as_bytes()).expect("io error");
        file.write("\n".as_bytes()).expect("write");
        self.update_journal(name);
    }

    /// Write to _log
    pub fn log(&mut self, level: &str, message: &str) {
        self._append("log", &format!("{} [{}] {}", make_timestamp_now(), level, message))
    }

    pub fn log_time(&mut self, message: &str) {
        self.log("time", message)
    }

    pub fn alarm(&mut self, message: &str) {
        self._append("alarm", &format!("{} {}", make_timestamp_now(), message))
    }

    pub fn assert(&mut self, message: &str) {
        self._append("assert", &format!("{} {}", make_timestamp_now(), message))
    }

    /// Write finalized _jobinfo data
    pub fn update_jobinfo(&mut self) {
        let mut jobinfo = self.read_json_obj("jobinfo");

        jobinfo.insert("cwd".to_string(), Json::String(self.files_path.clone()));
        // jobinfo.insert("host", socket.gethostname());
        jobinfo.insert("pid".to_string(), Json::I64(unsafe { getpid() } as i64));
        let exe = env::current_exe().expect("current_exe").to_str().expect("exe").to_string();
        jobinfo.insert("rust_exe".to_string(), Json::String(exe));
        // jobinfo.insert("rust_version", sys.version);

        let get_env = |k| {
            let v = env::var(k);
            match v {
                Ok(s) => Json::String(s),
                Err(_) => Json::Null,
            }
        };

        match env::var("SGE_ARCH") {
            Ok(_) => {
                let mut d = BTreeMap::new();
                {
                    let mut ins = |n: &str, v| {
                        d.insert(n.to_string(), v);
                    };
                    ins("root", get_env("SGE_ROOT"));
                    ins("cell", get_env("SGE_CELL"));
                    ins("queue", get_env("QUEUE"));
                    ins("jobid", get_env("JOB_ID"));
                    ins("jobname", get_env("JOB_NAME"));
                    ins("sub_host", get_env("SGE_O_HOST"));
                    ins("sub_user", get_env("SGE_O_LOGNAME"));
                    ins("exec_host", get_env("HOSTNAME"));
                    ins("exec_user", get_env("LOGNAME"));
                }
                jobinfo.insert("sge".to_string(), Json::Object(d));
            }
            Err(_) => (),
        }

        self.write_json_obj("jobinfo", &jobinfo);
        self.jobinfo = jobinfo;
    }

    /// Check whether this job has requested memory monitoring (aka kneecapping) 
    pub fn monitor_memory(&self) -> bool {
        let jobinfo = self.read_json_obj("jobinfo");
        jobinfo.get("monitor_flag") == Some(&Json::String("monitor".to_string()))
    }

    /// Completed successfully
    pub fn complete(&mut self) {
        self.write_raw("complete", make_timestamp_now());
        self.shutdown();
    }

    /// Shutdown this chunk -- update jobinfo with final stats
    pub fn shutdown(&mut self) {
        self.log_time("__end__");

        // Common to fail() and complete()
        let endtime = Local::now();

        let mut jobinfo = self.read_json_obj("jobinfo");

        let mut wall_clock = BTreeMap::new();
        wall_clock.insert("start".to_string(),
                          Json::String(make_timestamp(self.start_time)));
        wall_clock.insert("end".to_string(), Json::String(make_timestamp(endtime)));
        wall_clock.insert("duration_seconds".to_string(),
                          Json::I64((endtime - self.start_time).num_seconds()));
        jobinfo.insert("wallclock".to_string(), Json::Object(wall_clock));

        let mut rusage = BTreeMap::new();
        rusage.insert("self".to_string(), get_rusage_self());
        rusage.insert("children".to_string(), get_rusage_child());
        jobinfo.insert("rusage".to_string(), Json::Object(rusage));

        self.write_json_obj("jobinfo", &jobinfo);

        // sys.exit does not actually exit the process but only exits the thread.
        // If this thread is not the main thread, use os._exit. This won't call
        // cleanup handlers, flush stdio buffers, etc. But calling done() from
        // another thread means the process exited with an error so this is okay.
        // if isinstance(threading.current_thread(), threading._MainThread)   {
        //    sys.exit(0);
        // }
        // else {
        //    os._exit(0);
        // }
    }
}

/// Shortcut function to decode a JSON `&str` into an object
pub fn obj_decode<T: Decodable>(s: JsonDict) -> T {
    let mut decoder = json::Decoder::new(Json::Object(s));
    Decodable::decode(&mut decoder).unwrap()
}

/// Shortcut function to decode a JSON `&str` into an object
pub fn json_decode<T: Decodable>(s: Json) -> T {
    let mut decoder = json::Decoder::new(s);
    Decodable::decode(&mut decoder).unwrap()
}

/// Shortcut function to decode a JSON `&str` into an object
pub fn obj_encode<T: ToJson>(v: &T) -> Json {
    v.to_json()
}

pub fn encode_to_json<T: Encodable>(v: &T) -> Json {

    let mut buf = String::new();
    {
        let mut encoder = json::Encoder::new(&mut buf);
        Encodable::encode(v, &mut encoder);
    }

    Json::from_str(&buf).unwrap()
}


pub struct Resource {
    pub __mem_gb: Option<f64>,
    pub __threads: Option<usize>,
}

impl Resource {
    pub fn none() -> Resource {
        Resource { 
            __mem_gb: None,
            __threads: None,
        }
    }
}


pub trait MartianStage {
    fn split(&self, args: JsonDict) -> JsonDict;
    fn main(&self, args: JsonDict, outs: JsonDict) -> JsonDict;
    fn join(&self, args: JsonDict, outs: JsonDict, chunk_defs: Vec<JsonDict>, chunk_outs: Vec<JsonDict>) -> JsonDict;
}


pub fn initialize(args: Vec<String>) -> Metadata {
    let mut md = Metadata::new(args);
    println!("got metadata: {:?}", md);
    md.update_jobinfo();

    md.log_time("__start__");

    md.update_journal("stdout");
    md.update_journal("stderr");

    // # Increase the maximum open file descriptors to the hard limit
    //
    // let _, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    // try:
    // resource.setrlimit(resource.RLIMIT_NOFILE, (hard, hard))
    // except Exception as e:
    // # Since we are still initializing, do not allow an unhandled exception.
    // # If the limit is not high enough, a preflight will catch it.
    // metadata.log("adapter", "Adapter could not increase file handle ulimit to %s: %s" % (str(hard), str(e)))
    // pass
    //

    // # Cache invocation and version JSON.
    // invocation = jobinfo["invocation"]
    // version = jobinfo["version"]

    md
}

pub fn do_split(stage: &MartianStage, mut md: Metadata)
{
    let args = md.read_json_obj("args");
    let stage_defs = stage.split(args);
    md.write_json_obj("stage_defs", &stage_defs);
    md.complete();
}

pub fn do_main(stage: &MartianStage, mut md: Metadata)
{
    let args = md.read_json_obj("args");
    let outs = md.read_json_obj("outs");

    let outs = stage.main(args, outs);

    md.write_json_obj("outs", &outs);
    md.complete();
}


pub fn do_join(stage: &MartianStage, mut md: Metadata)
{
    let args = md.read_json_obj("args");
    let outs = md.read_json_obj("outs");
    let chunk_defs = md.read_json_obj_array("chunk_defs");
    let chunk_outs = md.read_json_obj_array("chunk_outs");

    let outs = stage.join(args, outs, chunk_defs, chunk_outs);

    md.write_json_obj("outs", &outs);
    md.complete();
}

/// Log a panic to the martian output machinery
pub fn log_panic(md: &mut Metadata, panic: &panic::PanicInfo) {

    let payload =
        match panic.payload().downcast_ref::<String>() {
            Some(as_string) => format!("{}", as_string),
            None => format!("{:?}", panic.payload())
        };

    let loc = panic.location().expect("location");
    let msg = format!("{}: {}\n{}", loc.file(), loc.line(), payload);
    md.write_raw("errors", msg);
}

pub fn setup_logging(md: &Metadata)
{
    let level = log::LogLevelFilter::Debug;
    let log_path = md.make_path("log");

    let logger_config = fern::DispatchConfig {
        format: Box::new(|msg: &str, level: &log::LogLevel, _location: &log::LogLocation| {
            let time_str = Local::now().format("%Y-%m-%d %H:%M:%S").to_string();
            format!("[{}][{}] {}", time_str, level, msg)
        }),
        output: vec![fern::OutputConfig::stdout(), fern::OutputConfig::file(&log_path)],
        level: level,
    };

    if let Err(e) = fern::init_global_logger(logger_config, level) {
        panic!("Failed to initialize global logger: {}", e);
    }
}



pub fn martian_main(args: Vec<String>, stage_map: HashMap<String, Box<MartianStage>>) {

    info!("got args: {:?}", args);

    // setup Martian metadata
    let md = initialize(args);

    let mem_limit_gb = {
        let args = md.read_json_obj("args");
        match args.get("__mem_gb") {
            Some(j) => j.as_f64().unwrap_or(6.0),
            None => 6.0,
        }
    };

    // Hook rust logging up to Martian _log file
    setup_logging(&md);

    // Get the stage implementation
    let stage = stage_map.get(&md.stage_name).expect("couldn't find requested stage");

    // Setup monitor thread -- this handles heartbeat & memory checking
    let stage_done = Arc::new(AtomicBool::new(false));
    let monitor_memory = md.monitor_memory();
    let mut md_monitor = md.clone();
    let stage_done_monitor = stage_done.clone();
    let monitor_handle = thread::spawn(move || {
        loop {

            // Write to the heartbeat file: so Martian knows we're still alive
            md_monitor.update_journal_main("heartbeat", true);
            let one_min = time::Duration::from_millis(60000);
            thread::park_timeout(one_min);

            if monitor_memory {
                // Check the total memory consumption of the stage. Kill ourself if we go over the limit
                let mut ru_self: rusage = default_rusage();
                unsafe { getrusage(0, &mut ru_self); }

                let mut ru_child: rusage = default_rusage();
                unsafe {  getrusage(1, &mut ru_child); }

                // maxrss is reported in kb
                let max_rss = max(ru_self.ru_maxrss, ru_child.ru_maxrss) * 1024;
                if max_rss > (mem_limit_gb * 1e9) as i64 {
                    // Shutdown the process due to memory
                    info!("Calling panic due to mem consumption");
                    panic!("Memory consumption exceed limit. Maxrss: {}, Limit: {}", max_rss, (mem_limit_gb * 1e9) as usize);
                }
            }

            if stage_done_monitor.load(Ordering::Relaxed) {
                break;
            }
        }
    });

    // Setup panic hook. If a stage panics, we'll shutdown cleanly to martian
    let p = panic::take_hook();
    let mut _panic_md = md.clone();
    panic::set_hook(Box::new(move |info| {
        let mut panic_md = _panic_md.clone();
        let backtrace = Backtrace::new();

        let thread = thread::current();
        let thread = thread.name().unwrap_or("unnamed");

        let msg = match info.payload().downcast_ref::<&'static str>() {
            Some(s) => *s,
            None => match info.payload().downcast_ref::<String>() {
                Some(s) => &**s,
                None => "Box<Any>",
            }
        };

        let msg =
            match info.location() {
                Some(location) => {
                    format!("thread '{}' panicked at '{}': {}:{}{:?}",
                           thread,
                           msg,
                           location.file(),
                           location.line(),
                           backtrace)
                }
                None => format!("thread '{}' panicked at '{}'{:?}", thread, msg, backtrace),
            };

        error!("{}", msg);
        panic_md.write_raw("errors", msg);
        p(info);
    }));


    if md.stage_type == "split"
    {
        do_split(stage.as_ref(), md);
    }
    else if md.stage_type == "main"
    {
        do_main(stage.as_ref(), md);
    }
    else if md.stage_type == "join"
    {
        do_join(stage.as_ref(), md);
    }
    else
    {
        panic!("Unrecognized stage type");
    };

    stage_done.store(true, Ordering::Relaxed);
    monitor_handle.thread().unpark();
    monitor_handle.join().unwrap();
}


/// Helper function to bump up file handle limit if a stage requires it.
pub fn set_file_handle_limit(request: usize) -> isize {
    let mut req = rlimit {
        rlim_cur: request as c_ulong,
        rlim_max: request as c_ulong,
    };

    unsafe {
        libc::getrlimit(libc::RLIMIT_NOFILE, &mut req);
        req.rlim_cur = min(req.rlim_max, request as u64);
        libc::setrlimit(libc::RLIMIT_NOFILE, &mut req) as isize
    }
}
