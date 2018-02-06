# Profiling Rust code

## Profiling a full executable
1. Enable debug symbols in the release build by adding the following lines to the Cargo.toml
```rust
[profile.release]
debug = true
```
2. Compile and generate the executable
3. Run the executable through perf
```bash
perf record -g [executable]
```
4. Clone the flamegraph repo and add it to yout PATH (<https://github.com/brendangregg/FlameGraph>)
5. Build the flamegraph
```bash
perf script | stackcollapse-perf.pl | flamegraph.pl > flame.svg
```
6. Flame.svg will contain a handy visualization. (I have noticed that some symbols are not read consistently leading to lots of [unknown] blocks in the flamegraph.)
 
## Profiling a benchmark [Nightly feature]

Often you might want to profile a specific function/set of functions instead of the whole executable. A good way top do this is rust is using benchmark tests (https://doc.rust-lang.org/1.12.1/book/benchmark-tests.html)
1. Write a Rust benchmark test . Call it my_bench_test()
2. Enable debug symbols in the release build by adding the following lines to the Cargo.toml
```rust
[profile.release]
debug = true
```
3. Build the benchmark executable. It would be in target/release/
```bash
cargo bench --no-run
```
4. Run the executable through perf
```bash
perf record -g target/release/[executable] --bench my_bench_test
```
Here I have used my_bench_test to filter the benchmark I am interested from all the benchmark tests.
5. Build the flamegraph
```bash
perf script | stackcollapse-perf.pl | flamegraph.pl > flame.svg
```

If you want to stick to the Rust stable release, use the bluss/bencher(https://github.com/bluss/bencher) for writing benchmark tests

### Reference
http://blog.adamperry.me/rust/2016/07/24/profiling-rust-perf-flamegraph/
