// Export versions of __builtin_cpu_supports for external calls.
// This is a tiny library used for testing in preflights.

int cpu_supports_sse4_2() {
  return __builtin_cpu_supports("sse4.2");
}

int cpu_supports_popcnt() {
  return __builtin_cpu_supports("popcnt");
}

int cpu_supports_avx() {
  return __builtin_cpu_supports("avx");
}
