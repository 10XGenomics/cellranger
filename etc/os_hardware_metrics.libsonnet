// OS and hardware telemetry metrics shared across products.
{
  local os_group = ['os_stats'],
  local hw_group = ['hardware_stats'],
  local with_values(name, groups, values, extra={}) = { name: name, groups: groups, values: values } + extra,
  local semver_metric(name, special, thresholds) = with_values(name, os_group, { semver_thresholds: thresholds }, { special: special }),

  // Build the OS/hardware block; product-specific command groups are injected.
  os_hardware_metrics(command_groups)::
    [
      with_values('event_count', os_group + hw_group + command_groups, { exponential: { max: 768, min: 3, multiple: 2 } }, { recent_count: 'os_stats' }),
    ] +
    [
      semver_metric(item.name, item.special, item.thresholds)
      for item in [
        { name: 'glibc', special: 'glibc', thresholds: ['2.12', '2.17', '2.28', '2.34'] },
        { name: 'kernel', special: 'kernel', thresholds: ['3.10', '4.4', '4.15', '4.18', '5.4', '5.14', '5.15'] },
      ]
    ] +
    [
      { groups: os_group, name: 'distro', special: 'distro' },
      { groups: os_group, name: 'container', special: 'container' },
      { groups: os_group, name: 'nextflow_task_workdir_set', env: 'NXF_TASK_WORKDIR', values: { boolean: true } },
      { groups: os_group, name: 'nextflow_assets_set', env: 'NXF_ASSETS', values: { boolean: true } },
      with_values('cgroup_memory_soft_limit', os_group, { exponential: { min: 1073741824, multiple: 2 } }, { file_content: { cgroup: 'memory', line: '[0-9]+', path: ['memory.soft_limit_in_bytes', 'memory.memsw.soft_limit_in_bytes', 'memory.high'] } }),
      with_values('cgroup_memory_hard_limit', os_group, { exponential: { min: 1073741824, multiple: 2 } }, { file_content: { cgroup: 'memory', line: '[0-9]+', path: ['memory.limit_in_bytes', 'memory.memsw.limit_in_bytes', 'memory.max'] } }),
      { filesystem: { property: 'type', which: 'bins' }, groups: os_group, name: 'bin_filesystem', preset: 'filesystem' },
      { filesystem: { property: 'type', which: 'pipestance' }, groups: command_groups, name: 'pipestance_filesystem', preset: 'filesystem' },
    ] +
    [
      with_values(item.name, os_group, { thresholds: [1024, 2048, 65536, 131072] }, { rlimit: { hard: false, resource: item.resource } })
      for item in [
        { name: 'process_ulimit', resource: 'nproc' },
        { name: 'vmem_ulimit', resource: 'as' },
      ]
    ] +
    [
      with_values('transparent_huge_pages', os_group, { match: [{ exp: 'always' }, { exp: 'madvise' }, { exp: 'never' }] }, { file_content: { line: '\\[(\\w+)\\]', path: ['/sys/kernel/mm/transparent_hugepage/enabled'] } }),
      with_values('cpu_mfg', hw_group, { match: [{ exp: '(?i:Intel)' }, { exp: '(?i:AMD)' }] }, { file_content: { line: 'vendor_id\\s*:\\s*(\\S*)', path: ['/proc/cpuinfo'] } }),
      with_values('avx', hw_group, { match: [{ exp: '\\bavx512' }, { exp: '\\bavx2\\b' }, { exp: '\\bavx\\b' }] }, { file_content: { line: 'flags\\s*:\\s*(.*)', path: ['/proc/cpuinfo'] } }),
      with_values('processor_count', hw_group + command_groups, { thresholds: [7, 15, 23, 31, 47, 65, 95, 127] }, { file_content: { last_line: 'processor\\s*:\\s*(\\d+)', path: ['/proc/cpuinfo'] } }),
      with_values('socket_count', hw_group, { thresholds: [0, 1, 3] }, { file_content: { last_line: 'physical id\\s*:\\s*(\\d+)', path: ['/proc/cpuinfo'] } }),
      with_values('cores_per_socket', hw_group, { thresholds: [4, 8, 12, 16, 24, 32] }, { file_content: { line: 'cpu cores\\s*:\\s*(\\d+)', path: ['/proc/cpuinfo'] } }),
      with_values('memory', hw_group + command_groups, { exponential: { min: 8388608, multiple: 2 } }, { file_content: { line: '^MemTotal:\\s*([0-9]+)', path: ['/proc/meminfo'] } }),
    ],
}
