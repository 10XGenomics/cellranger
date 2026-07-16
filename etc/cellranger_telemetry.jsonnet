local os_hw = import 'os_hardware_metrics.libsonnet';
local shared = import 'shared.libsonnet';

local count_multi_vdj_groups = [
  'count_success',
  'multi_success',
  'vdj_success',
];
local failure_groups = [
  'aggr_failure',
  'annotate_failure',
  'count_failure',
  'mkfastq_failure',
  'mkref_failure',
  'mkvdjref_failure',
  'multi_failure',
  'reanalyze_failure',
  'vdj_failure',
];
local sample_metrics_groups = [
  'sample_metrics',
];

local named_re_matcher = shared.named_re_matcher;
local path_metric_float = shared.path_metric_float;
local lib_metrics_float(name, groups=count_multi_vdj_groups) = {
  local metric = path_metric_float(name, '${lib_metrics}', groups),
  file_content: metric.file_content,
  groups: metric.groups,
  name: metric.name,
  values: metric.values,
};
local lib_metrics_preset(name, json_path, preset, groups=count_multi_vdj_groups) = {
  file_content: {
    json_path: [json_path],
    path: ['${lib_metrics}'],
  },
  groups: groups,
  name: name,
  preset: preset,
};
local sample_tag_codes =
  [
    std.format('%s-%s%02d', [plate, row, col])
    for plate in ['A', 'B', 'C', 'D']
    for row in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    for col in std.range(1, 12)
  ] +
  [
    std.format('BC%03d', [idx])
    for idx in std.range(1, 16)
  ] +
  [
    std.format('AB%03d', [idx])
    for idx in std.range(1, 16)
  ] +
  [
    std.format('OB%d', [idx])
    for idx in std.range(1, 4)
  ] +
  [
    std.format('CMO%03d', [idx])
    for idx in std.range(301, 312)
  ];
local sample_tag_matches = [
  {
    exp: '^$',
    key: 'None',
  },
  {
    exp: '\\|',
    key: 'multiple_tags',
  },
] + [
  {
    exp: code,
  }
  for code in sample_tag_codes
];
local sample_metrics_float(name) = path_metric_float(name, '${sample_metrics}', sample_metrics_groups);

local cell_annotation_families = [
  { model: 'Azimuth', key: 'azimuth_v1_top_fine_cell_types', preset: 'azimuth_fine_cell_types', species: 'human' },
  { model: 'Azimuth', key: 'azimuth_v1_top_coarse_cell_types', preset: 'azimuth_coarse_cell_types', species: 'human' },
  { model: '10x_Cloud', key: 'human_pca_v1_beta_top_fine_cell_types', preset: 'pca_v1_beta_fine_cell_types', species: 'human' },
  { model: '10x_Cloud', key: 'human_pca_v1_beta_top_coarse_cell_types', preset: 'pca_v1_beta_coarse_cell_types', species: 'human' },
  { model: '10x_Cloud', key: 'mouse_pca_v1_beta_top_fine_cell_types', preset: 'pca_v1_beta_fine_cell_types', species: 'mouse' },
  { model: '10x_Cloud', key: 'mouse_pca_v1_beta_top_coarse_cell_types', preset: 'pca_v1_beta_coarse_cell_types', species: 'mouse' },
];
local cell_annotation_metrics(path_var, groups) = std.flattenArrays([
  [
    {
      file_content: {
        json_path: ['cell_annotation_analysis_metrics', f.model, f.key, std.toString(rank)],
        path: [path_var],
      },
      groups: groups,
      name: f.key + '.' + f.species + '.' + std.toString(rank),
      preset: f.preset,
    }
    for rank in std.range(1, 10)
  ]
  for f in cell_annotation_families
]);
{
  config_version: 1,
  product: 'cellranger',
  version: '10.1.0',

  groups: shared.groups {
    aggr_success: {
      globs: {
        aggr_metrics: 'SC_RNA_AGGREGATOR_CS/COUNT_AGGR/SC_RNA_AGGREGATOR/SUMMARIZE_AGGREGATED_REPORTS/fork0/files/summary.json',
      },
      subcommands: [
        'aggr',
      ],
      when: 'success',
    },
    annotate_failure: {
      subcommands: [
        'annotate',
      ],
      when: 'failure',
    },
    annotate_success: {
      globs: {
        lib_metrics: 'CELLRANGER_ANNOTATE_CS*/MERGE_ALL_CELL_TYPE_METRICS/fork0/files/merged_cell_annotation_metrics.json',
      },
      subcommands: [
        'annotate',
      ],
      when: 'success',
    },
    count_success: {
      globs: {
        lib_metrics: 'SC_*/SC_MULTI_CORE/MERGE_METRICS/fork0/files/summary.json',
      },
      subcommands: [
        'count',
      ],
      when: 'success',
    },
    library_metrics_vdj: {
      globs: {
        vdj_lib_metrics: 'SC_*/SC_MULTI_CORE/MULTI_GEM_WELL_PROCESSOR/VDJ_*_GEM_WELL_PROCESSOR/LIB_VDJ_ANALYZER/MERGE_METRICS/fork0/files/summary.json',
      },
      subcommands: [
        'multi',
        'vdj',
      ],
      when: 'success',
    },
    mkvdjref_failure: {
      subcommands: [
        'mkvdjref',
      ],
      when: 'failure',
    },
    mkvdjref_success: {
      subcommands: [
        'mkvdjref',
      ],
      when: 'success',
    },
    multi_failure: {
      subcommands: [
        'multi',
      ],
      when: 'failure',
    },
    multi_success: {
      globs: {
        lib_metrics: 'SC_*/SC_MULTI_CORE/MERGE_METRICS/fork0/files/summary.json',
        vdj_lib_metrics: 'SC_*/SC_MULTI_CORE/MULTI_GEM_WELL_PROCESSOR/VDJ_*_GEM_WELL_PROCESSOR/LIB_VDJ_ANALYZER/MERGE_METRICS/fork0/files/summary.json',
      },
      subcommands: [
        'multi',
      ],
      when: 'success',
    },
    sample_metrics: {
      globs: {
        lib_metrics: 'SC_*/SC_MULTI_CORE/MERGE_METRICS/fork0/files/summary.json',
        sample_metrics: 'SC_*/SC_MULTI_CORE/SANITIZE_MAP_CALLS/fork0/files/*_summary.json',
      },
      subcommands: [
        'multi',
      ],
      when: 'success',
    },
    sample_metrics_vdj: {
      globs: {
        vdj_sample_metrics: 'SC_*/SC_MULTI_CORE/MULTI_GEM_WELL_PROCESSOR/VDJ_*_GEM_WELL_PROCESSOR/PER_SAMPLE_VDJ_ANALYZER/MERGE_METRICS/fork*/files/summary.json',
      },
      subcommands: [
        'multi',
      ],
      when: 'success',
    },
    vdj_failure: {
      subcommands: [
        'vdj',
      ],
      when: 'failure',
    },
    vdj_success: {
      globs: {
        lib_metrics: 'SC_*/SC_MULTI_CORE/MERGE_METRICS/fork0/files/summary.json',
        vdj_lib_metrics: 'SC_*/SC_MULTI_CORE/MULTI_GEM_WELL_PROCESSOR/VDJ_*_GEM_WELL_PROCESSOR/LIB_VDJ_ANALYZER/MERGE_METRICS/fork0/files/summary.json',
      },
      subcommands: [
        'vdj',
      ],
      when: 'success',
    },
  },

  metrics:
    os_hw.os_hardware_metrics([
      'count_success',
      'multi_success',
      'vdj_success',
      'aggr_success',
      'annotate_success',
      'reanalyze_success',
      'mkvdjref_success',
      'mkfastq_success',
      'mkref_success',
      'count_failure',
      'multi_failure',
      'vdj_failure',
      'aggr_failure',
      'annotate_failure',
      'reanalyze_failure',
      'mkvdjref_failure',
      'mkfastq_failure',
      'mkref_failure',
    ]) +
    [
      {
        groups: [
          'count_failure',
          'multi_failure',
          'vdj_failure',
          'aggr_failure',
          'annotate_failure',
          'reanalyze_failure',
          'mkvdjref_failure',
          'mkfastq_failure',
          'mkref_failure',
        ],
        name: 'failed_stage',
        special: 'failed_stage',
      },
      {
        file_content: {
          last_line: std.join('|', [x.exp for x in shared.error_matches]),
          path: ['_log'],
        },
        groups: failure_groups,
        name: 'log_error',
        values: { match: shared.error_matches },
      },
      {
        file_content: {
          last_line: std.join('|', [x.exp for x in shared.error_matches]),
          path: ['_errors'],
          stage_relative: true,
        },
        groups: failure_groups,
        name: 'error',
        values: { match: shared.error_matches },
      },
      {
        groups: [
          'count_success',
          'multi_success',
          'vdj_success',
          'aggr_success',
          'annotate_success',
          'reanalyze_success',
          'mkvdjref_success',
          'mkfastq_success',
          'mkref_success',
          'count_failure',
          'multi_failure',
          'vdj_failure',
          'aggr_failure',
          'annotate_failure',
          'reanalyze_failure',
          'mkvdjref_failure',
          'mkfastq_failure',
          'mkref_failure',
        ],
        name: 'subcommand',
        special: 'subcommand',
      },
      {
        flag: '--expect-cells',
        groups: [
          'count_success',
          'count_failure',
        ],
        name: 'expect_cells',
        values: {
          float: true,
        },
      },
      {
        flag: '--force-cells',
        groups: [
          'count_success',
          'count_failure',
        ],
        name: 'force_cells',
        values: {
          float: true,
        },
      },
      {
        flag: '--create-bam',
        groups: [
          'count_success',
          'count_failure',
        ],
        name: 'create_bam',
        values: {
          boolean: true,
        },
      },
      {
        flag: '--no-bam',
        groups: [
          'count_success',
          'count_failure',
        ],
        name: 'no_bam',
        values: {
          boolean: true,
        },
      },
      {
        flag: '--nosecondary',
        groups: [
          'count_success',
          'count_failure',
        ],
        name: 'nosecondary',
        values: {
          boolean: true,
        },
      },
      {
        flag: '--r1-length',
        groups: [
          'count_success',
          'count_failure',
        ],
        name: 'r1_length',
        values: {
          float: true,
        },
      },
      {
        flag: '--r2-length',
        groups: [
          'count_success',
          'count_failure',
        ],
        name: 'r2_length',
        values: {
          float: true,
        },
      },
      {
        flag: '--include-introns',
        groups: [
          'count_success',
          'count_failure',
        ],
        name: 'include_introns',
        values: {
          boolean: true,
        },
      },
      {
        flag: '--chemistry',
        groups: [
          'count_success',
          'count_failure',
        ],
        name: 'chemistry',
        preset: 'chemistry',
      },
      {
        flag: '--no-libraries',
        groups: [
          'count_success',
          'count_failure',
        ],
        name: 'no_libraries',
        values: {
          boolean: true,
        },
      },
      {
        flag: '--check-library-compatibility',
        groups: [
          'count_success',
          'count_failure',
        ],
        name: 'check_library_compatibility',
        values: {
          boolean: true,
        },
      },
      {
        flag: '--filter-aggregates',
        groups: [
          'count_success',
          'count_failure',
        ],
        name: 'filter_aggregates',
        values: {
          boolean: true,
        },
      },
      {
        flag: '--min-crispr-umi',
        groups: [
          'count_success',
          'count_failure',
        ],
        name: 'min-crispr-umi',
        values: {
          float: true,
        },
      },
      {
        flag: '--denovo',
        groups: [
          'vdj_success',
          'vdj_failure',
        ],
        name: 'denovo',
        values: {
          boolean: true,
        },
      },
      {
        flag: '--skip-clonotyping',
        groups: [
          'vdj_success',
          'vdj_failure',
        ],
        name: 'skip_clonotyping',
        values: {
          boolean: true,
        },
      },
      {
        flag: '--chain',
        groups: [
          'vdj_success',
          'vdj_failure',
        ],
        name: 'chain',
        values: {
          match: [
            named_re_matcher('auto'),
            named_re_matcher('TR'),
            named_re_matcher('IG'),
            named_re_matcher('TR_GD'),
          ],
        },
      },
      {
        flag: '--memgb',
        groups: [
          'mkref_success',
          'mkref_failure',
          'mkvdjref_success',
          'mkvdjref_failure',
        ],
        name: 'memgb',
        values: {
          float: true,
        },
      },
      {
        flag: '--nthreads',
        groups: [
          'mkref_success',
          'mkref_failure',
          'mkvdjref_success',
          'mkvdjref_failure',
        ],
        name: 'nthreads',
        values: {
          float: true,
        },
      },
      {
        flag: '--normalize',
        groups: [
          'aggr_success',
          'aggr_failure',
        ],
        name: 'normalize',
        values: {
          match: [
            named_re_matcher('mapped'),
            named_re_matcher('none'),
          ],
        },
      },
      {
        flag: '--normalize',
        groups: [
          'aggr_success',
          'aggr_failure',
        ],
        name: 'nosecondary',
        values: {
          boolean: true,
        },
      },
      {
        flag: '--min-crispr-umi',
        groups: [
          'aggr_success',
          'aggr_failure',
        ],
        name: 'min-crispr-umi',
        values: {
          float: true,
        },
      },
      {
        flag: '--enable-tsne',
        groups: [
          'aggr_success',
          'aggr_failure',
        ],
        name: 'enable-tsne',
        values: {
          boolean: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'gene_expression',
            'emptydrops_minimum_umis',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'emptydrops_minimum_umis',
        values: {
          float: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'gene_expression',
            'global_minimum_umis',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'global_minimum_umis',
        values: {
          float: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'gene_expression',
            'r1_length',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'gene_expression_r1_length',
        values: {
          float: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'gene_expression',
            'r2_length',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'gene_expression_r2_length',
        values: {
          float: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'gene_expression',
            'chemistry',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'chemistry',
        preset: 'chemistry',
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'gene_expression',
            'expect_cells',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'expect_cells',
        values: {
          float: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'gene_expression',
            'force_cells',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'force_cells',
        values: {
          float: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: ['gene_expression', 'cell_calling_mode'],
        },
        groups: ['multi_success', 'multi_failure'],
        name: 'cell_calling_mode',
        values: {
          match: [
            named_re_matcher('gradient'),
            named_re_matcher('ordmag_nonambient'),
            named_re_matcher('ordmag'),
          ],
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'gene_expression',
            'nosecondary',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'nosecondary',
        values: {
          boolean: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'gene_expression',
            'include_introns',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'include_introns',
        values: {
          boolean: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'gene_expression',
            'check_library_compatibility',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'check_library_compatibility',
        values: {
          boolean: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'gene_expression',
            'aligner',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'aligner',
        values: {
          match: [
            named_re_matcher('hurtle'),
            named_re_matcher('star'),
          ],
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'gene_expression',
            'create_bam',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'create_bam',
        values: {
          boolean: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'gene_expression',
            'filter_probes',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'filter_probes',
        values: {
          boolean: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'gene_expression',
            'filter_high_occupancy_gems',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'filter_high_occupancy_gems',
        values: {
          boolean: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'gene_expression',
            'min_assignment_confidence',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'min_assignment_confidence',
        values: {
          float: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'gene_expression',
            'skip_cell_annotation',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'skip_cell_annotation',
        values: {
          boolean: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'gene_expression',
            'detect_chemistry_sample_reads',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'detect_chemistry_sample_reads',
        values: {
          float: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'gene_expression',
            'detect_chemistry_total_reads',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'detect_chemistry_total_reads',
        values: {
          float: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'gene_expression',
            'umi_min_read_length',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'umi_min_read_length',
        values: {
          float: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'gene_expression',
            'min_major_probe_bc_frac',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'min_major_probe_bc_frac',
        values: {
          float: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'feature',
            'r1_length',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'feature_r1_length',
        values: {
          float: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'feature',
            'r2_length',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'feature_r2_length',
        values: {
          float: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'feature',
            'filter_aggregates',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'filter_aggregates',
        values: {
          boolean: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'feature',
            'min_crispr_umi',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'min_crispr_umi',
        values: {
          float: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'feature',
            'crispr_probe_barcode_location',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'crispr_probe_barcode_location',
        values: {
          float: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'vdj',
            'r1_length',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'r1_length',
        values: {
          float: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'vdj',
            'r2_length',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'r2_length',
        values: {
          float: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'vdj',
            'multiplet_filter',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'multiplet_filter',
        values: {
          boolean: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'vdj',
            'shared_contig_filter',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'shared_contig_filter',
        values: {
          boolean: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'vdj',
            'umi_baseline_filter',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'umi_baseline_filter',
        values: {
          boolean: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'vdj',
            'min_contig_length',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'min_contig_length',
        values: {
          float: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'vdj',
            'skip_clonotyping',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'skip_clonotyping',
        values: {
          boolean: true,
        },
      },
      {
        file_content: {
          converter: 'multicsv',
          flag: '--csv',
          json_path: [
            'vdj',
            'max_reads_per_barcode',
          ],
        },
        groups: [
          'multi_success',
          'multi_failure',
        ],
        name: 'max_reads_per_barcode',
        values: {
          float: true,
        },
      },
    ] + shared.martian_metrics(
      [
        'count_success',
        'multi_success',
        'vdj_success',
        'aggr_success',
        'annotate_success',
        'reanalyze_success',
        'mkvdjref_success',
        'mkfastq_success',
        'mkref_success',
      ],
      [
        'count_failure',
        'multi_failure',
        'vdj_failure',
        'aggr_failure',
        'annotate_failure',
        'reanalyze_failure',
        'mkvdjref_failure',
        'mkfastq_failure',
        'mkref_failure',
      ],
      no_perf_success_groups=['mkref_success', 'mkvdjref_success'],
    ) + [
      {
        file_content: { json_path: ['samples_count'], path: ['${aggr_metrics}'] },
        groups: ['aggr_success'],
        name: 'samples_count',
        values: { float: true },
      },
      {
        file_content: { json_path: ['reference_genomes'], path: ['${lib_metrics}'] },
        groups: ['count_success', 'multi_success', 'vdj_success'],
        name: 'reference_genomes',
        preset: 'reference_genomes',
      },
      {
        file_content: { json_path: ['vdj_reference_genomes'], path: ['${vdj_lib_metrics}'] },
        groups: ['library_metrics_vdj'],
        name: 'vdj_reference_genomes',
        preset: 'reference_genomes',
      },
      {
        file_content: { json_path: ['target_panel_hash'], path: ['${lib_metrics}'] },
        groups: ['multi_success'],
        name: 'target_panel_hash',
        preset: 'target_panel_hash',
      },
      {
        file_content: { json_path: ['target_panel_name'], path: ['${lib_metrics}'] },
        groups: ['multi_success'],
        name: 'target_panel_name',
        preset: 'target_panel_name',
      },
      {
        file_content: { json_path: ['chemistry_name'], path: ['${lib_metrics}'] },
        groups: ['count_success', 'multi_success', 'vdj_success'],
        name: 'chemistry_name',
        preset: 'chemistry',
      },
    ] +
    [
      {
        file_content: { json_path: ['chemistry_defs', x[1], 'name'], path: ['${lib_metrics}'] },
        groups: ['count_success', 'multi_success', 'sample_metrics'],
        name: x[0] + '_chemistry',
        preset: 'chemistry',
      }
      for x in [
        ['ANTIBODY', 'Antibody Capture'],
        ['ANTIGEN', 'Antigen Capture'],
        ['CRISPR', 'CRISPR Guide Capture'],
        ['Custom', 'Custom'],
        ['GEX', 'Gene Expression'],
        ['MULTIPLEXING', 'Multiplexing Capture'],
        ['VDJ', 'VDJ'],
        ['VDJ_B', 'VDJ-B'],
        ['VDJ_T', 'VDJ-T'],
        ['VDJ_T_GD', 'VDJ-T-GD'],
      ]
    ] +
    [
      lib_metrics_float(name)
      for name in [
        'barcode_subseq_0_tot_count',
        'barcode_subseq_0_obs_count',
        'barcode_subseq_0_exp_count',
        'barcode_subseq_0_logp',
        'barcode_subseq_0_ge_0_count',
        'barcode_subseq_0_ge_1_count',
        'barcode_subseq_0_ge_2_count',
        'barcode_subseq_0_ge_3_count',
        'barcode_subseq_0_ge_4_count',
        'barcode_subseq_0_ge_5_count',
        'barcode_subseq_1_tot_count',
        'barcode_subseq_1_obs_count',
        'barcode_subseq_1_exp_count',
        'barcode_subseq_1_logp',
        'barcode_subseq_1_ge_0_count',
        'barcode_subseq_1_ge_1_count',
        'barcode_subseq_1_ge_2_count',
        'barcode_subseq_1_ge_3_count',
        'barcode_subseq_1_ge_4_count',
        'barcode_subseq_1_ge_5_count',
        'barcode_subseq_sum_logp',
      ]
    ] +
    [
      lib_metrics_float(name)
      for name in [
        'ANTIBODY_bc_bases_with_q30_frac',
        'ANTIBODY_bc_bases_with_q30_in_gel_bead_frac',
        'ANTIBODY_bc_bases_with_q30_in_probe_frac',
        'ANTIBODY_feature_reads_in_cells',
        'ANTIBODY_filtered_bcs_transcriptome_union',
        'ANTIBODY_frac_feature_reads_usable',
        'ANTIBODY_good_bc_frac',
        'ANTIBODY_good_bc_in_gel_bead_frac',
        'ANTIBODY_good_bc_in_probe_frac',
        'ANTIBODY_good_umi_frac',
        'ANTIBODY_multi_cdna_pcr_dupe_reads_frac',
        'ANTIBODY_read_bases_with_q30_frac',
        'ANTIBODY_reads_lost_to_aggregate_GEMs',
        'ANTIBODY_reads_per_cell',
        'ANTIBODY_recognized_feature_bc_frac',
        'ANTIBODY_total_read_pairs',
        'ANTIBODY_umi_bases_with_q30_frac',
        'ANTIBODY_unprocessed_read_pairs',
      ]
    ] +
    [
      lib_metrics_float(name)
      for name in [
        'ANTIGEN_bc_bases_with_q30_frac',
        'ANTIGEN_bc_bases_with_q30_in_gel_bead_frac',
        'ANTIGEN_bc_bases_with_q30_in_probe_frac',
        'ANTIGEN_feature_reads_in_cells',
        'ANTIGEN_good_bc_frac',
        'ANTIGEN_good_umi_frac',
        'ANTIGEN_multi_cdna_pcr_dupe_reads_frac',
        'ANTIGEN_read_bases_with_q30_frac',
        'ANTIGEN_reads_per_cell',
        'ANTIGEN_total_read_pairs',
        'ANTIGEN_umi_bases_with_q30_frac',
        'ANTIGEN_unprocessed_read_pairs',
      ]
    ] +
    [
      lib_metrics_float(name)
      for name in [
        'CRISPR_bc_bases_with_q30_frac',
        'CRISPR_bc_bases_with_q30_in_gel_bead_frac',
        'CRISPR_bc_bases_with_q30_in_probe_frac',
        'CRISPR_feature_bc_extracted_frac',
        'CRISPR_feature_reads_in_cells',
        'CRISPR_frac_feature_reads_usable',
        'CRISPR_good_bc_frac',
        'CRISPR_good_bc_in_gel_bead_frac',
        'CRISPR_good_bc_in_probe_frac',
        'CRISPR_good_umi_frac',
        'CRISPR_multi_cdna_pcr_dupe_reads_frac',
        'CRISPR_read_bases_with_q30_frac',
        'CRISPR_reads_per_cell',
        'CRISPR_recognized_feature_bc_frac',
        'CRISPR_total_read_pairs',
        'CRISPR_umi_bases_with_q30_frac',
        'CRISPR_unprocessed_read_pairs',
        'CRISPR_unrecognized_feature_bc_frac',
      ]
    ] +
    [
      lib_metrics_float(name)
      for name in [
        'Custom_bc_bases_with_q30_frac',
        'Custom_bc_bases_with_q30_in_gel_bead_frac',
        'Custom_bc_bases_with_q30_in_probe_frac',
        'Custom_feature_reads_in_cells',
        'Custom_frac_feature_reads_usable',
        'Custom_good_bc_frac',
        'Custom_good_bc_in_gel_bead_frac',
        'Custom_good_bc_in_probe_frac',
        'Custom_good_umi_frac',
        'Custom_multi_cdna_pcr_dupe_reads_frac',
        'Custom_read_bases_with_q30_frac',
        'Custom_reads_per_cell',
        'Custom_recognized_feature_bc_frac',
        'Custom_total_read_pairs',
        'Custom_umi_bases_with_q30_frac',
        'Custom_unprocessed_read_pairs',
      ]
    ] +
    [
      lib_metrics_float(name)
      for name in [
        'MULTIPLEXING_bc_bases_with_q30_frac',
        'MULTIPLEXING_bc_bases_with_q30_in_gel_bead_frac',
        'MULTIPLEXING_bc_bases_with_q30_in_probe_frac',
        'MULTIPLEXING_feature_reads_in_cells',
        'MULTIPLEXING_good_bc_frac',
        'MULTIPLEXING_good_umi_frac',
        'MULTIPLEXING_multi_cdna_pcr_dupe_reads_frac',
        'MULTIPLEXING_read_bases_with_q30_frac',
        'MULTIPLEXING_reads_per_cell',
        'MULTIPLEXING_total_read_pairs',
        'MULTIPLEXING_umi_bases_with_q30_frac',
        'MULTIPLEXING_unprocessed_read_pairs',
      ]
    ] +
    [
      lib_metrics_float(name)
      for name in [
        'bc_bases_with_q30_frac',
        'bc_bases_with_q30_in_gel_bead_frac',
        'bc_bases_with_q30_in_probe_frac',
        'estimated_gdna_content',
        'estimated_gdna_unspliced_threshold',
        'filtered_bcs_mean_counts',
        'filtered_bcs_median_counts',
        'filtered_bcs_transcriptome_union',
        'filtered_reads_per_filtered_bc',
        'good_bc_frac',
        'good_umi_frac',
        'multi_cdna_pcr_dupe_reads_frac',
        'multi_filtered_bcs_conf_mapped_barcoded_reads_cum_frac',
        'multi_transcriptome_conf_mapped_reads_frac',
        'multi_transcriptome_half_mapped_reads_frac',
        'multi_transcriptome_mapped_reads_frac',
        'multi_transcriptome_split_mapped_reads_frac',
        'multi_transcriptome_targeted_conf_mapped_reads_frac',
        'read_bases_with_q30_frac',
        'reads_per_cell',
        'total_cell_associated_partitions',
        'total_read_pairs',
        'total_singlets',
        'umi_bases_with_q30_frac',
        'unprocessed_read_pairs',
      ]
    ] +
    [
      lib_metrics_float(name, ['multi_success'])
      for name in [
        'samples_with_any_singlets',
        'sample_tags_count',
      ]
    ] +
    [
      {
        file_content: { json_path: ['sample_multiplexing_method'], path: ['${lib_metrics}'] },
        groups: ['multi_success'],
        name: 'sample_multiplexing_method',
        values: {
          match: [
            { exp: 'CMO' },
            { exp: 'Hashtag' },
            { exp: 'None' },
            { exp: 'OH' },
            { exp: 'RTL' },
          ],
        },
      },
    ] +
    [
      lib_metrics_float(name, ['count_success', 'multi_success'])
      for name in [
        'raw_matrix_barcodes_count',
        'raw_matrix_features_count',
        'raw_matrix_nonzero_entries_count',
        'filtered_matrix_barcodes_count',
        'filtered_matrix_features_count',
        'filtered_matrix_nonzero_entries_count',
        'raw_matrix_feature_type_count',
        'filtered_matrix_feature_type_count',
        'ANTIBODY_raw_matrix_feature_type_count',
        'ANTIBODY_filtered_matrix_feature_type_count',
        'ANTIGEN_raw_matrix_feature_type_count',
        'ANTIGEN_filtered_matrix_feature_type_count',
        'CRISPR_raw_matrix_feature_type_count',
        'CRISPR_filtered_matrix_feature_type_count',
        'Custom_raw_matrix_feature_type_count',
        'Custom_filtered_matrix_feature_type_count',
        'MULTIPLEXING_raw_matrix_feature_type_count',
        'MULTIPLEXING_filtered_matrix_feature_type_count',
        'raw_probe_matrix_features_count',
        'raw_probe_matrix_nonzero_entries_count',
        'raw_probe_matrix_gapfill_probes_count',
        'detected_hashtag_feature_count',
        'hashtag_reference_version',
      ]
    ] +
    [
      path_metric_float(name, '${vdj_lib_metrics}', [
        'library_metrics_vdj',
      ])
      for name in [
        'IGH_vdj_recombinome_mapped_reads_frac',
        'IGK_vdj_recombinome_mapped_reads_frac',
        'IGL_vdj_recombinome_mapped_reads_frac',
        'TRA_vdj_recombinome_mapped_reads_frac',
        'TRB_vdj_recombinome_mapped_reads_frac',
        'VDJ_bc_bases_with_q30_frac',
        'VDJ_read_bases_with_q30_frac',
        'VDJ_total_read_pairs',
        'VDJ_umi_bases_with_q30_frac',
        'VDJ_unprocessed_read_pairs',
        'multi_vdj_recombinome_mapped_reads_frac',
        'vdj_assemblable_read_pairs_per_filtered_bc',
        'vdj_filtered_bcs',
        'vdj_filtered_bcs_cum_frac',
        'vdj_good_bc_frac',
        'vdj_total_raw_read_pairs_per_filtered_bc',
      ]
    ] +
    [
      sample_metrics_float(name)
      for name in [
        'ANTIBODY_feature_reads_in_cells',
        'ANTIBODY_multi_filtered_bcs',
        'ANTIBODY_multi_filtered_bcs_mean_counts',
        'ANTIBODY_multi_filtered_bcs_median_counts',
        'ANTIBODY_multi_usable_reads_per_filtered_bc',
        'ANTIBODY_reads_lost_to_aggregate_GEMs',
        'ANTIBODY_recognized_feature_bc_frac_in_filtered_barcodes',
        'ANTIBODY_total_read_pairs_in_filtered_barcodes',
        'ANTIGEN_feature_reads_in_cells',
        'ANTIGEN_multi_filtered_bcs',
        'ANTIGEN_multi_filtered_bcs_mean_counts',
        'ANTIGEN_multi_filtered_bcs_median_counts',
        'ANTIGEN_multi_usable_reads_per_filtered_bc',
        'ANTIGEN_total_read_pairs_in_filtered_barcodes',
        'CRISPR_feature_bc_extracted_frac_in_filtered_barcodes',
        'CRISPR_feature_reads_in_cells',
        'CRISPR_frac_cells_with_multiple_protospacer',
        'CRISPR_frac_cells_with_protospacer',
        'CRISPR_multi_filtered_bcs',
        'CRISPR_multi_filtered_bcs_mean_counts',
        'CRISPR_multi_filtered_bcs_median_counts',
        'CRISPR_multi_usable_reads_per_filtered_bc',
        'CRISPR_recognized_feature_bc_frac_in_filtered_barcodes',
        'CRISPR_total_read_pairs_in_filtered_barcodes',
        'CRISPR_unrecognized_feature_bc_frac_in_filtered_barcodes',
        'Custom_feature_reads_in_cells',
        'Custom_multi_filtered_bcs',
        'Custom_multi_filtered_bcs_mean_counts',
        'Custom_multi_filtered_bcs_median_counts',
        'Custom_multi_usable_reads_per_filtered_bc',
        'Custom_recognized_feature_bc_frac_in_filtered_barcodes',
        'Custom_total_read_pairs_in_filtered_barcodes',
        'MULTIPLEXING_multi_filtered_bcs',
        'MULTIPLEXING_multi_filtered_bcs_mean_counts',
        'MULTIPLEXING_multi_filtered_bcs_median_counts',
        'MULTIPLEXING_multi_usable_reads_per_filtered_bc',
        'estimated_gdna_content',
        'estimated_gdna_unspliced_threshold',
        'filtered_bcs_mean_counts',
        'filtered_bcs_median_counts',
        'filtered_reads_per_filtered_bc',
        'median_genes_per_singlet',
        'median_umi_per_singlet',
        'multi_filtered_bcs_conf_mapped_barcoded_reads_cum_frac',
        'multi_transcriptome_conf_mapped_reads_frac_in_filtered_barcodes',
        'multi_transcriptome_half_mapped_reads_frac_in_filtered_barcodes',
        'multi_transcriptome_mapped_reads_frac_in_filtered_barcodes',
        'multi_transcriptome_split_mapped_reads_frac_in_filtered_barcodes',
        'multi_transcriptome_targeted_conf_mapped_reads_frac_in_filtered_barcodes',
        'sample_tags_count',
        'singlets_assigned_to_other_samples',
        'singlets_assigned_to_this_sample',
        'total_genes_detected',
        'total_read_pairs_in_filtered_barcodes',
        'total_singlets',
      ]
    ] +
    [
      {
        file_content: {
          json_path: [
            'sample_multiplexing_method',
          ],
          path: [
            '${sample_metrics}',
          ],
        },
        groups: [
          'sample_metrics',
        ],
        name: 'sample_multiplexing_method',
        values: {
          match: [
            {
              exp: 'None',
            },
            {
              exp: 'RTL',
            },
            {
              exp: 'OH',
            },
            {
              exp: 'CMO',
            },
            {
              exp: 'Hashtag',
            },
          ],
        },
      },
      {
        file_content: {
          json_path: [
            'sample_tags',
          ],
          path: [
            '${sample_metrics}',
          ],
        },
        groups: [
          'sample_metrics',
        ],
        name: 'sample_tags',
        values: {
          match: sample_tag_matches,
        },
      },
    ] +
    [
      path_metric_float(name, '${vdj_sample_metrics}', [
        'sample_metrics_vdj',
      ])
      for name in [
        'IGH_vdj_assembly_prod_cdr_bc_frac',
        'IGH_vdj_assembly_umis_per_cell_median',
        'IGK_IGH_vdj_assembly_contig_pair_productive_full_len_bc_frac',
        'IGK_vdj_assembly_prod_cdr_bc_frac',
        'IGK_vdj_assembly_umis_per_cell_median',
        'IGL_IGH_vdj_assembly_contig_pair_productive_full_len_bc_frac',
        'IGL_vdj_assembly_prod_cdr_bc_frac',
        'IGL_vdj_assembly_umis_per_cell_median',
        'TRA_TRB_vdj_assembly_contig_pair_productive_full_len_bc_frac',
        'TRA_vdj_assembly_prod_cdr_bc_frac',
        'TRA_vdj_assembly_umis_per_cell_median',
        'TRB_vdj_assembly_prod_cdr_bc_frac',
        'TRB_vdj_assembly_umis_per_cell_median',
        'TRD_vdj_assembly_prod_cdr_bc_frac',
        'TRD_vdj_assembly_umis_per_cell_median',
        'TRG_TRD_vdj_assembly_contig_pair_productive_full_len_bc_frac',
        'TRG_vdj_assembly_prod_cdr_bc_frac',
        'TRG_vdj_assembly_umis_per_cell_median',
        'multi_raw_vdj_paired_clonotype_diversity',
        'multi_vdj_assembly_contig_pair_productive_full_len_bc_count',
        'multi_vdj_assembly_contig_pair_productive_full_len_bc_frac',
        'vdj_filtered_bcs',
      ]
    ] +
    cell_annotation_metrics('${lib_metrics}', ['count_success', 'annotate_success']) +
    cell_annotation_metrics('${sample_metrics}', ['sample_metrics']) +
    // metrics related to the sequencing platform(s)
    [
      {
        file_content: {
          json_path: ['fastq_id_components'],
          path: ['${lib_metrics}'],
        },
        groups: count_multi_vdj_groups,
        name: 'GEX_fastq_id_components',
        values: { float: true },
      },
    ] +
    [
      lib_metrics_float(name)
      for name in [
        'ANTIBODY_fastq_id_components',
        'ANTIGEN_fastq_id_components',
        'CRISPR_fastq_id_components',
        'MULTIPLEXING_fastq_id_components',
        'Custom_fastq_id_components',
      ]
    ] +
    [
      path_metric_float(name, '${vdj_lib_metrics}', ['multi_success', 'vdj_success'])
      for name in [
        'VDJ_fastq_id_components',
      ]
    ] +
    [
      path_metric_float(name, '${vdj_lib_metrics}', ['multi_success'])
      for name in [
        'VDJ_B_fastq_id_components',
        'VDJ_T_fastq_id_components',
        'VDJ_T_GD_fastq_id_components',
      ]
    ] +
    // sequencer platform inferred from FASTQ headers
    [
      lib_metrics_preset('GEX_inferred_platform', 'inferred_platform', 'sequencer_platform'),
    ] +
    [
      lib_metrics_preset(name, name, 'sequencer_platform')
      for name in [
        'ANTIBODY_inferred_platform',
        'ANTIGEN_inferred_platform',
        'CRISPR_inferred_platform',
        'MULTIPLEXING_inferred_platform',
        'Custom_inferred_platform',
      ]
    ] +
    [
      {
        file_content: {
          json_path: [name],
          path: ['${vdj_lib_metrics}'],
        },
        groups: ['multi_success', 'vdj_success'],
        name: name,
        preset: 'sequencer_platform',
      }
      for name in [
        'VDJ_inferred_platform',
      ]
    ] +
    [
      {
        file_content: {
          json_path: [name],
          path: ['${vdj_lib_metrics}'],
        },
        groups: ['multi_success'],
        name: name,
        preset: 'sequencer_platform',
      }
      for name in [
        'VDJ_B_inferred_platform',
        'VDJ_T_inferred_platform',
        'VDJ_T_GD_inferred_platform',
      ]
    ],

  presets: shared.presets + shared.cr_sr_presets + {
    chemistry: {
      match: [
        named_re_matcher(k, insensitive=true)
        for k in [
          'auto',
          'custom',
          'threeprime',
          'fiveprime',
          'sc3pv1',
          'sc3pv2',
          'sc3pv3-polya',
          'sc3pv3-cs1',
          'sc3pv4-polya',
          'sc3pv4-cs1',
          'sc3pv4-polya-ocm',
          'sc3pv4-cs1-ocm',
          'sc3pv3ht-polya',
          'sc3pv3ht-cs1',
        ]
      ] + [
        {
          exp: '(?i)^%s$' % std.strReplace(k, '-', '-?'),
          key: k,
        }
        for k in [
          'sc5p-pe',
          'sc5p-pe-v3',
          'sc5p-pe-ocm-v3',
          'sc5p-r2',
          'sc5p-r2-v3',
          'sc5p-r2-ocm-v3',
          'sc5pht',
          'sc-fb',
        ]
      ] + [
        named_re_matcher(k, insensitive=true)
        for k in [
          'sfrp',
          'sfrp-no-trim-r2',
          'mfrp-rna',
          'mfrp-ab',
          'mfrp-crispr',
          'mfrp-47',
          'mfrp-96-r1',
          'mfrp-96-rna-r2',
          'mfrp-uncollapsed',
          'mfrp-rna-r1',
          'mfrp-ab-r1',
          'mfrp-ab-r2pos50',
          'mfrp-r1-48-uncollapsed',
          'mfrp-r1-no-trim-r2',
          'flex-v2-96-r1',
          'flex-v2-96-rna-r2',
          'flex-v2-r1',
          'flex-v2-rna-r2',
          'flex-v2-ab-r2:45',
          'flex-v2-ab-r2:64',
          'flex-v2-crispr-r2:1',
          'flex-v2-singleplex',
          'arc-v1',
          'scvdj',
          'scvdj-v3',
          'scvdj-v3-ocm',
          'scvdj-r2',
          'scvdj-r2-v3',
          'scvdj-r2-ocm-v3',
          'mfrp',
          'mfrp-r1',
          'flex-v2-r2',
          'sc3pv3',
          'sc3pv3ht',
          'sc3pv4',
          'sc3pv4ocm',
        ]
      ] + [
        {
          exp: '.*',
          key: 'unknown',
        },
      ],
    },
  },
}
