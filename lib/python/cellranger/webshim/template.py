#!/usr/bin/env python

import argparse
import jinja2
import json
import os
import sys
import tenkit.safe_json as tk_safe_json
import cellranger.webshim.lz_string as lz_string
import cellranger.webshim.constants.shared as constants

DEFAULT_TEMPLATE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "template")

def get_template_for_pipeline(pipeline, data):
    if pipeline == constants.PIPELINE_COUNT:
        is_barnyard = len(data['genomes']) > 1
        if is_barnyard: return 'count-multi-genome.html'
        else: return 'count-single-genome.html'
    elif pipeline == constants.PIPELINE_AGGR:
        is_barnyard = len(data['genomes']) > 1
        if is_barnyard: return 'aggr-multi-genome.html'
        else: return 'aggr-single-genome.html'
    elif pipeline == constants.PIPELINE_REANALYZE:
        return 'reanalyze.html'
    elif pipeline == constants.PIPELINE_VDJ:
        return 'vdj.html'
    else:
        raise ValueError("Invalid pipeline type: %s" % pipeline)

def convert_webshim_json_to_html(data, pipeline, template_dir=None):
    if template_dir is None:
        template_dir = DEFAULT_TEMPLATE_DIR
    loader = jinja2.FileSystemLoader(template_dir)
    env = jinja2.Environment(loader=loader, trim_blocks=True, lstrip_blocks=True,
                             variable_start_string='[[',
                             variable_end_string=']]')
    env.globals['include_file'] = lambda name: loader.get_source(env, name)[0]
    template_html = get_template_for_pipeline(pipeline, data)
    template = env.get_template(template_html)
    compressed_data = lz_string.compressToEncodedURIComponent(json.dumps(tk_safe_json.json_sanitize(data)))
    return template.render(data=data, js_compressed_data=compressed_data).encode('utf-8')

def main():
    parser = argparse.ArgumentParser(description='Build webshim html')
    parser.add_argument('--template-dir', type=str, help='template HTML dir')
    parser.add_argument('--pipeline-type', type=str, help='pipeline type')
    args = parser.parse_args()

    data = json.loads(sys.stdin.readline().strip())

    print convert_webshim_json_to_html(data, pipeline=args.pipeline_type, template_dir=args.template_dir)

if __name__ == '__main__':
    main()
