# Looks at the sample def and determines which feature-counter calls can be disabled

__MRO__ = """
stage DISABLE_FEATURE_STAGES(
    in  map[]  sample_def,
    out bool  disable_crispr,
    out bool  disable_antibody,
    src py     "stages/common/disable_feature_stages",
)
"""

def main(args, outs):
    sample_def = args.sample_def
    library_types = [x.get('library_type') for x in sample_def
                                if x.get('library_type') is not None]

    found_crispr = ('CRISPR Guide Capture' in library_types) or ('Gene Expression and CRISPR Guide Capture' in library_types)
    found_antibody = ('Antibody Capture' in library_types) or ('Gene Expression and Antibody Capture' in library_types)

    outs.disable_crispr = not(found_crispr)
    outs.disable_antibody = not(found_antibody)


