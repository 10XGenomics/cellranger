#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

from collections import namedtuple, OrderedDict
import h5py
import cellranger.io as cr_io

FeatureDef = namedtuple('FeatureDef', ['index', 'id', 'name', 'feature_type', 'tags'])

# Required HDF5 datasets
REQUIRED_DATASETS = ['id', 'name', 'feature_type']

class FeatureDefException(Exception):
    pass

class FeatureReference(object):
    '''Store a list of features (genes, antibodies, etc).'''

    def __init__(self, feature_defs, all_tag_keys):
        '''Create a FeatureReference.

        Args:
            feature_defs (list of FeatureDef): All feature definitions.
            all_tag_keys (list of str): All optional tag keys.
        '''
        self.feature_defs = feature_defs
        self.all_tag_keys = all_tag_keys

        # Assert uniqueness of feature IDs
        id_map = {}
        for fdef in self.feature_defs:
            if fdef.id in id_map:
                this_fd_str = 'ID: %s; name: %s; type: %s' % (fdef.id, fdef.name, fdef.feature_type)
                seen_fd_str = 'ID: %s; name: %s; type: %s' % (id_map[fdef.id].id, id_map[fdef.id].name, id_map[fdef.id].feature_type)
                raise FeatureDefException('Found two feature definitions with the same ID: (%s) and (%s). All feature IDs must be distinct.' %
                                          (this_fd_str, seen_fd_str))
            id_map[fdef.id] = fdef

        self.id_map = id_map

    @staticmethod
    def addtags(feature_ref, new_tags, new_labels=None):
        '''Add new tags and corresponding labels to existing feature_ref
           If new labels are None, empty strings are supplied by default
        
        Args:
            feature_ref: a FeatureReference instance
            new_tags: a list of new tags
            new_labels: per feature list of label values corresponding to the new tags
        '''
        assert len(new_tags) > 0
        for tag in new_tags:
            assert tag not in feature_ref.all_tag_keys

        use_labels = []
        if new_labels is not None:
            assert len(feature_ref.feature_defs) == len(new_labels)
            for labels in new_labels:
                assert len(labels) == len(new_tags)
            use_labels = new_labels
        else:
            # initialize to empty
            for i in range(len(feature_ref.feature_defs)):
                use_labels += [[''] * len(new_tags)]
        assert len(feature_ref.feature_defs) == len(use_labels)

        augmented_features = []
        for fd, newvals in zip(feature_ref.feature_defs, use_labels):
            A = {a:b for a,b in zip(new_tags, newvals)}
            A.update(fd.tags)
            augmented_features.append(FeatureDef(index=fd.index,
                                      id=fd.id,
                                      name=fd.name,
                                      feature_type=fd.feature_type,
                                      tags=A))

        return FeatureReference(feature_defs=augmented_features,
                                all_tag_keys=feature_ref.all_tag_keys + new_tags)

    @staticmethod
    def join(feature_ref1, feature_ref2):
        '''Concatenate two feature references, requires unique ids and identical tags'''
        assert feature_ref1.all_tag_keys == feature_ref2.all_tag_keys
        feature_defs1 = feature_ref1.feature_defs
        feature_defs2 = feature_ref2.feature_defs
        return FeatureReference(feature_defs=feature_defs1 + feature_defs2, all_tag_keys=feature_ref1.all_tag_keys)

    @classmethod
    def empty(cls):
        return cls(feature_defs=[], all_tag_keys=[])

    def get_num_features(self):
        return len(self.feature_defs)

    def to_hdf5(self, group):
        '''Write to an HDF5 group.'''

        # Write required datasets
        for col in REQUIRED_DATASETS:
            data = [getattr(f, col) for f in self.feature_defs]
            cr_io.create_hdf5_string_dataset(group, col, data)

        # Write tag datasets
        for col in self.all_tag_keys:
            # Serialize missing data as empty unicode string
            data = [f.tags.get(col, '') for f in self.feature_defs]
            cr_io.create_hdf5_string_dataset(group, col, data)

        # Record names of all tag columns
        cr_io.create_hdf5_string_dataset(group, '_all_tag_keys', self.all_tag_keys)

    @classmethod
    def from_hdf5(cls, group):
        '''Load from an HDF5 group.

        Args:
            group (h5py.Dataset): Group to load from.
        Returns:
            feature_ref (FeatureReference): New object.
        '''
        data = OrderedDict()

        # Load HDF5 datasets
        # FIXME: ordering may not be guaranteed in python3
        for name in group:
            node = group[name]
            if node.shape is None:
                # Empty datset
                data[name] = []
            elif isinstance(node, h5py.Dataset):
                if node.dtype.char == 'S':
                    data[name] = cr_io.read_hdf5_string_dataset(node)
                else:
                    data[name] = node[:]

        all_tag_keys = data['_all_tag_keys']

        # Build FeatureDefs
        feature_defs = []
        num_features = len(data[REQUIRED_DATASETS[0]])

        for i in xrange(num_features):
            fdef_dict = {col: data[col][i] for col in REQUIRED_DATASETS}
            tags = {k: data[k][i] for k in all_tag_keys if len(data[k][i]) > 0}
            fdef_dict['index'] = i
            fdef_dict['tags'] = tags
            feature_defs.append(FeatureDef(**fdef_dict))

        return cls(feature_defs, all_tag_keys=all_tag_keys)
