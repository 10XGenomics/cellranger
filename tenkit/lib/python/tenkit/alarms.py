#!/usr/bin/env python
#
# Copyright (c) 2014 10X Genomics, Inc. All rights reserved.
#
# Utilities for reporting alarms based on sequencing data
#

import os
import json
from tenkit.constants import ALARMS_LOCATION


def eval_rule(alarm, metrics):
    '''Evaluate an alarm, against a dict of metrics'''

    alarm = alarm.copy()

    try:

        # Raw value of metric -- loupe may want to do custom formatting
        # raw_value must be supplied
        alarm['raw_value'] = eval(alarm['raw_value'], {}, metrics)

        # None values of metrics cannot cause an alarm
        # to be raised -- the condition must be detected
        # by another metric
        if alarm['raw_value'] is None:
            alarm['raised'] = False
            return alarm

        # Test if we've tripped the alarm
        raised = eval(alarm["test"], {}, metrics)

        # Format the message
        message = alarm["message"] % metrics

        # Custom formatting - either supply a code snippet, a format string, or just get the default python conversion
        if alarm.has_key("formatted_value"):
            alarm['formatted_value'] = eval(alarm['formatted_value'], {}, metrics)
        elif alarm.has_key('format'):
            alarm['formatted_value'] = alarm['format'] % alarm['raw_value']
        else:
            alarm['formatted_value'] = str(alarm['raw_value'])

        alarm["raised"] = raised
        alarm["message"] = message

    except Exception as e:
        print "error evaluating alarm: %s" % str(e)
        alarm["raised"] = False

    return alarm


def evaluate_alarms(alarms, metrics):
    ''' Evaluate the alarms against the metrics, prune alarms that are superseded by another
        alarm, and return the list of 'active' alarms '''

    raised = []
    alarm_dict = {x['id']:x for x in alarms}

    def alarm_all_parents(name):
        all_parents = []
        parent = alarm_dict[name].get('parent', None)

        if parent is not None:
            if type(parent) is list:
                my_parents = parent
            else:
                my_parents = [parent]

            for p in my_parents:
                all_parents.append(p)
                all_parents.extend(alarm_all_parents(p))

        return all_parents

    for a in alarms:
        a = eval_rule(a, metrics)
        if a["raised"]:
            raised.append(a)

        # Check that all the parents exist
        alarm_all_parents(a['id'])


    raised_ids = [r['id'] for r in raised]

    # filter out alarms superseded by other alarms
    filtered_raised = []

    for r in raised:
        parents = alarm_all_parents(r['id'])

        # Check if any transitive parents were raised
        # If so, don't report this alarm
        if any(p in raised_ids for p in parents):
            continue
        else:
            filtered_raised.append(r)

    return filtered_raised


def load_rules(targets):
    ''' Select the appropriate alarm rules to use, and load '''

    rules = []

    common_fn = os.path.join(ALARMS_LOCATION, "alarms-common.json")
    with open(common_fn) as common_file:
        rules.extend(json.load(common_file))

    if targets is None:
        alarms_fn = "alarms-wgs.json"
    else:
        alarms_fn = "alarms-targeted.json"

    fn = os.path.join(ALARMS_LOCATION, alarms_fn)
    if not os.path.exists(fn):
        raise NameError("unable to find alarms file: %s" % fn)

    with open(fn) as rules_file:
        rules.extend(json.load(rules_file))

    return rules
