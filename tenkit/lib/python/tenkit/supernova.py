#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
# Support code for Supernova stages
#
# "Pseudo" Martian API for reporting errors and warnings from C++
# code prior to deployment of a formal adapter.
#
# Code to plot molecule length histogram and the kmer spectrum

import os
import json
import martian
from constants import ALARMS_LOCATION

def check_alert(stage, alert):
    keys = ["action", "metric", "compare", "threshold", "message"]
    for key in keys:
        if not alert.has_key(key):
           print key, " is missing in "
           print alert
           martian.throw("incorrectly formatted alert, see stdout.")
    if not ( alert["compare"] == "<" or alert["compare"] == ">" ):
        print alert
        martian.throw("invalid value for compare in alert")
    if not ( type(alert["threshold"]) == int or type(alert["threshold"]) == float ):
        martian.throw("%s: invalid type for threshold" % type(alert["threshold"]))

def load_alerts():
    alerts_file = os.path.join(ALARMS_LOCATION, "alarms-supernova.json")
    json_string = open( alerts_file, "r" ).read()
    try:
        alerts = json.loads(json_string)
    except:
        martian.throw("Incorrectly formatted alarms-supernova.json file.")
    for stage, alert_list in alerts.iteritems():
        for alert in alert_list:
            check_alert(stage, alert)
    return alerts

def write_stage_alerts(stage, path, alerts_file="alerts.list"):
    alerts = load_alerts()
    out_file = os.path.join(path, alerts_file)
    if not os.path.exists(path):
        os.makedirs(path)
    out_handle = open(out_file, "w")
    keys = ["metric", "threshold", "compare", "action", "message"]
    if not alerts.has_key(stage):
        martian.throw("No alerts found for stage %s" % stage)
    for alert in alerts[stage]:
        out_handle.write("#\n")
        out_handle.write(alert["metric"]+"\n")
        out_handle.write(str(alert["threshold"])+"\n")
        out_handle.write(alert["compare"]+"\n")
        out_handle.write(alert["action"]+"\n")
        out_handle.write(alert["message"]+"\n") 
    out_handle.close()

class AlertLogger:
    def __init__(self, stage):
        alerts_all = load_alerts()
        self.alerts = alerts_all[stage]

    def issue(self, metric, value, format_string=""):
        for alert in self.alerts:
            ## find the right metric
            if alert["metric"] == metric:
                ## should we trigger?
                if (alert["compare"] == ">") ^ (value < alert["threshold"]):
                    ## optional formatting of alert message with format_string or value
                    if len(format_string) == 0:
                        format_string = str(value)
                    message = alert["message"].replace("{}", format_string)
                    ## issue an alert
                    if alert["action"] == "alarm":
                        martian.alarm(message)
                    elif alert["action"] == "exit":
                        martian.exit(message)

class SupernovaAlarms:
    SN_STAGE_ALARMS  = "martian_alerts.json"
    SN_ROLLUP_ALARMS = "alerts_rollup.txt"
    SN_EXIT  = u"exit"
    SN_ALARM = u"alarm"
    SN_ALARM_HEAD      = "The following warning(s) were issued prior to encountering an error:"
    SN_UNEXPECTED_TEXT = "An unexpected error has occurred."
    SN_ALERT_HANDLERS={ SN_ALARM    : martian.alarm, \
                        u"log_info" : martian.log_info, \
                        u"log_warn" : martian.log_warn, \
                        u"throw"    : martian.throw, \
                        SN_EXIT     : martian.exit }

    def __init__(self, base_dir, 
                 alarms_file = SN_STAGE_ALARMS, rollup_file = SN_ROLLUP_ALARMS, 
                 delete = True ):
        self._alarms_file=os.path.join(base_dir, alarms_file)
        self._rollup_file=os.path.join(base_dir, rollup_file)
        self._posted=False
        if delete: self.check_delete()

    def exit(self, msg=None):
        exit_handler=self.SN_ALERT_HANDLERS[self.SN_EXIT]
        if msg is None:
            full_msg = self.SN_UNEXPECTED_TEXT
        else:
            full_msg = msg
        if os.path.exists(self._rollup_file):
            rollup = open(self._rollup_file, "r").read()
            full_msg += "\n\n"
            full_msg += self.SN_ALARM_HEAD
            full_msg += "\n\n"
            full_msg += rollup
        exit_handler(full_msg)

    def post(self):
        self._posted=True
        
        ## if alarm file does not exist, do nothing
        if os.path.exists(self._alarms_file):
            with open(self._alarms_file, 'r') as fp:
                alerts=json.loads(fp.read())
        else:
            return

        meta_alarm=[]
        exit_str=''

        for k,v in alerts.iteritems():
            if not k in self.SN_ALERT_HANDLERS:
                meta_alarm.append("unknown key {} in {} (BUG)".format(k,self._alarms_file))
            elif k == self.SN_EXIT:
                exit_str=';'.join(v)
            else:
                handler=self.SN_ALERT_HANDLERS[k]
                for post in v:
                    handler( post )

        for alarm in meta_alarm:
            martian.alarm(alarm)
        
        if len(exit_str) > 0:
            self.exit(exit_str)

    def __del__(self):
        if not self._posted: 
            martian.alarm("C++ alarms were not posted, but left in file (BUG).")
        else:
            self.check_delete()

    def check_delete(self, filename = None):
        if filename == None: 
            filename = self._alarms_file
        if os.path.isfile( filename ):
            os.unlink(filename)


## this function plots a normalized histogram
def plot_norm_pct_hist( plt, values, binsize, start, **plt_args ):
    x = start
    xvals = []
    yvals = []
    norm = 0.0
    for v in values:
        xvals.append(x)
        yvals.append(v)
        xvals.append(x+binsize)
        norm += v
        yvals.append(v)
        x += binsize
    for i in xrange (len(yvals)):
        yvals[i] = yvals[i]/norm*100.0
    plt.plot( xvals, yvals, **plt_args)
    plt.xlim( start, x )

def try_plot_molecule_hist(args):
    try:
        import matplotlib as mpl
        mpl.use('Agg')
        import matplotlib.pyplot as plt
        import json
        import numpy as np
        from statsmodels.nonparametric.smoothers_lowess import lowess

        ## final plot name
        plot_name = os.path.join( args.parent_dir, "stats", 
                                  "molecule_lengths.pdf" )
        mol_fn = os.path.join( args.parent_dir, "stats", 
                               "histogram_molecules.json" )
        mol_dict = json.loads(open( mol_fn, "r" ).read())
        if mol_dict["numbins"] > 0:
            xmin = 0     ## min length kb
            xmax = 300   ## max length kb
            binsize=1    ## bin size of hist in kb
            ## compute length-weighted histogram
            lwhist = []
            for x, v in enumerate( mol_dict["vals"] ):
               lwhist.append( (x + 0.5)* v )
            ## truncate
            lwhist = lwhist[:xmax]
            ## normalize
            norm = sum(lwhist)
            lwhist = np.array([100*x/norm for x in lwhist])
            ## lowess smoothing
            xvalues = np.arange(0, xmax, 1) + 0.5
            newhist = lowess(lwhist, xvalues, frac=0.1, return_sorted=False)
            ## do the plotting
            fig, ax = plt.subplots()
            ## set up axis, grid
            ax.grid(True)
            plt.locator_params( 'x', nbins=10 )
            plt.locator_params( 'y', nbins=10 )
            plt.plot ( newhist, **{"lw": 2, "ls": "-", "color": "blue"} ) 
            plt.xlim( xmin, xmax )
            plt.ylim( 0, None )
            plt.xlabel ( "Inferred molecule length (kb)")
            plt.ylabel ( "Percent of total DNA mass (1 kb bins)")
            plt.title ("Length-weighted histogram (LOWESS-smoothed)")
            plt.savefig( plot_name )
            plt.close()
    except ImportError, e:
        martian.log_info( "Error importing libraries for plotting" )
        martian.log_info( str(e) )
    except KeyError, e:
        martian.log_info( "Invalid key in json while plotting" )
        martian.log_info( str(e) )
    except IOError, e:
        martian.log_info( "Could not find the molecule json for plotting" )
        martian.log_info( str(e) )
    except Exception as e:
        martian.log_info( "Unexpected error while plotting molecule length")
        martian.log_info( str(e) )

## take a set of (distinct) numbers and format them into a list of strings
## where we have an integer followed by a suffix
## picks out the representation that has the shortest length
def nice_labels ( numbers ):
    suffixes = ['', 'K', 'M', 'G']
    suff_len = []
    ## figure out which suffix gives us the shortest label length
    for i, suff in enumerate( suffixes ):
        test   = [float(y)/(1000.0**i) for y in numbers]
        labels = ["%d%s"% (int(y), suff) for y in test]
        ## make sure that in the new representation there are no
        ## degenerate cases
        if len(set(labels)) == len(labels):
            suff_len.append( (sum(map(len, labels)), i) )
    ## if we fail to find any satisfactory suffixes, just use defaults
    if len(suff_len) == 0:
        return map(str, numbers), 0
    else:
        suff_len.sort()
        i = suff_len[0][1]
        labels = ["%d%s"% (int(float(y)/(1000.0**i)), suffixes[i]) for y in numbers]
    return labels, i

def try_plot_kmer_spectrum(args):
    try:
        import matplotlib as mpl
        mpl.use('Agg')
        import matplotlib.pyplot as plt
        import json

        ## final plot name
        plot_name = os.path.join( args.parent_dir, "stats", 
                                  "kmer_spectrum.pdf" )
        kmer_fn = os.path.join( args.parent_dir, "stats", 
                               "histogram_kmer_count.json" )
        data = json.loads(open( kmer_fn, "r" ).read())
        ## max k-mer multiplicity
        MAX = 100 
        if len( data["vals"] ) == 0:
            martian.log_info ("No kmer data to plot.")
            return
        elif len( data["vals"] ) < MAX:
            martian.log_info ("Warning: no kmers with multiplicity >= %d"%MAX )
            MAX = len( data["vals"] )
        fig, ax = plt.subplots()
        #plt.setp(ax.get_yticklabels(), visible=False)
        ## set mode to 1.0
        xvals = range(MAX)
        yvals = data["vals"][:MAX]
        ax.plot (xvals, yvals, lw = 2.0, color="blue" )
        ## plot tail
        tail_height = float(sum(data["vals"][MAX:]))
        _, ymax = plt.ylim()
        plt.axvspan (xmin=MAX-1, xmax=MAX, ymin=0, ymax=tail_height/ymax, ls=None, color="blue")
        ## set up axis, grid
        ax.grid(True)
        plt.locator_params( 'x', nbins=10 )
        plt.locator_params( 'y', nbins=10 )
        plt.xlim(0, MAX)
        yt = ax.get_yticks()
        ylabels, yexp = nice_labels( yt )
        plt.yticks ( yt, ylabels, rotation=45 )
        plt.xlabel ( "filtered kmer multiplicity" )
        plt.ylabel ( "counts" )
        plt.savefig (plot_name)
        plt.close()
    except ImportError, e:
        martian.log_info( "Error importing libraries for plotting" )
        martian.log_info( str(e) )
    except KeyError, e:
        martian.log_info( "Invalid key in json while plotting" )
        martian.log_info( str(e) )
    except IOError, e:
        martian.log_info( "Could not find the molecule json for plotting" )
        martian.log_info( str(e) )
    except Exception as e:
        martian.log_info( "Unexpected error while plotting molecule length")
        martian.log_info( str(e) )


