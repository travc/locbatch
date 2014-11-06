#!/usr/bin/python

__doc__ = """
    @TCC -- More usage documentation goes here"""

import sys, os
import csv
import logging
from collections import OrderedDict
from copy import deepcopy

### module level globals  ###################################
# setup logger in module level scope
logging.basicConfig(level=logging.INFO)
LOG = logging.getLogger(__name__)


class cEventDict(OrderedDict):
    """An OrderedDict which can be 'locked' so it throws error when trying to set non-existent attr
    Not terribly efficient, but useful for debugging / code maintenance"""
    _block_new_fields = False
    def __setitem__(self, name, value):
        if( self._block_new_fields and not name in self ):
            raise RuntimeError("setitem blocked for '{}'".format(name))
        else:
            super(cEventDict, self).__setitem__(name, value)
            
    def __setattr__(self, name, value):
        if( self._block_new_fields and not name in self.__dict__ ):
            raise RuntimeError("setattr blocked for '{}'".format(name))
        else:
            super(cEventDict, self).__setattr__(name, value)


### event_factory class ###################################

class event_factory(object):

    default_fields = OrderedDict([ # name, cast_function
        ['name', str],
        ['node_id', str],
        ['data_path', str],
        ['soundfile', str], # if time_is_global is True, then soundfile should be a metafile
        ['soundfile_pattern', str], # used to generate a rec_session if time_is_global == False
        ['channel', int], # should always be 1 for wired arrays
        ['start_time', float],
        ['duration', float],
        ['low_freq', float],
        ['high_freq', float],
        ['create_filter_order', int],
        ['smooth_cenv_window', float],
        ['temperature', float],
        ['RH', float],
        ['mics_filename', str],
        ['sum_extent', str],
        ['sum_resolution', str],
        ['refine_extent', str],
        ['refine_resolution', str],
        ['exclude_nodes', str],
        ['sum_computation_method', int],
        ['gcc_key_max_c_val_threshold', float],
        ['error_est_threshold', float],
        # values below should NOT be set in config
        ['idx', int],
        ['annotation_type', str],
        ['time_is_global', bool],
        ['rec_session_identifier', str],
        ])

    def __init__(self, defaults_dict, overrides_dict):
        self.defaults = cEventDict([[k, None] for k in self.default_fields])
        # update defaults with passed in defaults_dict
        self.defaults.update([ (k, self.default_fields[k](v))
                                if k in self.default_fields
                                else (k,v)
                                for k,v in dict(defaults_dict).iteritems() ])
        self.overrides = OrderedDict(overrides_dict)
        # @TCC debugging output
        if( False ):
            print "-- defaults",'-'*30
            for k,v in self.defaults.iteritems():
                print "{} : {} {}".format(k, v, type(v))
            print "-- overrides",'-'*30
            for k,v in self.overrides.iteritems():
                print "{} : {} {}".format(k, v, type(v))
            print '-'*45
            sys.exit(0)
        
    def apply_overrides(self, event):
        event.update([ (k, self.default_fields[k](v))
                            if k in self.default_fields
                            else (k,v)
                            for k,v in self.overrides.iteritems() ])


    ### Loading events from files #####################################

    def read_syrinx_annotations_file(self, events_filename, additional_defaults={}):
        """Read events from Syrinx annotations file"""
        LOG.info("Loading events from Syrinx file: '{}'".format(events_filename))
        events = []
        with open(events_filename, 'rU') as events_file:
            reader = csv.DictReader(events_file, delimiter='\t')
            for row in reader:
                # make fieldnames case insensitive
                fieldnames = [x.lower() for x in reader.fieldnames]
                tmp = {}
                for k,v in row.iteritems():
                    if( k.lower() in tmp ):
                        raise RuntimeError("Duplicate fieldname '{}' in annotation file '{}'\nNOTE: fieldnames should not be case-sensitive".format(k.lower(), events_filename))
                    tmp[k.lower()] = v
                row = tmp
                # syrinx annotation fields to event fields mapping
                ff2ef = {
                        'soundfile':'soundfile',
                        # channel -> node_id takes special treatment
                        'lefttimesec':'start_time',
                        # righttimesec -> duration takes special treatment
                        'bottomfreqhz':'low_freq',
                        'topfreqhz':'high_freq',
                        'rh':'RH',
                        }
                # fill out the event
                e = deepcopy(self.defaults)
                e._block_new_fields = True
                # apply additional defaults
                e.update(additional_defaults)
                e['annotation_type'] = 'syrinx'
                e['time_is_global'] = True # implies soundfile is really a metafile
                e['idx'] = len(events)+1 # 1-based index
                e['name'] = os.path.basename(events_filename)+' : '+str(e['idx'])
                # righttimesec -> duration takes special treatment
                e['duration'] = float(row['righttimesec'])-float(row['lefttimesec'])
                del row['righttimesec']
                # channel -> node_id takes special treatment
                e['node_id'] = str(int(row['channel'])).zfill(2) # channel column is really node id
                del row['channel']
                # channel is fixed at 1
                e['channel'] = 1
                # all other values in the row
                for ff in reader.fieldnames: # use fieldnames to keep in order
                    if( ff not in row ): # skip if we've already dealt with and removed
                        continue
                    v = row[ff]
                    if( v != '' ): # skip blank fields
                        if( ff in ff2ef ): # change the fieldname?
                            ef = ff2ef[ff]
                        else:
                            ef = ff
                        if( ef in self.default_fields ): # cast?
                            e[ef] = self.default_fields[ef](v) # cast type
                            if( e[ef] is None ):
                                raise RuntimeError("Failed to cast field {} value {} using {}".format(ef, v, self.default_fields[ef]))
                        else:
                            e._block_new_fields = False
                            e[ef] = v
                            e._block_new_fields = True
                # join to data_path if there is a data_path set
                if( e['data_path'] is not None ):
                    e['soundfile'] = os.path.normpath(os.path.join(e['data_path'], e['soundfile']))
                e['rec_session_identifier'] = row['soundfile'] # events with same rec_session_identifier use same instance
                # apply overrides
                self.apply_overrides(e)
                
                events.append(e)
        return events

    def read_raven_annotations_file(self, events_filename, additional_defaults={}):
        """Read events from Raven annotations (SelectionTable) file"""
        LOG.info("Loading events from Raven file: '{}'".format(events_filename))
        events = []
        with open(events_filename, 'rU') as events_file:
            reader = csv.DictReader(events_file, delimiter='\t')
            for row in reader:
                # make fieldnames case insensitive
                fieldnames = [x.lower() for x in reader.fieldnames]
                tmp = {}
                for k,v in row.iteritems():
                    if( k.lower() in tmp ):
                        raise RuntimeError("Duplicate fieldname '{}' in annotation file '{}'\nNOTE: fieldnames should not be case-sensitive".format(k.lower(), events_filename))
                    tmp[k.lower()] = v
                row = tmp
                # skip some entries
                if( not row['view'].lower().startswith('spectrogram') ): # only use spectrogram annotations
                    continue
                if( row['begin time (s)'] == row['end time (s)'] ): # skip point selections
                    continue
                # Raven annotation fields to event fields mapping
                ff2ef = {
                        'begin time (s)':'start_time',
                        # end time (s) -> duration takes special treatment
                        'low freq (hz)':'low_freq',
                        'high freq (hz)':'high_freq',
                        'rh':'RH',
                        }
                # fill out the event
                e = deepcopy(self.defaults)
                e._block_new_fields = True
                # apply additional defaults
                e.update(additional_defaults)
                e['idx'] = len(events)+1 # 1-based index
                e['name'] = os.path.basename(events_filename)+' : '+str(row['selection'])
                e['annotation_type'] = 'raven'
                e['time_is_global'] = False
                # end time (s) -> duration takes special treatment
                e['duration'] = float(row['end time (s)'])-float(row['begin time (s)'])
                del row['end time (s)']
                # fill in rest of fields
                for ff in fieldnames: # use fieldnames to keep in order
                    if( ff not in row ): # skip if we've already dealt with and removed
                        continue
                    v = row[ff]
                    if( v != '' ): # skip blank fields
                        if( ff in ff2ef ): # change the fieldname?
                            ef = ff2ef[ff]
                        else:
                            ef = ff
                        if( ef in self.default_fields ): # cast?
                            e[ef] = self.default_fields[ef](v) # cast type
                        else: # field is not something we use, but add it anyway for output
                            e._block_new_fields = False
                            e[ef] = v
                            e._block_new_fields = True
                # if soundfile is default, then prefix the path so rec_session works properly
                if( not 'soundfile' in reader.fieldnames and e['data_path'] is not None ):
                    e['soundfile'] = os.path.realpath(os.path.join(e['data_path'], e['soundfile']))
                # fill out meta fields based on info read in and defaults
                e['rec_session_identifier'] = e['soundfile'] # events with same rec_session_identifier can use same rec_session instance
                # apply overrides
                self.apply_overrides(e)
                events.append(e)
                assert e['mics_filename'] is not None
        return events
