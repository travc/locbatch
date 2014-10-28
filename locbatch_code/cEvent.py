#!/usr/bin/python

__doc__ = """cEvent class for holding annotation information
    and methods to read and write them from files
    @TCC -- More usage documentation goes here"""

import sys, os
import copy
import copy
import csv # for reading event/annotation files
import ConfigParser
import logging

# Do we have OrderedDict?
try:
    from ordereddict import OrderedDict
    HAVE_OrderedDict = True
except ImportError:
    HAVE_OrderedDict = False

import mppool # simplified multiprocessing interface
from mppool import mp_pool
from misc_utils import verbose, ConfigGet, StrOrNone, CastOrNone
import rec_session # for looking up timing if needed

### module level globals  ###################################

# a ConfigParser for default values...
# must be set elsewhere
config = None

#### setup logger in module level scope ####
logging.basicConfig(level=logging.INFO)
LOG = logging.getLogger(__name__)

### Event class, just a container  ###################################

class cEvent():
    ## necessiary variables ##
    annotation_type = None # 'syrinx' or eventually 'raven', 'praat', ect.
    time_is_global = True  # is the start_time relative to start of rec_session or file
    rec_session_identifier = None # events with same rec_session_identifier can use same rec_session instance

    name = None # just a label
    idx = None # numeric identifier

    start_time = None # start time of event in seconds (assuming nominal FS)
    duration = None # duration in seconds (assuming nominal FS)

    node_id = None
    channel = None

    filter_desc = None # None is acceptable value
    nodeloc_filename = None # file with microphone locations

    temperature = None # temperature in C required for computing speed of sound
    RH = None  # relative humidity in % required for computing speed of sound

    # any addtional fields/values included in the orginal annotation file (a dictionary)
    values_from_file = {}


    def __str__(self):
        """just print all the members, for debugging really"""
        state = ["%s=%r" % (attribute, value) for (attribute, value) in self.__dict__.items()]
        return '\n'.join(state)


    def Get(self, attr_name, default=True):
        """Get the attribute or if attribute is none (default) return default"""
        rv = None
        try:
            rv = getattr(self,attr_name)
        except AttributeError:
            pass
        if( rv == None ):
            rv = self.GetDefaultForAttr(attr_name)
        return rv

    def GetDefaultForAttr(self, attr_name, config_section='Event defaults'):
        assert( config != None )
        rv = None
        if( attr_name == 'filter_desc' ):
            rv = self.GenerateDefaultFilterDesc()
        elif( attr_name == 'node_id' ):
            rv = config.get('Event defaults','node_id')
        elif( attr_name == 'channel' ):
            rv = config.getint('Event defaults','channel')
        elif( attr_name == 'temperature' ):
            rv = config.getfloat('Event defaults','temperature')
        elif( attr_name == 'RH' ):
            rv = config.getfloat('Event defaults','RH')
        elif( attr_name == 'nodeloc_filename' ):
            rv = config.get('Event defaults','nodeloc_filename')
        elif( attr_name == 'exclude_nodes' ):
            rv = config.get('Event defaults','exclude_nodes')
        return rv

    def GenerateDefaultFilterDesc(self):
        if( config.getboolean('Event defaults', 'create_filter_from_freq_range') ):
            order = config.getint('Event defaults', 'create_filter_order')
            if( order <= 0 ):
                LOG.error("create_filter_order must be > 0")
                sys.exit(2)
            return '*butter {0:d} band {1:.1f} {2:.1f}'.format(order, self.low_freq, self.high_freq)
        else:
            return StrOrNone(config_get_type(config, 'Event defaults', 'default_filter'))

    def FillInDefaults(self, config):
        #name = None  # no default, not required
        #start_time = None # no default, required
        assert( self.start_time != None )
        #duration = None  # no default_ required
        assert( self.duration != None )
        #node_id = None # default, required
        self.node_id = self.Get('node_id')
        assert( self.node_id != None )
        #channel = None # default, required
        self.channel = self.Get('channel')
        assert( self.channel != None )
        #filter_desc = None # default, not required
        self.filter_desc = self.Get('filter_desc')
        #nodeloc_filename = None # default, required
        self.nodeloc_filename = self.Get('nodeloc_filename')
        assert( self.nodeloc_filename != None )
        #temperature = None # default, required
        self.temperature = self.Get('temperature')
        assert( self.temperature != None )
        #RH = None  # default, required
        self.RH = self.Get('RH')
        assert( self.RH != None )
        self.exclude_nodes = self.Get('exclude_nodes')
        ## overrides from config
        if( config.has_option('Event overrides', 'min_freq') ):
            value = config.getfloat('Event overrides', 'min_freq')
            if( self.low_freq < value ):
                self.low_freq = value
        if( config.has_option('Event overrides', 'max_freq') ):
            value = config.getfloat('Event overrides', 'max_freq')
            if( self.high_freq > value ):
                self.high_freq = value


### Loading events from files #####################################

def ReadSyrinxAnnotationsFile(cfg, events_filename):
    """Read events from Syrinx annotations file"""
    global config
    config = cfg
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
            # fill out the event
            e = cEvent()
            e.values_from_file = {}
            e.annotation_type = 'syrinx'
            e.time_is_global = True
            e.rec_session_identifier = row['soundfile'] # events with same rec_session_identifier can use same rec_session instance
            e.idx = len(events)+1 # 1-based index
            e.name = os.path.basename(events_filename)+' : '+str(e.idx)
            e.metafile = row['soundfile']
            del row['soundfile']
            e.start_time = float(row['lefttimesec'])
            del row['lefttimesec']
            e.duration = float(row['righttimesec'])-e.start_time
            del row['righttimesec']
            e.node_id = str(int(row['channel'])).zfill(2) # channel column is really node id
            del row['channel']
            e.channel = 1 # @TCC channel fixed at 1
            e.low_freq = float(row['bottomfreqhz'])
            del row['bottomfreqhz']
            e.high_freq = float(row['topfreqhz'])
            del row['topfreqhz']
            #nodeloc_filename = None
            if( 'temperature' in fieldnames ):
                e.temperature = row['temperature']
                del row['temperature']
            if( 'rh' in fieldnames ):
                e.RH = row['rh']
                del row['rh']
            # any addtional fields/values included in the orginal annotation file (a dictionary)
            e.values_from_file.update(row)
            # fill in defaults
            e.FillInDefaults(config)
            # create the filter
            e.filter_desc = e.GenerateDefaultFilterDesc() # default filter based on high and low freq
            events.append(e)
    return events

def ReadRavenAnnotationsFile(cfg, events_filename):
    """Read events from Raven selection table file"""
    # expected columns: Selection       View    Channel Begin Time (S)  End Time (S)    Low Freq (Hz)   High Freq (Hz)
    global config
    config = cfg
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
            # fill out the event
            e = cEvent()
            e.values_from_file = {}
            e.annotation_type = 'raven'
            e.time_is_global = False
            e.rec_session_identifier = cfg.get('Files','soundfile_pattern') # events with same rec_session_identifier can use same rec_session instance
            e.soundfile_pattern = cfg.get('Files','soundfile_pattern')
            e.idx = len(events)+1 # 1-based index
            e.name = os.path.basename(events_filename)+' : '+str(row['selection'])
            # optional important columns
            if( 'node_id' in fieldnames ):
                e.node_id = row['node_id']
                del row['node_id']
            else:
                e.node_id = cfg.get('Files','node_id')
            if( 'soundfile' in fieldnames ):
                e.soundfile = row['soundfile']
                del row['soundfile']
            else:
                e.soundfile = cfg.get('Files','soundfile')
            if( 'temperature' in fieldnames ):
                e.temperature = row['temperature']
                del row['temperature']
            if( 'rh' in fieldnames ):
                e.RH = row['rh']
                del row['rh']
            # required columns
            e.start_time = float(row['begin time (s)'])
            del row['begin time (s)']
            e.duration = float(row['end time (s)'])-e.start_time
            del row['end time (s)']
            e.channel = int(row['channel'])
            assert e.channel == 1, "Only 1 channel recordings supported for now"
            del row['channel']
            e.low_freq = float(row['low freq (hz)'])
            del row['low freq (hz)']
            e.high_freq = float(row['high freq (hz)'])
            del row['high freq (hz)']
            #nodeloc_filename = None
            # any addtional fields/values included in the orginal annotation file (a dictionary)
            e.values_from_file.update(row)
            # fill in defaults
            e.FillInDefaults(config)
            # create the filter
            e.filter_desc = e.GenerateDefaultFilterDesc() # default filter based on high and low freq
            events.append(e)
    return events



#####################################################################
### DEPRICATED ###
#    def FreqsFromFilterDesc(self, FS,
#            freq_range_lower_padding=0, freq_range_upper_padding=0, default=False):
#        freq_range = [0.0, FS/2.0-1]
#        if( self.filter_desc ):
#            filter_desc = self.filter_desc
#        elif( default ): # used default if default flag is set
#            filter_desc = self.GenerateDefaultFilterDesc()
#        if( filter_desc ):
#            # partially parse the filter_desc, @TCC should probably move into a loc_misc function
#            if( filter_desc[0] == '*' ): # auto generated filter, can get type and freq
#                l = filter_desc.split()
#                type = l[2]
#                freq = [ float(v) for v in l[3:] ]
#                if( type == 'band' ): freq_range = [freq[0]-freq_range_lower_padding, freq[1]+freq_range_upper_padding]
#                elif( type == 'low' ): freq_range = [0, freq[0]+freq_range_upper_padding]
#                elif( type == 'high' ): freq_range = [freq[0]-freq_range_lower_padding, FS/2.0]
#                if( freq_range[0] < 0 ): freq_range[0] = 0
#                if( freq_range[1] > FS/2.0 ): freq_range[1] = FS/2.0
#        return freq_range
