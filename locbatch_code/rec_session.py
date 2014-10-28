#!/usr/bin/python

__doc__ = 'Some usage documentation goes here'

import os,sys
from sys import stderr,stdout
import subprocess
import glob
import logging

import wave

import numpy as NP # for array methods returning data, could use pylab instead of numpy
import string # for array methods returning data
import re # for parsing filenames
from collections import OrderedDict

#### setup logger in module level scope ####
logging.basicConfig(level=logging.INFO)
LOG = logging.getLogger(__name__)

### global constants (wrapped in simple class) @TCC-- Move into appropriate class? ###
class sConstants:
    FILE_TYPE_RAW = 0
    FILE_TYPE_WAV = 1


### AUDIO FILE UTILITY FUNCTIONS ############################################
def GetAudioFileParams(filename, config):
    """Using sox...
        returns (file_type, FS, bytes_per_sample, num_channels, num_samps)
        num_samps is aka 'frames', 1 samp includes all channels
        bytes_per_sample is for 1 sample from 1 channel, not a full samp/frame"""
    if( filename[-4:].lower() == '.wav' ): # Try native handling for wav files
        return GetAudioFileParams_native(filename, config)
    output_unhandled = False
    SOX_EXEC = config.get('System',sys.platform+' sox exec')
    startupinfo = None
    if( sys.platform == 'win32' ):
        startupinfo = subprocess.STARTUPINFO()
        startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
    info = subprocess.Popen([SOX_EXEC,'--i',filename], stdout=subprocess.PIPE, startupinfo=startupinfo).communicate()[0]
    reading_comments_flag = False
    for l in info.splitlines():
        if( not l.strip() ):
            continue
        if( reading_comments_flag and l.strip() ):
            if( comments ):
                comments += '\n'
            comments += l
        else:
            if( l.startswith('Input File') ):
                input_file = l.split(':',1)[1].strip()[1:-1]
            elif( l.startswith('Channels') ):
                num_channels = int(l.split(':',1)[1].strip())
            elif( l.startswith('Sample Rate') ):
                sample_rate = int(l.split(':',1)[1].strip())
            elif( l.startswith('Precision') ):
                bits_per_sample = int(l.split(':',1)[1].strip()[0:-4])
            elif( l.startswith('Duration') ):
                tmp = l.split(':',1)[1].strip()
                tmp = tmp.split('=',1)
                duration_time = tmp[0]
                duration_samples = int(tmp[1].split(None,1)[0])
            elif( l.startswith('Sample Encoding') ):
                encoding = l.split(':',1)[1].strip()
            elif( l.startswith('Comments') ):
                comments = ''
                reading_comments_flag = True
            else:
                try:
                    other += '\n'+l
                except NameError:
                    other = l
                if( output_unhandled ):
                    LOG.warn("Unhandled parameter: {}".format(l))
                pass
    # some file formats handled internally... (rest handled via sox)
    file_type = 'SOX'
    if( str(encoding) == '16-bit Signed Integer PCM' and filename[-4:].lower() == '.wav' ):
        file_type = sConstants.FILE_TYPE_WAV
    # return
    return (file_type, sample_rate, bits_per_sample/8, num_channels, duration_samples)

def GetAudioFileParams_native(filename, config):
    """returns (file_type, FS, bytes_per_sample, num_channels, num_samps)
        num_samps is aka 'frames', 1 samp includes all channels
        bytes_per_sample is for 1 sample from 1 channel, not a full samp/frame
        file_type is from global constants"""

    global constants

    (path,base) = os.path.split(filename)
    (base,ext) = os.path.splitext(base)

    if( ext == '.raw' ): # raw PCM data
        file_type = sConstants.FILE_TYPE_RAW
        # we need num_channels and bytes_per_sample from config
        num_channels = int(config.get('Recording Data', 'num_channels'))
        bytes_per_sample = int(config.get('Recording Data', 'bytes_per_sample'))
        FS = int(config.get('Recording Data', 'nominal_FS'))
        # get file length
        file_size = os.path.getsize(filename)
        num_samps = float(file_size) / float(bytes_per_sample) / float(num_channels)
        if( num_samps != int(num_samps) ):
            LOG.error("Problem reading file '{}': fractional number of samples? (bytes_per_sample, num_channels, or file type probably wrong)".format(filename))
            sys.exit(2)
        num_samps = int(num_samps)

    elif( ext == '.wav' ): # @TCC Need to verify (also good place to get channels and sample rate)
        file_type = sConstants.FILE_TYPE_WAV
        tmp = wave.open(filename,'r')
        num_channels = tmp.getnchannels()
        bytes_per_sample = tmp.getsampwidth()
        FS = tmp.getframerate()
        num_samps = tmp.getnframes()
        tmp.close()

    else:
        # try using sox: @TCC TODO need to enable and test/debug
        LOG.error("Unknown audio filetype for '{}'".format(filename))
        sys.exit(2)

    return (file_type, FS, bytes_per_sample, num_channels, num_samps)


### Recording session UTILITY FUNCTIONS #################################

### create (and return) the session ########################################
def CreateRecSession(config, msync_filename):
    if( msync_filename and msync_filename.lower() != 'none' ):
        (foo, ext) = os.path.splitext(msync_filename)
        if( ext.lower() == '.msync' ):
            return cVoxRecSession(config, msync_filename)
        else:
            return cWiredRecSession(config, msync_filename)


### cRecSession Class (virtual base class )########################################

class cRecSession(object) :
    def __init__(self, config, fn=None):
        raise RuntimeError("THIS SHOULD NOT BE INSTANSIATED")


### cWiredRecSession Class ########################################

class cWiredRecSession(cRecSession) :
    """This is for recording sessions where all nodes/mics are already synchronized recordings.
        Each node is a single channel continuous (not segmented) file."""

    def __init__(self, config):
        self.config = config
        self.nominal_FS = None
        self.num_samps = None
        self.filegroups = None
        self.session_identifier = None

    def init_from_metafile(self, metafile):
        LOG.info("Creating Rec Session from metafile '{}'".format(metafile))
        self.node_id_list = []
        self.metafile = metafile
        self.session_identifier = metafile
        self.nominal_FS = None
        self.num_samps = None
        # filegroup based organization....
        current_filegroup = []
        self.filegroups = OrderedDict()
        cur_node_id = 1
        base_path = os.path.dirname(metafile)
        with open(metafile, 'rU') as f:
            for l in f:
                l = l.strip()
                if( not l ):
                    continue
                if( l.startswith("annotationfile=") ):
                    self.annotationfile = l.split('=',1)[1]
                elif( l == 'filegroup:' ):
                    if( len(current_filegroup) > 0 ):
                        #self.filegroups.append(current_filegroup)
                        self.filegroups[str(cur_node_id).zfill(2)] = current_filegroup
                        cur_node_id += 1
                    current_filegroup = []
                else:
                    current_filegroup.append({'file':os.path.join(base_path,l)})
                    # @TCC maybe read params here too?
        # make sure the last filegroup is appended to the list
        if( len(current_filegroup) > 0 ):
            #self.filegroups.append(current_filegroup)
            self.filegroups[str(cur_node_id).zfill(2)] = current_filegroup
        self._init_file_params()

    def init_from_filename_pattern(self, filename_pattern, node_ids):
        LOG.info("Creating Rec Session from pattern '{}'".format(filename_pattern))
        self.filegroups = OrderedDict()
        self.filename_pattern = filename_pattern
        self.nominal_FS = None
        self.num_samps = None
        if( self.filename_pattern and self.filename_pattern.lower() != 'none' ):
            # determine the filenames from filename_pattern
            for nid in node_ids:
                self.filegroups[nid] = []
                fnp = self.filename_pattern.format(node_id=nid)
                tmp = fnp.find('*')
                if( tmp == -1 ):
                    raise RuntimeError("finename pattern must have '*' wildcard character")
                flist = sorted(glob.glob(fnp))
                for f in flist:
                    #ftmp = list(f)
                    #nid = ftmp[tmp:-(len(self.filename_pattern)-tmp-1)]
                    #nid = "".join(nid)
                    self.filegroups[nid].append({'file':f})
                LOG.debug("filegroup['{}'] {}".format(nid, self.filegroups[nid]))
            if( len(self.filegroups) == 0 ):
                raise RuntimeError("No data files found matching filename pattern '%s'"%self.filename_pattern)
        else:
            raise RuntimeError("finename pattern must be provided")
        self._init_file_params()

    def _init_file_params(self):
        # read and check the parameters for each file
        self.nominal_FS = None
        self.num_samps = None
        self.bytes_per_sample = None
        for k,v in self.filegroups.iteritems():
            for e in v:
                if( not os.path.isfile(e['file']) ):
                    raise IOError("File not found: '"+str(e['file']))
                (file_type, FS, bytes_per_sample, num_channels, num_samps) = \
                        GetAudioFileParams(e['file'], self.config)
                e.update({  'file_type':file_type, 'FS':FS,
                            'bytes_per_sample':bytes_per_sample,
                            'num_channels':num_channels,
                            'sample_rate':FS,
                            'num_samples':num_samps
                        })
                if( self.nominal_FS is None ):
                    self.nominal_FS = FS
                if( self.bytes_per_sample is None ):
                    self.bytes_per_sample = bytes_per_sample
                assert FS == self.nominal_FS, "Error: '"+str(e['file'])+"' All audio files must have same sample rate"
                assert bytes_per_sample == self.bytes_per_sample, "Error: '"+str(e['file'])+"' All audio files must have same bytes_per_sample"
                assert num_channels == 1, "Error: '"+str(e['file'])+"' Number of channels per file must be 1"
        self.AssertWiredFiles()

    def AssertWiredFiles(self):
        """Make sure the files (in self.filegroups) match across nodes.
        All nodes must have same # of files
        Number of samples in 1st file for each node must match, same for 2nd file, 3rd, ... N-1 (each recording segement must be same length, except last)
        """
        # make sure every node has the same number of files
        num_files = None
        for node_id in self.filegroups:
            if( num_files is None ):
                num_files = len(self.filegroups[node_id])
            elif( len(self.filegroups[node_id]) != num_files ):
                LOG.error("Node id {} has different number of sound files {}, expecting {}".format(node_id, len(self.filegroups[node_id]), num_files))
                sys.exit(2)
        file_lens = None
        # make sure all the soundfiles for each node (at the same segment) are the same lengths (except last files)
        for node_id in self.filegroups:
            fg = self.filegroups[node_id]
            if( file_lens is None ):
                file_lens = [ f['num_samples'] for f in fg ]
            else:
                for i,f in enumerate(fg):
                    if( i+1 < len(file_lens) and f['num_samples'] != file_lens[i] ): # last file being diffrent is ok
                        LOG.error("File lengths do not match across nodes. '{}' ({} samples), expecting {} samples".format(f['file'], f['num_samples'], file_lens[i]))
                        sys.exit(2)
        # @TCC should add more checks, but they are currently done on initialization (_init_file_params)


    def GetNodeIDs(self):
        return self.filegroups.keys()

    def GetNominalFS(self):
        return float(self.nominal_FS)

    def GetNumChannels(self, node_id):
        return 1 # must be one channel per node for wired

    def Samp2GTime(self, node_id, samp):
        """sample in recording session to global time"""
        return float(samp)/self.GetNominalFS()

    def GTime2Samp(self, node_id, gtime):
        """Global time to sample in recording session"""
        return int(NP.ceil(float(gtime)*self.GetNominalFS()))

    def FTime2GTime(self, node_id, soundfile, ftime):
        """Return the global time which corresponds to the time (ftime) in the specific file (soundfile)"""
        gtime = 0.0
        found = False
        for f in self.filegroups[node_id]:
            if( f['file'] == soundfile ):
                gtime += ftime
                found = True
                break
            else:
                gtime += f['num_samples'] / float(f['sample_rate'])
        if( not found ):
            LOG.error("Did not find file '{}' in recording session".format(soundfile))
            sys.exit(2)
        LOG.debug('node {}, file {}, ftime {} -> gtime {}'.format(node_id, soundfile, ftime, gtime))
        return gtime

    def ReadDataAtGTime(self, node_id, gtime, length_time, channels=None, center_flag=True, length_in_samps=False):
        """Read a block of audio data starting at gtime of length length_time (assuming nominal_FS)
             gtime is global timesync time
             length_time is assuming nominal sample rate, not measured sample rate (which varies)
                 this means that reads with the same length_time will return the same number of samples
             returns pylab/numpy array with a column per channel
             channels is list or single value (not a tuple) of channles to return
                    must be 1 (for now)
                    channels numbering starts with 1
                    channles == None returns all channels
             center_flag is ignored
             length_in_samps == True makes 'length_time' be interpreted as number of samples
            """
        assert( channels is None or channels == 1 or channels == [1] ), "Wired mode requires one channel audio files"
        if( length_in_samps ):
            length_samp = int(length_time)
            length_time = float(length_samp / self.GetNominalFS())
        else:
            length_samp = int(length_time * self.GetNominalFS() + .5) # rounding length up
            length_time = float(length_time)
        assert( length_samp > 0 )
        if( not node_id in self.filegroups ):
            LOG.error("Node ID '%s' not found"%node_id)
            sys.exit(2)
        start_samp = int(NP.round(gtime*self.nominal_FS))
        # if start is before 0, pre-pad with 0s
        if( start_samp < 0 ):
            LOG.warn("Requesting data before start of recording")
            if( hasattr(channels, "__len__") ):
                data = NP.zeros((-1*fsamp, len(channels)))
            else:
                data = NP.zeros(-1*fsamp)
            start_samp = 0
        file_idx = 0
        foffset = 0
        data = NP.array([])
        while( data.shape[0] < length_samp ):
            if( file_idx < len(self.filegroups[node_id]) ):
                e = self.filegroups[node_id][file_idx]
                fsamp = start_samp+data.shape[0] - foffset
                if( fsamp < e['num_samples'] ):
                    data = NP.concatenate((data, self._ReadDataFromFile(e, fsamp, length_samp-data.shape[0], channels=channels)))
                foffset += e['num_samples']
                file_idx += 1
            else:
                LOG.warn("Requesting data after end of recording")
                # post-pad with 0s
                if( hasattr(channels, "__len__") ):
                    data = NP.concatenate((data, NP.zeros((length_samp-data.shape[0], len(channels)))))
                else:
                    data = NP.concatenate((data, NP.zeros(length_samp-data.shape[0])))
        return data


    def _ReadDataFromFile(self, e, fsamp, length_samp, channels=None):
        # e is filegroups entry specifying a single file
        # fsamp is start sample number in file (inclusive)
        # lenght_samp is number of samples to read
        # @TCC FIXME Only tested on single channel files
        data = NP.array([])
        fn = e['file']
        file_type = e['file_type']
        bytes_per_sample = e['bytes_per_sample']
        file_len_samp = e['num_samples']
        Nc = e['num_channels']
        LOG.debug("ReadData: {} samples (offset={}, filelen={}) from {}".format(length_samp, fsamp, file_len_samp, fn))
        assert e['sample_rate'] == self.nominal_FS
        assert e['bytes_per_sample'] == self.bytes_per_sample
        if( file_type == sConstants.FILE_TYPE_RAW ): # RAW file type
            f = open(fn,'rb')
            f.seek(fsamp * Nc*bytes_per_sample) # initial offset
            data = NP.concatenate((data, NP.fromfile(f, NP.dtype('<i%d'%bytes_per_sample), min((length_samp*Nc-len(data), (file_len_samp-fsamp)*Nc)) ) ))
            f.close()
        elif( file_type == sConstants.FILE_TYPE_WAV ): # WAV file type
            f = wave.open(fn,'r')
            f.setpos(fsamp) # initial offset
            nframes_to_read = int(min((length_samp-(len(data)/Nc), file_len_samp-fsamp)))
            tmp = f.readframes(nframes_to_read)
            data = NP.concatenate((data, NP.fromstring(tmp,'<i%d'%bytes_per_sample) ))
            del tmp # don't need this now, and it could be large
            f.close()
        else: # other file type... use sox
            nframes_to_read = int(min((length_samp-(len(data)/Nc), file_len_samp-fsamp)))
            SOX_EXEC = self.config.get('System',sys.platform+' sox exec')
            out_byps = 2
            cmd = [SOX_EXEC,
                    fn,                    # input filename
                    '-t','raw',            # output file type raw
                    '-e','signed-integer', # output encode as signed ints
                    '-L',                  # output little endin
                    '-b',str(out_byps*8),  # output bytes per sample
                    '-',                   # output to stdout
                    'trim',str(int(fsamp))+'s',str(nframes_to_read)+'s'] # only extract requested part
            startupinfo = None
            if( sys.platform == 'win32' ):
                startupinfo = subprocess.STARTUPINFO()
                startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
            #print cmd
            data = NP.fromstring(subprocess.Popen(cmd, stdout=subprocess.PIPE, startupinfo=startupinfo).communicate()[0],'<i%d'%(out_byps))
            data = data.reshape(len(data)/Nc, Nc) # make samples x channels... don't need for one channel
        return NP.squeeze(data) # if only one channel, remove the signleton dimension



