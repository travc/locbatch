#!/usr/bin/python

__doc__ = """Correlation Envelope Sum Localization
    @TCC -- More usage documentation goes here"""

import sys
import configparser
import time
import math
import logging

import numpy as NP
from scipy.signal import decimate
import filtfilt

## local imports
import mppool
from mppool import mp_pool
import loc_misc
import gxcorr


USE_MP_FOR_DIST_CALC = False # MP is currently slower for the sensor to point distance calculation
PROGRESS_PRINTING = False

#### setup logger in module level scope ####
logging.basicConfig(level=logging.INFO)
LOG = logging.getLogger(__name__)


### cCESLoc ##################################################################

class cGCCParams():
    """parameters / options for computing cross-correlations for CES localization"""

    INTRA_NODE_LOC_ERROR = 1 # meters
    INTER_NODE_LOC_ERROR = 0.01 # meters
    SPEED_SOUND_ERROR_FACTOR = 0.01 # faction of speed of sound
    EXCLUDE_INTERNODE = False
    EXCLUDE_INTRANODE = False
    MAX_NUM_PASSES = 0
    USE_MP = True # false means don't use mp_pool even if it exists

    def __init__(self, config=None):
        if( config ):
            self.SetFromConfig(config)

    def SetFromConfig(self, config):
        """Set the parameters from a ConfigParser
             if config is string, loads that as a configuration filename"""
        # if config is string, assume it is a filename and load it
        if( isinstance(config, (str, bytes)) ):
            config_filename = config
            if( not os.path.isfile(config_filename) ):
                raise RuntimeError("String '{0}' passed to SetFromConfig, but file doesn't exit".format(config_filename))
            config = configparser.ConfigParser()
            if( config_filename ):
                config.read(config_filename)
        # set the values from config
        self.INTRA_NODE_LOC_ERROR = config.getfloat('CES','intra_node_loc_error')
        self.INTER_NODE_LOC_ERROR = config.getfloat('CES','inter_node_loc_error')
        self.SPEED_SOUND_ERROR_FACTOR = config.getfloat('CES','speed_of_sound_error_factor')
        self.EXCLUDE_INTERNODE = config.getboolean('CES','GCC_exclude_internode')
        self.EXCLUDE_INTRANODE = config.getboolean('CES','GCC_exclude_intranode')
        self.MAX_NUM_PASSES = config.getint('CES','GCC_max_num_passes')
        if( config.has_option('CES','GCC_use_multiprocessing') ):
            self.USE_MP = config.getboolean('CES','GCC_use_multiprocessing')



class cCESLoc():

    # member variables ('declaring' to be explicit)
    event = None
    session = None
    sloc = None
    gcc_params = None # cGCCParams instance (created from config in DoCES)
    speed_of_sound = None # computed from temperature and RH values of event in CompGCCs and reused in CompSum

    # results of CompGCCs
    gccs = None # generalized cross-correlation values (a dict for each gcc computed)
    # results   of CompSum
    max_C_pt = None # point with maximum pseudo-likelihood result (in meters)
    max_C_val = None
    max_C_pt_idxs = None # index of max point in C (needed for slicing views)
    C = None # pseudo-likelihood volume (Z,Y,X dimension order)
    C_extent = None # bounding box for C (in meters)
    grid_shape = None # dimensions (in points/index) of C
    x_vals = None # x points (in meters) where C was computed
    y_vals = None # y ""
    z_vals = None # z ""


    def __init__(self, config):
        self.config = config # needed for initilizing cRecSession.. @TCC should remove that dependency


    def CompGCCs(self, event, sloc, session, params=None, progress_callback=None): #lambda x,y: sys.stderr.write(x+'\n') or True):
        """Compute generalized cross-correlation pairs between sensors
            sloc is sensor/mic locations
            session is rec_session
            params is a cGCCParams, a ConfigParser, a config filename, or None (defaults)
            keep_going = progress_callback(message, fraction_done)
            allows aborting... returns False if aborted, True otherwise"""
        self.sloc = sloc
        self.session = session
        filter_func = loc_misc.CreateFilterFunc(self.session.GetNominalFS(), loc_misc.GenerateFilterDescription(event))
        if( isinstance(params, cGCCParams) ): # all good
            pass
        elif( not params ): # default values
            params = cGCCParams()
        else: # assume config or config_filename
            params = cGCCParams(params)

        # convience
        FS = self.session.GetNominalFS()
        #dur_samps = int(math.ceil((event.end_time - event.start_time)*FS))
        dur_samps = int(event['duration']*FS)

        # config options / values
        INTRA_NODE_LOC_ERROR = params.INTRA_NODE_LOC_ERROR
        INTER_NODE_LOC_ERROR = params.INTER_NODE_LOC_ERROR
        SPEED_SOUND_ERROR_FACTOR = params.SPEED_SOUND_ERROR_FACTOR
        EXCLUDE_INTERNODE = params.EXCLUDE_INTERNODE
        EXCLUDE_INTRANODE = params.EXCLUDE_INTRANODE
        MAX_NUM_PASSES = params.MAX_NUM_PASSES
        USE_MP = params.USE_MP

        if( not mp_pool() ): USE_MP = False

        tic = time.time() # timer

        ## compute speed of sound from temperature and RH
        self.speed_of_sound = loc_misc.CalcSpeedOfSound(event['temperature'], event['RH'])

        if( not USE_MP ):
            # make sure global variable is assinged in this process if not using workers
            globals()['mppool'].g_session = self.session
            globals()['mppool'].g_sloc = sloc
        else:
            # respawn mp_pool setting global variables for worker functions
            mp_pool().Respawn(None, {'g_session':self.session, 'g_sloc':sloc})

        ### compute the max dist+error (r_ij) for each pair of sensors
        r = [0]*len(sloc) # preallocate
        for i in range(len(sloc)):
            r[i] = [0]*len(sloc) # preallocate
            for j in range(i): # symmetric, only need lower tri
                r[i][j] = loc_misc.CalcDistance(sloc[i][2:5], sloc[j][2:5]); # distance between nodes
                if( sloc[i][1] == sloc[j][1] ): # same node?
                    r[i][j] += INTRA_NODE_LOC_ERROR
                else: # different nodes
                    r[i][j] += INTER_NODE_LOC_ERROR
                r[i][j] /= self.speed_of_sound  # convert from meters to sec
                # extra padding if smoothing cenv values to eliminate any possible edge effects (assumes window is in ms)
                if( event['smooth_cenv_window'] != None and event['smooth_cenv_window'] > 0 ):
                    r[i][j] += event['smooth_cenv_window']/1000.0/2.0
                # additional padding proportional to distance
                r[i][j] *= (1+SPEED_SOUND_ERROR_FACTOR)
                r[j][i] = r[i][j]  # fill in symmeteric value

        key_node = event['node_id']
        self.key_node_ids_used = [event['node_id']]
        key_channel = event['channel']

        # need to determine if start_time is global time or time in (segmented) file
        if( event['time_is_global'] ):
            start_gtime = event['start_time']
        else:
            start_gtime = self.session.FTime2GTime(key_node, event['soundfile'], event['start_time'])

        # bookkeeping variables
        cor_pairs_done = NP.zeros((len(sloc), len(sloc)))
        done = False
        pass_num = 0

        actual_max_num_passes = MAX_NUM_PASSES if MAX_NUM_PASSES > 0 else len(sloc)

        total_num_cors_to_compute = len(sloc)*(len(sloc)-1)/2
        if( EXCLUDE_INTERNODE or EXCLUDE_INTRANODE ):
            # actually count up the total num cors to compute
            total_num_cors_to_compute = 0
            for si1 in range(len(sloc)):
                for si2 in range(si1):
                    s1 = sloc[si1]
                    s2 = sloc[si2]
                    if( s1[0] == key_node and s1[1] == key_channel ): # key
                        total_num_cors_to_compute += 1 # always compute
                        continue
                    if( s2[0] == key_node and s2[1] == key_channel ): # key (symmetric)
                        total_num_cors_to_compute += 1 # always compute
                        continue
                    if( EXCLUDE_INTERNODE and s1[0] != s2[0] ):
                        continue  # skip
                    if( EXCLUDE_INTRANODE and s1[0] == s2[0] ):
                        continue  # skip
                    total_num_cors_to_compute += 1
        LOG.info("CompGCCs: num of cors to compute = {}".format(total_num_cors_to_compute))

        timing_key_node = key_node
        timing_key_channel = key_channel
        key_max_c_val = 1 # set for real after first pass
        key_offset = 0 # difference between first start_time and start_time of current key sloc

        cor = [] # actual results
        max_c_list = [] # list used determining which node to use as key in subsequent passes

        while( not done ):
            pass_num += 1

            if( PROGRESS_PRINTING ):
                print("\b"*40, file=sys.stderr)
                print("\bCompGCCs: pass {0}/{1}".format(int(pass_num),actual_max_num_passes), file=sys.stderr)

            if( progress_callback != None ):
                keep_going = progress_callback("Pass {0:d} of {1:d}".format(pass_num, actual_max_num_passes),
                                        len(cor)/float(total_num_cors_to_compute))
                if( not keep_going ): return False

            # load the key data
            gtime = start_gtime+key_offset
            length_samps = dur_samps
            key_data = self.session.ReadDataAtGTime(key_node, gtime, length_samps,
                    channels=[key_channel], center_flag=True, length_in_samps=True)
            if( filter_func ): # filter if the coefficents are defined
                key_data = filter_func(key_data)
            key_sloc_idx = loc_misc.FindInSloc(key_node, key_channel, sloc)

            # decimate @TCC TESTING
            decimate_factor = self.config.getint('TESTING','decimate_factor')
            if( decimate_factor > 1 ):
                key_data = decimate(key_data,decimate_factor)
                dFS = FS / float(decimate_factor)
            else:
                dFS = FS

            # create list of sloc idxs to gxcorr key with
            targets = []
            for i in range(len(sloc)):
                if( i == key_sloc_idx ): continue # skip self
                if( cor_pairs_done[i][key_sloc_idx] or cor_pairs_done[key_sloc_idx][i] ): continue # skip done already
                if( pass_num > 1 and EXCLUDE_INTRANODE and sloc[i][0] == sloc[key_sloc_idx][0] ): continue
                if( pass_num > 1 and EXCLUDE_INTERNODE and sloc[i][0] != sloc[key_sloc_idx][0] ): continue
                targets.append(i)
                cor_pairs_done[key_sloc_idx][i] = True # done really means "in target list or done"

            # compute the cross correlations

            if( USE_MP ): # using mp_pool
                ## dispach to worker processes
                #mp_pool().Respawn(None, {'g_session':self.session, 'g_sloc':sloc, 'g_key_data':key_data}) # @TCC could respawn and pass key_data that way instead, but doesn't seem any faster for event duration ~= 2s
                result_list = []
                for si in targets:
                    target_padding = math.ceil(r[key_sloc_idx][si]*FS) # padding to add to target in samples
                    result = mp_pool().apply_async( \
                            WorkerCompGXCorrPair, [key_sloc_idx, si, key_data, target_padding, gtime, dur_samps, filter_func, decimate_factor] )
                    if( result is None ):
                        LOG.error("RESULT IS NONE "+str([key_sloc_idx, si, key_data, target_padding, gtime, dur_samps, filter_func, decimate_factor]))
                        sys.exit(1)
                    else:
                        result_list.append(result)
                ## collect results
                for result in result_list:
                        cor.append(result.get())

            else: # single process
                for si in targets:
                    target_padding = math.ceil(r[key_sloc_idx][si]*FS) # padding to add to target in samples
                    c = WorkerCompGXCorrPair(key_sloc_idx, si, key_data, target_padding, gtime, dur_samps, filter_func, decimate_factor)
                    cor.append(c)

            if( progress_callback != None ):
                keep_going = progress_callback("Pass {0:d} of {1:d}".format(pass_num, actual_max_num_passes),
                                        len(cor)/float(total_num_cors_to_compute))
                if( not keep_going ): return False

            if( pass_num == 1 ): # first pass
                # create and sort the max_c_list (largest correlation value will be at end)
                # used to determine key sesnor for subsequent passes
                for c in cor:
                    max_c_list.append([c['target_sloc'][0:2], c['max_c_val'], c['max_c_lag']])
                max_c_list.sort(key=lambda x : x[1])

            # all done?
            if( key_max_c_val < event['gcc_key_max_c_val_threshold'] \
                    or len(max_c_list) < 1 \
                    or (MAX_NUM_PASSES and pass_num >= MAX_NUM_PASSES) ):
                done = 1
            else:
                # prep for next pass...
                # pick the sens with max correlation vs the first key as the next key
                tmp = max_c_list.pop()
                key_max_c_val = tmp[1]
                key_node = tmp[0][0]
                key_channel = tmp[0][1]
                key_offset = float(tmp[2])/FS
                self.key_node_ids_used.append(key_node)

        #print("\b"*40, file=sys.stderr)
        #print("total cors computed =",len(cor),"(including ones later culled)", file=sys.stderr)

        # done with while loop, so cleanup
        # cull out any values we had to compute which we don't actually want (first pass values)
        keep = []
        for c in cor:
            keep.append(True)
            if( EXCLUDE_INTRANODE and c['key_sloc'][0] == c['target_sloc'][0] ):
                keep[-1] = False
            if( EXCLUDE_INTERNODE and c['key_sloc'][0] != c['target_sloc'][0] ):
                keep[-1] = False
        cor = [cor[i] for i in range(len(cor)) if keep[i]]
        # cleanup globals
        if( not USE_MP ): # if assinged locally (not workers)
            del globals()['mppool'].g_session
            del globals()['mppool'].g_sloc

        # save results
        self.gccs = cor

        if( PROGRESS_PRINTING ):
            print("\b"*40, file=sys.stderr)
            print("\bCompGCCs (and envelopes): {0:d} gccs computed in {1:.2f} sec".format(len(cor), time.time()-tic), file=sys.stderr)

        LOG.info("CompGCCs (and envelopes): {0:d} gccs computed in {1:.2f} sec".format(len(cor), time.time()-tic))

        # Optional additional smoothing of envelopes
        # useful for high freq sounds and/or error in node positions
        # @TCC -- move this to the worker function ??
        if( event['smooth_cenv_window'] != None and event['smooth_cenv_window'] > 0 ):
            LOG.info("Additional cenv smoothing: convolve with {} ms hanning window".format(event['smooth_cenv_window']))
            for gcc in self.gccs:
                w = NP.hanning(max((1, NP.round(event['smooth_cenv_window']*1000.0/self.speed_of_sound))))
                gcc['cenv'] = NP.convolve(w/w.sum(), gcc['cenv'], mode='same')
        return True


    def CompSum(self, FUSE_AXIS, XY_RESOLUTION=1, Z_RESOLUTION=1, COMP_METHOD=0, USE_MP=True):
        """The sum/accumulation/fuse step"""

        if( not mp_pool() ): USE_MP = False

        assert( self.sloc != None ) # must have set sloc

        # hanlde if FUSE_AXIS is not a full 6 value arary
        if( len(FUSE_AXIS) == 6 ):
            pass # all good
        elif( len(FUSE_AXIS) == 4 ): # no z value provided, assume 0
            FUSE_AXIS.extend([0, 0])
        elif( len(FUSE_AXIS) == 5 ): # only has one z value...
            FUSE_AXIS.append(FUSE_AXIS[4])
        else:
            raise RuntimeError("FUSE_AXIS should be: minx, maxx, miny, maxy, minz, maxz values\ngot:{0}".format(FUSE_AXIS))

        # convience
        gccs = self.gccs
        session = self.session
        speed_of_sound = self.speed_of_sound
        FS = session.GetNominalFS()

        # decimated data? @TCC TESTING
        decimate_factor = self.config.getint('TESTING','decimate_factor')
        if( decimate_factor > 1 ):
            FS = FS / float(decimate_factor)

        sloc = self.sloc # possibly modified (culled) in earlier steps

        tic = time.time()
        LOG.info("CompSum: Sum extent = {}".format(FUSE_AXIS))

        # construct the grid of values to compute at [X,Y,Z]
        X = []
        Y = []
        Z = []
        x_extent = FUSE_AXIS[0], FUSE_AXIS[1]+XY_RESOLUTION
        y_extent = FUSE_AXIS[2], FUSE_AXIS[3]+XY_RESOLUTION
        x_list = NP.arange(x_extent[0], x_extent[1], XY_RESOLUTION)
        y_list = NP.arange(y_extent[0], y_extent[1], XY_RESOLUTION)

        # ensure z is only one slice if that is what we want
        if( Z_RESOLUTION <= 1e-8 or abs(FUSE_AXIS[4]-FUSE_AXIS[5]) < Z_RESOLUTION ):
            z_extent = FUSE_AXIS[4], FUSE_AXIS[4]
            z_list = [FUSE_AXIS[4]]
        else:
            z_extent = FUSE_AXIS[4], FUSE_AXIS[5]+Z_RESOLUTION
            z_list = NP.arange(z_extent[0], z_extent[1], Z_RESOLUTION)

        grid_shape = (len(z_list), len(y_list), len(x_list)) # NOTE: Z,Y,X order
        for z in z_list:
            for y in y_list:
                for x in x_list:
                    X.append(x)
                    Y.append(y)
                    Z.append(z)
        # make X,Y,Z numpy arrays (keep flattened for convience)
        X = NP.array(X)
        Y = NP.array(Y)
        Z = NP.array(Z)

        # compute distance from each sensor to points in the grid to
        # @TCC -- THIS IS SLOW
        #tic2 = time.time()
        if( not USE_MP_FOR_DIST_CALC or not USE_MP ):  # MP (at least this form) doesn't provide significant speedup
            sensor_dists = []
            for i in range(len(sloc)):
                pt = sloc[i][2:5]
                sensor_dists.append( NP.sqrt( (X-pt[0])**2 + (Y-pt[1])**2 + (Z-pt[2])**2 ) )
        else:
            # use multiprocessing to compute dists... need X,Y,Z, and sloc in globals
            mp_pool().Respawn(None, \
                    {'g_X':X,
                    'g_Y':Y,
                    'g_Z':Z,
                    'g_sloc':sloc,
                    })
            result_list = []
            num_batches = min(len(sloc), mp_pool().num_processes)
            for i in range(len(sloc)): #range(num_batches):
                #idxs = range(i,len(sloc), num_batches)
                idxs = [i]
                result = mp_pool().apply_async( WorkerCompDists, [idxs] )
                result_list.append(result)
            ## collect results
            sensor_dists = [None]*len(sloc)
            for result in result_list:
                (idxs, Ds) = result.get()
                for i in range(len(idxs)):
                    sensor_dists[idxs[i]] = Ds[i]
        #print "DIST TIME :",(time.time()-tic2)

        # respawn the mp_pool to pass global variables (or set in this process)
        if( not USE_MP ):
            # make sure global variable is assinged in this process if not using workers
            globals()['mppool'].g_session = self.session
            globals()['mppool'].g_sloc = sloc
            globals()['mppool'].g_gccs = gccs
            globals()['mppool'].g_sensor_dists = sensor_dists
        else:
            # respawn mp_pool setting global variables for worker functions
            mp_pool().Respawn(None, \
                    {'g_session':self.session,
                    'g_gccs':gccs,
                    'g_sloc':sloc,
                    'g_sensor_dists':sensor_dists,
                    })

        ## compute the actual sums

        if( not USE_MP ): # single process
            gccs_idxs = range(len(gccs))
            C = WorkerCompSum(gccs_idxs, FS, speed_of_sound, COMP_METHOD)

        else: # distribute the computation to workers
            result_list = []
            num_batches = min(len(gccs), mp_pool().num_processes)
            for i in range(num_batches):
                gccs_idxs = range(i,len(gccs), num_batches)
                result = mp_pool().apply_async( \
                        WorkerCompSum, [gccs_idxs, FS, speed_of_sound, COMP_METHOD] )
                result_list.append(result)
            ## collect results
            C = NP.zeros(len(sensor_dists[0])) # allocate space for sum results
            for result in result_list:
                C += result.get()

        # cleanup globals
        if( not USE_MP ): # if assinged locally (not workers)
            del globals()['mppool'].g_session
            del globals()['mppool'].g_sloc
            del globals()['mppool'].g_gccs
            del globals()['mppool'].g_sensor_dists
        else:
            pass
            # @TCC could kill/respawn the mp_pool to save on memory

        # find the location of the max C value
        mi = NP.argmax(C)
        max_C_val = C[mi]
        max_C_pt = ( X[mi], Y[mi], Z[mi] )
        # what is the indexs of max point... need for finding the max slice planes
        max_C_pt_idxs = NP.unravel_index(mi, grid_shape)

        # save values as instance variables
        self.max_C_pt = max_C_pt
        self.max_C_val = max_C_val
        self.max_C_pt_idxs = max_C_pt_idxs
        self.C = C.reshape(grid_shape)
        self.grid_shape = grid_shape
        #self.C_extent = NP.array([x_extent[0], x_extent[1], y_extent[0], y_extent[1]]) - XY_RESOLUTION/2.0
        self.C_extent = NP.array([x_list[0], x_list[-1]+XY_RESOLUTION, y_list[0], y_list[-1]+XY_RESOLUTION]) - XY_RESOLUTION/2.0 # @TCC -- this does not seem correct, but apparently works
        #self.C_extent = NP.array([X[0], X[-1], Y[0], Y[-1]])
        self.x_vals = x_list
        self.y_vals = y_list
        self.z_vals = z_list

        LOG.info("Grid size (x y z) = "+str(grid_shape[::-1]))
        LOG.info("Max location = {0[0]:.3f}, {0[1]:.3f}, {0[2]:.3f}  val={1:.3f}".format(max_C_pt, max_C_val))
        LOG.info("CompSum: done in {0:.2f} sec".format(time.time()-tic))


### global scope function ##########################################

def WorkerCompGXCorrPair(key_sloc_idx, target_sloc_idx, key_data, r_offset_samp, start_gtime, dur_samps, filter_func, decimate_factor=0):
    # actually compute cross-correlation for a pair of sensors
    # worker function (suitable for use with mp_pool worker processes)

    # Worker process namespace is a bit odd
    session = globals()['mppool'].g_session
    sloc = globals()['mppool'].g_sloc
    #key_data = globals()['mppool'].g_key_data # @TCC could pass key_through global (respawn workers) instead

    key_sloc = sloc[key_sloc_idx]
    target_sloc = sloc[target_sloc_idx]
    FS = session.GetNominalFS() # convience
    # load the 'target' data
    gtime = start_gtime - float(r_offset_samp/FS)
    length_samps = dur_samps + 2*r_offset_samp
    data = session.ReadDataAtGTime(target_sloc[0], gtime, length_samps, channels=[target_sloc[1]], \
                                                                    center_flag=True, length_in_samps=True)
    if( filter_func ): # filter if defined
        data = filter_func(data)

    # decimate @TCC TESTING
    if( decimate_factor > 1 ):
        data = decimate(data,decimate_factor)

    # cross correlation
    c, cenv, maxlag = gxcorr.GXCorr(key_data, data, r_offset_samp, mode='coeff-real', comp_cenv=True)
    # pack results
    cor = {}
    cor['c'] = c
    cor['maxlag'] = maxlag
    cor['cenv'] = cenv
    cor['key_sloc_idx'] = key_sloc_idx
    cor['target_sloc_idx'] = target_sloc_idx
    cor['key_sloc'] = key_sloc
    cor['target_sloc'] = target_sloc
    #cor['key_max_c_val'] = key_max_c_val  # not really needed
    # pack max correlation value and lag at max value
    tmp = NP.argmax(c)
    if( tmp.shape ): tmp = tmp[0] # take first if multiple
    cor['max_c_val'] = c[tmp]
    lags = range(-cor['maxlag'], cor['maxlag']+1)
    assert( len(lags) == len(cor['c']) )
    cor['max_c_lag'] = lags[tmp]
    # pack max envelope value and lag at max value
    if( cenv is None ): # it is possible we didn't compute envelopes
        tmp = NP.argmax(cenv)
        if( tmp.shape ): tmp = tmp[0] # take first if multiple
        cor['max_cenv_val'] = cenv[tmp]
        cor['max_cenv_lag'] = lags[tmp]
    else:
        cor['max_cenv_val'] = None
        cor['max_cenv_lag'] = None
    return cor


def WorkerCompSum(gccs_idxs, FS, speed_of_sound, COMP_METHOD=0):
    sloc = globals()['mppool'].g_sloc
    gccs = globals()['mppool'].g_gccs
    sensor_dists = globals()['mppool'].g_sensor_dists
    oobv = -10000 # value interp returns for extrapolation (something obvious since it shouldn't ever appear)

    C = NP.zeros(len(sensor_dists[0])) # allocate space for sum results
    for gccs_idx in gccs_idxs:
        # convience/clarity variables
        key_sloc_idx = gccs[gccs_idx]['key_sloc_idx']
        target_sloc_idx = gccs[gccs_idx]['target_sloc_idx']
        key_sloc = gccs[gccs_idx]['key_sloc']
        target_sloc = gccs[gccs_idx]['target_sloc']
        lags = range(-gccs[gccs_idx]['maxlag'], gccs[gccs_idx]['maxlag']+1)
        # expected TDOA between sensors in samples
        D = (sensor_dists[target_sloc_idx]-sensor_dists[key_sloc_idx]) * FS / speed_of_sound

        if( COMP_METHOD == 1 ): # c for all
            C += NP.interp(D, lags, gccs[gccs_idx]['c'], left=oobv, right=oobv)
        elif( COMP_METHOD == 2 ): # cenv for all
            C += NP.interp(D, lags, gccs[gccs_idx]['cenv'], left=oobv, right=oobv)

        else: # cenv for internode, c for intranode (default)
            if( key_sloc[0] == target_sloc[0] ): # same node, use correlation values (c)
                C += NP.interp(D, lags, gccs[gccs_idx]['c'], left=oobv, right=oobv)
            else: # internode, use envelope of correlation (cenv)
                C += NP.interp(D, lags, gccs[gccs_idx]['cenv'], left=oobv, right=oobv)
    return C


def WorkerCompDists(idxs):
    sloc = globals()['mppool'].g_sloc
    X = globals()['mppool'].g_X
    Y = globals()['mppool'].g_Y
    Z = globals()['mppool'].g_Z
    Ds = []
    for i in idxs:
        pt = sloc[i][2:5]
        Ds.append( NP.sqrt( (X-pt[0])**2 + (Y-pt[1])**2 + (Z-pt[2])**2 ) )
    return (idxs, Ds)
