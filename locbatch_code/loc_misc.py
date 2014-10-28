#!/usr/bin/python

__doc__ = """Localization related miscellaneous utility functions and classes
@TCC--Some more usage documentation goes here
"""

import sys
from sys import stderr,stdout
import math 
import ConfigParser
import warnings

import numpy as NP
from scipy.signal import lfilter, butter, BadCoefficients # needed for filtering

# matplotlib, but only GUI/window-toolkit agnostic
import matplotlib
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText

from misc_utils import list_unique, issequence
import mppool



### verbose helper funciton (would be macro in C) ###
VERBOSE_LEVEL_ = 8

def verbose(lvl,*mesg):
    global VERBOSE_LEVEL_
    if( lvl <= VERBOSE_LEVEL_ ): sys.stderr.write(' '.join([str(x) for x in mesg]) +"\n")


#### Math general #####################################

def CalcDistance(x,y): # N-dimension
    assert( len(x) == len(y) )
    d = 0
    for i in range(len(x)):
        d += (x[i]-y[i]) * (x[i]-y[i])
    return math.sqrt(d)


def Cart2Pol(x,y): # 2D 
    assert(len(x) == len(y))
    r = NP.sqrt(x*x + y*y)
    th = NP.arctan2(y,x)
    return th,r


#### Acoustics related general #####################################

def CalcSpeedOfSound(T, Rh=50, P=101.325):
    """Compute the speed of sound
        T: temperature in degees C
        Rh: relative humidity in % (default 50)
        P: pressure in kPa (default 101.325)
    adapted from:
    http://resource.npl.co.uk/acoustics/techguides/speedair/
    "This calculation shows the speed of sound in humid air according to Owen Cramer, "JASA, 93, p. 2510, 1993", with saturation
    vapor pressure taken from Richard S. Davis, "Metrologia, 29, p. 67, 1992", and a mole fraction of carbon dioxide of 0.0004.
    The calculator is valid over the temperature range 0 to 30' C (273.15 - 303.15 K) and over the pressure range 75 to 102 kPa.
    In the region between the air pressures 95.000 und 104.000 kPa there is no noticeable changing of the speed of sound c.
    The standard airpressure is 101325 Pa = 101.325 kPa or 1013.25 hectopascal."
    """
    # constants (for convience and clarity)
    Kelvin = 273.15 # For converting to Kelvin
    e = math.e
    pi = math.pi
    # ensure all the inputs are floats
    T_kel = float(T) + Kelvin  # ambient temperature in Kelvin
    Rh = float(Rh)
    P = float(P) * 1000.0  # convert pressure from kPa to Pa
    # Molecular concentration of water vapour (ENH) calculated from Rh
    # using Giacomos method by Davis (1991) as implemented in DTU report 11b-1997
    ENH = pi*10**(-8)*P + 1.00062 + T**2*5.6*10**(-7)
    PSV1 = T_kel**2*1.2378847*10**(-5)-1.9121316*10**(-2)*T_kel 
    PSV2 = 33.93711047-6.3431645*10**(3)/T_kel 
    PSV = e**PSV1*e**PSV2 
    H = Rh*ENH*PSV/P  # molecular concentration of water vapour
    Xw = H/100.0  # Mole fraction of water vapour 
    Xc = 400.0*10**(-6)  # Mole fraction of carbon dioxide 
    # Speed calculated using the method of Cramer from
    # JASA vol 93 pg 2510
    C1 = 0.603055*T + 331.5024 - T**2*5.28*10**(-4) + (0.1495874*T + 51.471935 -T**2*7.82*10**(-4))*Xw
    C2 = (-1.82*10**(-7)+3.73*10**(-8)*T-T**2*2.93*10**(-10))*P+(-85.20931-0.228525*T+T**2*5.91*10**(-5))*Xc
    C3 = Xw**2*2.835149 + P**2*2.15*10**(-13) - Xc**2*29.179762 - 4.86*10**(-4)*Xw*P*Xc
    C = C1 + C2 - C3  # speed
    return C


def CreateFilter(FS, filter_desc=None):
    # create a filter ([a,b] coefficients) from the filter_desc string, 
    # or take filter_desc as a filename and load values
    if( filter_desc == None ):
        return None
    filter_desc = filter_desc.strip()
    if( isinstance(filter_desc, basestring) ):
        if( len(filter_desc) > 7 and filter_desc[0:7] == '*butter' ):
            # create butterworth filter from description string
            l = filter_desc.split()
            order = int(l[1])
            kind = l[2]
            freq = l[3:]
            filter = CreateButterFilter(FS, order, kind, freq)
        else:
            # assume filter_desc is a filter filename (2 lines of whitespace separated numbers sepcifing A and B)
            print >>sys.stderr, "Loading filter from '{0}'".format(filter_desc)
            f = open(filter_desc,'r')
            a = [float(x) for x in f.readline().split()]
            b = [float(x) for x in f.readline().split()]
            filter = [a,b]
    else:
        raise RuntimeError, "Unparseable filter description '{0}'".format(filter_desc)
    return filter


def CreateButterFilter(FS, order, kind, freq):
    """create a butterworth filter
         FS: sample frequency
         order: order of filter to create 
         kind: 'low', 'high', or 'band'
         freq: single cutoff freq for kind 'low' or 'high', list of 2 frequencies for kind 'band'
            """
    #print "CreateButterFilter:",FS, order, kind, freq
    if( not issequence(freq) ): freq = [freq]
    freq = [ float(f) for f in freq ] # make sure freqs are floats
    valid_kinds = ['low','high','band']
    if( any(kind == t for t in valid_kinds) == False ):
        raise RuntimeError, "Unkown filter kind '{0}'".format(kind)
    if( kind == 'band' and len(freq) < 2 ):
        raise RuntimeError, "CreateButterFilter requires 2 frequencies for kind='band'"
    if( kind == 'band' and freq[0] >= freq[1] ):
        raise RuntimeError, "Min frequency must be < max freq for band filter"
    if( freq[-1] >= float(FS/2) ):
        raise RuntimeError, "Frequency ({0}) must be < sample_rate/2 ({1}) for filter".format(freq[-1], FS/2.)
    tmp = [ 2*f/FS for f in freq ] # recale freq to 0, 1 (whre 1 is nyquist)
    # try to create the filter, if we get a BadCoefficients warning, reduce the order and try again
    order_backoff = 1 # how much to reduce the order when backing off
    # @TCC -- it may make more sense to assume odd and backoff by 2
    while( True ):
        warnings.filterwarnings('error',"",BadCoefficients,"",0)
        try:
            (b,a) = butter(order, tmp, kind)
        except BadCoefficients:
            if( order <= 1 ):
                raise RuntimeError, "Error creating filter with order 1... very wrong"
            else:
                print >>sys.stderr, "WARNING: CreateButterFilter\n\tBadCoefficents for order",order,"... Backing off to order",order-order_backoff
                order = order-order_backoff
        else: # good filter
             break # done with while loop
    return [a,b]


def FindInSloc(node_id, channel, sloc, raise_on_fail=True):
    sloc_idx = None
    for i in range(len(sloc)):
        if( sloc[i][0] == node_id and sloc[i][1] == channel ):
            sloc_idx = i
            break
    if( raise_on_fail and sloc_idx == None ):
        raise RuntimeError, "Cannont find sloc entry for node_id '{0}' channel '{1}'".format(node_id, channel)
    return sloc_idx


def GetExtentFromSloc(sloc, xpadding=0.0, ypadding=None, zpadding=0.0, zplane=None):
    """compute an extent (bounding box) around sensor corrodinates in sloc,
            optionally adding padding
            if ypadding==None, ypadding=xpadding (default)
            zpadding is ignored if zplane is specified
            zplane: None = compute z extent normally from sloc (default)
                            'mean' = z extent values are mean of sensors
                            'median' = z extent values are median of sensors
                            {float} = z extent values are set to value """
    if( ypadding == None ): ypadding = xpadding
    x = []
    y = []
    z = []
    for si in range(len(sloc)):
        x.append(sloc[si][2])
        y.append(sloc[si][3])
        z.append(sloc[si][4])

    extent = [  min(x)-xpadding, max(x)+xpadding,
                            min(y)-ypadding, max(y)+ypadding,
                            min(z)-zpadding, max(z)+zpadding ]
    # handle specified zplane
    if( not zplane == None ):
        if( isinstance(zplane, basestring) ):
            if( zplane.lower() == 'mean' or zplane.lower() == 'average' or zplane.lower() == 'ave' ):
                extent[4] = NP.mean(z)
            elif( zplane.lower() == 'median' ):
                extent[4] = NP.median(z)
            else:
                raise RuntimeError, "Unknown zplane specifier '{0}'".format(zplane)
        else:
            extent[4] = zplane
        extent[5] = extent[4]
    # return
    return extent


def CullSloc(sloc, key_node_id, key_channel, sensors_to_use=None):
    """return a culled sloc only including entries for sensors listed in sensors_to_use
         key_node_id and key_channel are added if not in sensors_to_use
        sensors_to_use format is:
            [ [node_id1, channel1, channel2, ...], [node_id2, ...] ]
            if node_id is None, will include specificed channels for all nodes
            if channel is None (or not specified), will include all channels
        examples:
            sesnors_to_use=[['100',1], ['113',1,3]] : node 100 channel 1 + node 113 channels 1 and 3
            sesnors_to_use=[['100',]] : include all channels from node 100
            sensors_to_use=[[None,1,2]] : include channels 1 and 2 from all nodes"""
    stu = sensors_to_use # convience
    if( stu == None ): return sloc

    # get list of unique node_ids
    node_ids = list_unique([ s[0] for s in sloc ])

    # expand multiple node specification in sensors_to_use  
    tmp = []
    for s in stu:
        if( s[0] == None ):
            if( len(s) < 2 or s[1] == None ): # include all, so just return sloc
                return sloc
            for n in node_ids:
                tmp.append([n])
                tmp[-1].extend(s[1:])
        else:
            tmp.append(s)
    stu = tmp

    # expand multiple channel specification in sensors_to_use   
    tmp = []
    for s in stu:
        if( len(s) == 1 or s[1] == None ): # include all channels
            for s2 in sloc:
                if( s2[0] == s[0] ):
                    tmp.append([s2[0], s2[1]])
        else: # include specified channels
            for c in s[1:]:
                tmp.append([s[0],c])
    stu = tmp

    # ensure key node and channel are include
    if FindInSloc(key_node_id, key_channel, stu, False) == None:
        print >>sys.stderr, "WARNING: key/annotated node_id ('{0}') and channel ({1}) must be included in sensors to use, adding it in".format(key_node_id, key_channel)
        stu.append([key_node_id, key_channel])
    
    # actually cull sloc
    sloc_out = []
    for l in sloc:
        if( any([ s[0]==l[0] and s[1]==l[1] for s in stu ]) ):
            sloc_out.append(l)
    return sloc_out


def RoundFuseAxisToResolution(fuse_axis, XY_resolution, Z_resolution=None):
        # round fuse axis based on resolution
        fuse_axis[0] = math.floor(fuse_axis[0] / XY_resolution) * XY_resolution
        fuse_axis[1] = math.ceil(fuse_axis[1] / XY_resolution) * XY_resolution
        fuse_axis[2] = math.floor(fuse_axis[2] / XY_resolution) * XY_resolution
        fuse_axis[3] = math.ceil(fuse_axis[3] / XY_resolution) * XY_resolution
        if( len(fuse_axis)>=6 and fuse_axis[4] != fuse_axis[5] ):
            fuse_axis[4] = math.floor(fuse_axis[4] / Z_resolution) * Z_resolution
            fuse_axis[5] = math.ceil(fuse_axis[5] / Z_resolution) * Z_resolution
        return fuse_axis



###### plotting time aligned spectrograms #####

def PlotMultiSpecs(fig, session, gtime, duration, node_ids, channels, NFFT=512, overlap_fraction=None,
                    freq_range=None, view_freq_range=None,
                    key_node_id=None, normalize_together=True, filter_desc=None, padding=None):
    """Plot a spectrogram for each node listed in node_ids
         channels -- either a list (one for each node_id) or a single value 
         overlap_fraction -- None defaults to .25
         freq_range -- a (min,max) tuple, defaults to all (0,FS/2)"""

    # @TCC -- MAKE ARGS
    if( padding is None ):
        padding = max((0.25*duration, 1))
    interpolation = None # Plot smoothing; None==interp, 'nearest'==flat

    if not type(channels) is list:
        channels = [channels]*len(node_ids)
    assert( len(channels) == len(node_ids) )

    fig.clear()

    # define the subaxis grid
    if( True ):
        sub_ax = []
        ax = None
        for i in range(len(node_ids)):
            sub_ax.append(fig.add_subplot(len(node_ids), 1, i+1, sharex=ax, sharey=ax))
            if len(sub_ax) == 1:
                ax = sub_ax[0]
    fig.subplots_adjust(hspace=0, bottom=.025, top=.99, left=.1, right=.99)

    if( overlap_fraction == None ):
        overlap_fraction = .25
    noverlap = int(NFFT*overlap_fraction)

    # create the filter
    filter = None
    if( filter_desc ):
        filter = CreateFilter(session.GetNominalFS(), filter_desc)

    results = []
    if( mppool.mp_pool() ): # using multiprocessing pool
        # respawn mp_pool setting global variables for worker functions
        mppool.mp_pool().Respawn(None, {'g_session':session}) 
        result_obj_list = [mppool.mp_pool().apply_async(specgram_func, [i, gtime, duration, padding, node_ids[i], channels[i], NFFT, noverlap, freq_range, filter]) for i in range(len(node_ids))]
        for r in result_obj_list:
            results.append(r.get())

    else: # single process 
        globals()['mppool'].g_session = session # @TCC REMOVE IF RUNNING IN WORKERS
        for i in range(len(node_ids)):
            r = specgram_func(i, gtime, duration, padding, node_ids[i], channels[i], NFFT, noverlap, freq_range, filter)
            results.append(r)

    if normalize_together:
        (idx, Pxx, freqs, times) = results[0] # just to get intial min and max values
        min_val = Pxx.flatten()[0]
        max_val = Pxx.flatten()[0]
        for r in results:
            (idx, Pxx, freqs, times) = r
            m = min(Pxx.flatten())
            if( m < min_val ): min_val = m
            m = max(Pxx.flatten())
            if( m > max_val ): max_val = m
        norm = matplotlib.colors.Normalize(vmin=0, vmax=max_val, clip=False)
    else:
        norm = None # default normalization

    # actual spectrogram plots
    for r in results:
        (idx, Pxx, freqs, times) = r
        extent = times[0], times[-1], freqs[0], freqs[-1]
        ax = sub_ax[idx]
        ax.imshow(Pxx, cmap=matplotlib.cm.jet, extent=extent, origin='lower', norm=norm, interpolation=interpolation)
        ax.grid(False)
        ax.axis('tight')
        if( view_freq_range != None ):
            ax.set_ylim(view_freq_range)
        SetAxisFontSize(ax,'xx-small')
        #ax.set_ylabel('Frequency', fontsize='small')

        # annotate the node_id and channel (using nifty anchored text)
        if( key_node_id and key_node_id == node_ids[idx] ):
            color = 'r'
        else:
            color = 'w'
        at = AnchoredText(str(node_ids[idx])+"."+str(channels[idx]),
                    prop=dict(size=8), frameon=True, loc=2,)
        #at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
        at.patch.set_boxstyle("square,pad=0.")
        at.patch.set_alpha(0.5)
        at.patch.set_color(color)
        ax.add_artist(at)

        # draw lines showing actual timing edges
        #ax.axvline(gtime, color='w')
        if( key_node_id and key_node_id == node_ids[idx] ):
            color = 'r'
        else:
            color = 'w'
        t0 = times[0]-(NFFT/2/session.GetNominalFS())
        t_offset = (t0 - (gtime - padding))
        ax.axvline(gtime + t_offset, color=color, alpha=.5, linewidth=2)
        ax.axvline(gtime + duration + t_offset, color=color, alpha=.5, linewidth=2)


def specgram_func(i, gtime, duration, padding, node_id, channel, NFFT, noverlap, freq_range, filter=None):
    """Worker process to compute specgrams"""
    session = globals()['mppool'].g_session # Worker process namespace is a bit odd

    # actually read the audio data
    data = session.ReadDataAtGTime(node_id, gtime-padding, duration+2*padding, channel)
    if( filter ): # filter if the coefficents are defined
        data = lfilter(filter[1], filter[0], data)

    # determine some parameters from the data
    FS = session.nominal_FS
    if( data.ndim > 1 ):
        (total_num_samples, num_channels) = data.shape
    else:
        (total_num_samples,) = data.shape
        num_channels = 1
    if( num_channels > 1 ):
        raise RuntimeError("Can only compute specgram of one channel at a time")
    num_fft = (total_num_samples / noverlap ) - 2

    # spectrogram computation
    Pxx, freqs, times = matplotlib.mlab.specgram(data, NFFT=NFFT, Fs=FS, noverlap=noverlap, scale_by_freq=True)

    # cull down to the freq range (@TCC -- this seems horribly inefficient/overcomplicated)
    if( freq_range != None ):
        freqs_mask = NP.ones_like(freqs)
        freqs_mask[freqs < freq_range[0]] = 0
        freqs_mask[freqs > freq_range[1]] = 0
        freqs_mask = NP.nonzero(freqs_mask)[0]
        Pxx = Pxx[freqs_mask,:]
        freqs = freqs[freqs_mask]

    Pxx = 10. * NP.log10(Pxx) # convert to dB
    return (i, Pxx, freqs, times)


########### Matplotlib stuff ######

# Offsets in pixels for things specified in data coordinates
# New enough versions have offset_copy by Eric Firing:
import matplotlib
import matplotlib.transforms
if 'offset_copy' in dir(matplotlib.transforms):
    from matplotlib.transforms import offset_copy
    def MPL_offset_transform(ax, x, y):
        return offset_copy(ax.transData, x=x, y=y, units='dots')
else: # Without offset_copy we have to do some black transform magic
    from matplotlib.transforms import blend_xy_sep_transform, identity_transform
    def MPL_offset_transform(ax, x, y):
        # This trick makes a shallow copy of ax.transData (but fails for polar plots):
        trans = blend_xy_sep_transform(ax.transData, ax.transData)
        # Now we set the offset in pixels
        trans.set_offset((x,y), identity_transform())
        return trans


def SetAxisFontSize(ax, fontsize):
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)


class RectangleInteractor:
    """
    A create a rectangle with draggable corners 
    GetExtent() method returns data coordinates as (xmin, xmax, ymin, ymax)
    """

    def __init__(self, 
            ax, # the axes 
            extent, # (xmin,xmax, ymin,ymax) 
            release_callback=None, # called on button release as release_callback(self)
            epsilon=5, # max pixel distance to count as a vertex hit
            facecolor = 'none',
            filled = False,
            linestyle = '--',
            linecolor = 'r',
            linewidth = 1,
            linealpha = None,
            marker = 'o',
            markersize = 5,
            markeredgecolor = 'k',
            markerfacecolor = 'w',
            edgecolor = 'k',
            alpha = .75 ):

        if( linealpha == None ): linealpha = alpha

        self.ax = ax
        canvas = ax.figure.canvas
        self.epsilon = epsilon
        self.release_callback = release_callback

        xy = self._XYFromExtent(extent)

        # the ploygon
        self.poly = matplotlib.patches.Polygon(xy, animated=True, closed=True,
                fc=facecolor, ec=edgecolor, fill=filled, alpha=alpha )
        ax.add_patch(self.poly)

        x, y = zip(*self.poly.xy)
        self.line = matplotlib.lines.Line2D(x,y,
                linestyle=linestyle,
                color=linecolor,
                linewidth=linewidth,
                marker=marker, 
                markersize=markersize, 
                markeredgecolor=markeredgecolor, 
                markerfacecolor=markerfacecolor, 
                alpha=linealpha,
                animated=True)
        self.ax.add_line(self.line)

        cid = self.poly.add_callback(self.poly_changed) 
        self._ind = None # the active vert

        canvas.mpl_connect('draw_event', self.draw_callback)
        canvas.mpl_connect('button_press_event', self.button_press_callback)
        canvas.mpl_connect('button_release_event', self.button_release_callback)
        canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
        self.canvas = canvas


    def GetExtent(self):
        """returns sorted extent (xmin, xmax, ymin, ymax)"""
        e = self._GetExtent()
        if( e[0] > e[1] ):
            e[0], e[1] = e[1], e[0]
        if( e[2] > e[3] ):
            e[2], e[3] = e[3], e[2]
        return e


    def _GetExtent(self):
        x0 = self.poly.xy[0][0]
        x1 = self.poly.xy[2][0]
        y0 = self.poly.xy[0][1]
        y1 = self.poly.xy[1][1]
        return [x0, x1, y0, y1]

    def _XYFromExtent(self, (x0,x1,y0,y1)):
        xs = [x0, x0, x1, x1, x0]
        ys = [y0, y1, y1, y0, y0]
        return zip(xs,ys)

    def draw_callback(self, event):
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.ax.draw_artist(self.poly)
        self.ax.draw_artist(self.line)
        self.canvas.blit(self.ax.bbox)

    def poly_changed(self, poly):
        'this method is called whenever the polygon object is called'
        # only copy the artist props to the line (except visibility)
        vis = self.line.get_visible()
        matplotlib.artist.Artist.update_from(self.line, poly)
        self.line.set_visible(vis)  # don't use the poly visibility state

    def get_ind_under_point(self, event):
        'get the index of the vertex under point if within epsilon tolerance'
        # get display coords
        xy = NP.asarray(self.poly.xy)
        xyt = self.poly.get_transform().transform(xy)
        xt, yt = xyt[:, 0], xyt[:, 1]
        d = NP.sqrt((xt-event.x)**2 + (yt-event.y)**2)
        indseq = NP.nonzero(NP.equal(d, NP.amin(d)))[0]
        ind = indseq[0]
        if d[ind]>=self.epsilon:
            ind = None
        return ind

    def button_press_callback(self, event):
        'whenever a mouse button is pressed'
        if event.inaxes==None: return
        if event.button != 1: return
        self._ind = self.get_ind_under_point(event)

    def button_release_callback(self, event):
        'whenever a mouse button is released'
        if event.button != 1: return
        self._ind = None
        if( self.release_callback != None ):
            self.release_callback(self)

    def motion_notify_callback(self, event):
        'on mouse movement'
        if self._ind is None: return
        if event.inaxes is None: return
        if event.button != 1: return
        x,y = event.xdata, event.ydata

        x0, x1, y0, y1 = self._GetExtent()
        # indicies (initally) start from bottom left and go clockwise
        if( self._ind == 0 or self._ind == 1 ):
            x0 = x
        else:
            x1 = x
        if( self._ind == 0 or self._ind == 3 ):
            y0 = y
        else:
            y1 = y
        xy = self._XYFromExtent((x0,x1,y0,y1))

        self.poly.xy = xy #zip(xs,ys)
        self.line.set_data(zip(*self.poly.xy))

        self.canvas.restore_region(self.background)
        self.ax.draw_artist(self.poly)
        self.ax.draw_artist(self.line)
        self.canvas.blit(self.ax.bbox)

    def set_extent(self, extent):
        """programmatically move the selector"""
        self.extent = extent
        xy = self._XYFromExtent(self.extent)
        self.poly.xy = xy #zip(xs,ys)
        self.line.set_data(zip(*self.poly.xy))
        self.canvas.restore_region(self.background)
        self.ax.draw_artist(self.poly)
        self.ax.draw_artist(self.line)
        self.canvas.blit(self.ax.bbox)

