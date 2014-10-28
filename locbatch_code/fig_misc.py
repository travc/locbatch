#!/usr/bin/python

__doc__ = """Localization figure related miscellaneous stuff
@TCC--Some more usage documentation goes here
"""

import sys
import math 
import ConfigParser

import numpy as NP
from scipy.signal import lfilter, butter, BadCoefficients # needed for filtering

# matplotlib, but only GUI/window-toolkit agnostic
import matplotlib
import matplotlib.transforms
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText

from misc_utils import list_unique, issequence
import mppool


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

