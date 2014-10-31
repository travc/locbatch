#!/usr/bin/python
"""localize events (aka. selections) from an annotation file

Options (all optional):
    --help, -h   : display this usage info
    --config, -c : the configuration file to use (required)
"""

import sys
import os
import time
import argparse
import glob
import ConfigParser
import subprocess
import errno
import csv
import logging

import numpy as NP
from scipy.signal import lfilter

# graphics
import matplotlib as mpl
#mpl.use('agg') # no X (so show won't work)
import matplotlib.patheffects
import matplotlib.pyplot as plt

# local
import mppool
from mppool import mp_pool
import misc_utils
import rec_session
import cEvent
import CES_localize
import loc_misc
import fig_misc # needed for PlotMultiSpecs and offset transform... TCC should move or at least cleanup

################################################################################

SCRIPT_PATH = os.path.split(os.path.realpath(__file__))[0]
CODE_PATH = SCRIPT_PATH
DEFAULT_CONFIG = os.path.join(SCRIPT_PATH,'default.cfg')

#### setup logger in module level scope ####
logging.basicConfig(level=logging.INFO)
LOG = logging.getLogger(__name__)

################################################################################
## Figure parameters (matplotlib settings)
mpl.rcdefaults()
params = {
    'figure.figsize': (10, 7.5),
    'figure.facecolor': 'white',
    #'image.cmap': 'YlGnBu', # set below with value from config

    'font.size': 9,
    'axes.titlesize': 9,
    'axes.labelsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,

    'xtick.major.size': 2.5,
    'ytick.major.size': 3,

    'legend.fontsize': 11,
    'legend.borderpad': 0.33,
    'legend.labelspacing': 0.33,
    'legend.handletextpad': .3,
    'legend.frameon': True,
    'legend.fancybox': True,

    'lines.markersize': 4.5,
    'lines.linewidth': 1,
    'lines.dash_capstyle': 'round',

    'axes.linewidth': 0.75,
    'patch.linewidth': 0.5,

    'mathtext.fallback_to_cm': False,
    'mathtext.fontset': 'custom',
    }
mpl.rcParams.update(params)

################################################################################

def CreateSLocFromPatNodeloc(nodeloc_filename):
    # create sloc (sensor locations) from a Patricelli style mic locations file
    sloc = []
    with open(nodeloc_filename, 'rU') as f:
        for l in f:
            l = l.strip().split('#',1)[0]
            if( l ):
                l = l.split()
                # node id (first column or index if only 2 columns)
                if( len(l) < 3 ): # if only x y values, node id is index
                    sloc.append([str(len(sloc)+1).zfill(2), 1])
                else:
                    sloc.append(l[0])
                    del l[0]
                # position
                if( len(l) < 3 ):
                    l.extend([0]*(3-len(l)))
                sloc[-1].extend([float(x) for x in l])
    return sloc

### Main ######################################################################
def exit(retcode=0):
    # cleanup and exit
    global g_start_tic
    print >>sys.stderr, "END {:.2f} secs elapsed".format(time.time()-g_start_tic)
    sys.exit(retcode)

def Main():
    print >>sys.stderr, "START"
    global g_start_tic
    g_start_tic = time.time()

    ### Config and Arguments Parsing ###
    # parse cfg_file argument
    conf_parser = argparse.ArgumentParser(description=__doc__)
    conf_parser.add_argument('-c', '--cfg-file', type=argparse.FileType('r'), required=True,
            help="Config file")
    args = conf_parser.parse_args(sys.argv[1:])

    # build the config (read config files)
    if args.cfg_file:
        cfg = ConfigParser.SafeConfigParser()
        cfg.optionxform=str # make the ConfigParser options case sensitive
        cfg.readfp(open(DEFAULT_CONFIG,'rU'))
        cfg.readfp(args.cfg_file)

    # set verbosity levels (in modules too)
    verbose_lvl = cfg.getint('Main','verbose')
    # @TCC TODO set log level in all the other modules
    tmp = [LOG, rec_session.LOG, cEvent.LOG]
    if( verbose_lvl == 0 ):
        map(lambda x: x.setLevel(logging.WARNING), tmp)
    elif( verbose_lvl == 1 ):
        map(lambda x: x.setLevel(logging.INFO), tmp)
    else:
        map(lambda x: x.setLevel(logging.DEBUG), tmp)
    del tmp

    ### ### ##########################################

    # more argument/config handling... filenames and paths
    data_path = cfg.get('Files', 'data_path')
    cfg_filename = args.cfg_file.name
    if( data_path is None or data_path == '' or data_path.strip().lower() == 'none' ):
        data_path = os.path.split(cfg_filename)[0]
    LOG.info("data path = '%s'",data_path)
    # set outdir
    outdir = cfg.get('Files', 'outdir')
    if( outdir is None or outdir == '' or outdir.lower() == 'none' ):
        # default outdir is the basename of the cfg_file (under the same directory as the cfg_file)
        tmp = os.path.split(os.path.abspath(cfg_filename))[0]
        outdir = os.path.splitext(os.path.basename(cfg_filename))[0]
        outdir = os.path.join(tmp, outdir+'_lbresults')
    LOG.info("outdir = '%s'",outdir)
    # create outdir if needed
    try:
        os.makedirs(outdir)
    except OSError as err:
        if( err.errno == errno.EEXIST and os.path.isdir(outdir) ):
            pass
        else:
            raise
    # output filenames (not including the directory part which is outdir)
    out_filename = cfg.get('Files', 'out_filename')
    if( out_filename is None or out_filename == '' or out_filename.lower() == 'none' ):
        if( outdir == '.' or outdir == '.\\' or outdir == './' ):
            out_filename = os.path.splitext(os.path.basename(cfg_filename))[0]+'_results_'+cfg.get('Files','events').strip().replace(' ','')+'.csv'
        else:
            out_filename = 'results_'+cfg.get('Files','events').strip().replace(' ','')+'.csv'

    # have to add to "Event defaults" section of config so it can be filled into each event object
    cfg.set('Event defaults','data_path', data_path)
    cfg.set('Event defaults','nodeloc_filename', os.path.realpath(os.path.join(data_path,cfg.get('Files','mics_filename'))))
#    cfg.set('Event defaults','soundfile', cfg.get('Files','soundfile'))
#    cfg.set('Event defaults','soundfile pattern', cfg.get('Files','soundfile pattern'))

    # matplotlib colormap
    mpl.rcParams['image.cmap'] = cfg.get('Fig','colormap')
    fig_highlight_color = cfg.get('Fig','highlight_color')

    ### ### ##########################################

    ## read microphone/node/sensor locations file ##
    base_sloc = CreateSLocFromPatNodeloc(os.path.join(data_path, cfg.get('Files','mics_filename')))

    ## Read the annotation file ##
    events_filename = cfg.get('Files', 'events_filename')
    events_mode = cfg.get('Files', 'mode')
    rs = None # rec_session will get loaded per event if needed
    if( events_mode.lower() == 'syrinx' ):
        events = cEvent.ReadSyrinxAnnotationsFile(cfg, os.path.join(data_path,events_filename))
    elif( events_mode.lower() == 'raven' ):
        events = cEvent.ReadRavenAnnotationsFile(cfg, os.path.join(data_path,events_filename))
    else:
        raise RuntimeError("Unknown annotation file mode '{}'".format(events_mode))
    assert len(events) > 0, "No events specified in annotations file '{}'".format(os.path.join(data_path,events_filename))

    ## parse config events (index) list (the list of event idxs to localize)
    event_idxs = misc_utils.range_string_to_indicies(cfg.get('Files','events'), len(events), inclusive=True, sort_unique=True, ones_based=True)
    LOG.info("Event indicies to localize = %s",str(event_idxs))
    # error check event_idxs and output informative error if needed
    for x in event_idxs:
        if( x < 1 ):
            LOG.error("Event indexes must be >0 : offending value = %s",str(x))
            exit(2)
        if( x > len(events) ):
            LOG.error("Event index '{}' > number of events ({}) in annotations file '{}'".format(x, len(events), os.path.join(data_path,events_filename)))
            exit(2)

    ## open the csv writer
    csvfile = open(os.path.join(outdir, out_filename),'wb')
    csvwriter = csv.writer(csvfile, dialect='excel')
    csvwriter.writerow(['idx', 'soundfile', 'node', 'start time', 'duration', 'low freq', 'high freq', 'X', 'Y', 'Z', 'C', 'err est', 'warning'])

    # Create Pool of worker processes
    mppool.CreateProcessPool(cfg.get('System','num_processes'))

    ## localization for each event in event_idxs
    for event_idx in event_idxs:
        e = events[event_idx-1]

        fig_out_filename_base = "event_{:04d}_{:010.4f}_{:s}".format(e.idx, e.start_time, e.node_id)
        LOG.info("-"*10+" Processing: "+str(fig_out_filename_base)+" "+"-"*10)
        LOG.debug(str(e))

        # create the rec_session if needed
        if( rs is None or os.path.join(data_path, e.rec_session_identifier) != rs.session_identifier ):
            rs = rec_session.cWiredRecSession(cfg)
            if( hasattr(e,'metafile') ):
                assert not hasattr(e,'soundfile_pattern'), "give a metafile OR a soundfile_pattern"
                rs.init_from_metafile(os.path.join(data_path, e.metafile))
            elif( hasattr(e,'soundfile_pattern') ):
                assert not hasattr(e,'metafile'), "give a metafile OR a soundfile_pattern"
                rs.init_from_filename_pattern(os.path.join(data_path, e.soundfile_pattern), node_ids=[x[0] for x in base_sloc])
            else:
                raise RuntimeError("Neither a metafile or a soundfile_pattern found")

        # optional output of multi-spectrogram
        if( cfg.getboolean('Main', 'multi_specs') ):
            freq_padding = 2000
            #e.filter_desc = '*butter 3 band 100 1000'
            fig = plt.figure(figsize=(8, 14*1.5))
            fig_misc.PlotMultiSpecs(fig, rs, e.start_time, e.duration, [str(x+1).zfill(2) for x in range(14)], 1,
                                    NFFT=1024, overlap_fraction=None,
                                    view_freq_range=[max((0, e.low_freq-freq_padding)), min((240000, e.high_freq+freq_padding))],
                                    #view_freq_range=[0,4000],
                                    key_node_id=e.node_id,
                                    normalize_together=True,
                                    filter_desc=e.filter_desc,
                                    padding=0.5)
            fig.suptitle("{0}\n{1}".format(e.name,e.filter_desc))
            fig.subplots_adjust(hspace=0, bottom=.025, top=.95, left=.05, right=.99)
            if( cfg.getboolean('Main','interactive') ):
                plt.show()
            else:
                LOG.warn("Saving multi-spec figure... This can be rather slow")
                plt.savefig(os.path.join(outdir,fig_out_filename_base+'_mspec.'+cfg.get('Fig','format')),
                            dpi=cfg.getint('Fig','dpi'))
            plt.close(fig)

        # setup figure
        fig = plt.figure()

        # Spectrogram of key event
        ax = fig.add_subplot(2,2,2)
        PlotKeySpectrogram(e, rs, fig, ax,
                padding=cfg.getfloat('Fig','spec_time_padding'),
                freq_padding=cfg.getfloat('Fig','spec_freq_padding'),
                NFFT=cfg.getint('Fig','spec_NFFT'),
                overlap_fraction=cfg.getfloat('Fig','spec_overlap_fraction'),
                vmin=cfg.getfloat('Fig','spec_vmin'),
                selection_color=fig_highlight_color)

        # optionally exclude some nodes by removing them from sloc
        sloc = list(base_sloc)
        if( hasattr(e,'exclude_nodes') and e.exclude_nodes != '' and e.exclude_nodes.lower() != 'none' ):
            tmp_list = [x.strip() for x in e.exclude_nodes.split(',')]
            for tmp in tmp_list:
                if( tmp == '' ):
                    continue
                if( tmp == e.node_id ):
                    LOG.warn("Cannot exclude key node '{}' (exclude_nodes={})".format(e.node_id, e.exclude_nodes))
                    continue
                try:
                    i = [x[0] for x in sloc].index(tmp)
                    del sloc[i]
                except ValueError:
                    LOG.warn("Node id '{}' not found in list of nodes (sloc)".format(tmp))

        # create CES loc object which does the real work
        cesloc = CES_localize.cCESLoc(cfg)

        ## compute cross correlations ##
        cesloc.CompGCCs(e, sloc, rs, params=cfg)

        ## sum step ##
        extent_X_padding = cfg.getfloat('CES','extent_X_padding')
        extent_Y_padding = cfg.getfloat('CES','extent_Y_padding')
        limits = loc_misc.GetExtentFromSloc(sloc, xpadding=extent_X_padding, ypadding=extent_Y_padding, zplane='mean')
        XY_RESOLUTION = cfg.getfloat('CES','sum_XY_resolution')
        Z_RESOLUTION = cfg.getfloat('CES','sum_Z_resolution')
        COMP_METHOD = cfg.getint('CES','sum_computation_method')
        # actually compute
        cesloc.CompSum(limits, XY_RESOLUTION=XY_RESOLUTION, Z_RESOLUTION=Z_RESOLUTION, \
                    COMP_METHOD=COMP_METHOD, USE_MP=True)
        # plot
        ax1 = fig.add_subplot(2,2,3)
        PlotPseudoLikelihoodSlice(cesloc, ax=ax1, show_title=True)
        PlotNodeLocations(ax1, sloc, key_node_list=cesloc.key_node_ids_used, key_color=fig_highlight_color)

        ## refine sum step ##
        # determine limits and resolution
        p = cesloc.max_C_pt # convience
        refine_XY_padding = cfg.getfloat('CES','refine_XY_padding')
        refine_Z_padding = cfg.getfloat('CES','refine_Z_padding')
        refine_XY_resolution = cfg.getfloat('CES','refine_XY_resolution')
        refine_Z_resolution = cfg.getfloat('CES','refine_Z_resolution')
        refine_limits = [
                        p[0]-refine_XY_padding,
                        p[0]+refine_XY_padding,
                        p[1]-refine_XY_padding,
                        p[1]+refine_XY_padding,
                        p[2]-refine_Z_padding,
                        p[2]+refine_Z_padding ]
        # make the refine limits fall on resolution increments
        refine_limits = loc_misc.RoundFuseAxisToResolution(refine_limits, refine_XY_resolution, refine_Z_resolution)
        # actually compute
        cesloc.CompSum(refine_limits, XY_RESOLUTION=refine_XY_resolution, Z_RESOLUTION=refine_Z_resolution, \
                    COMP_METHOD=COMP_METHOD, USE_MP=True)
        # plot
        ax2 = fig.add_subplot(2,2,4)
        PlotPseudoLikelihoodSlice(cesloc, ax=ax2, show_title=True)
        PlotNodeLocations(ax2, sloc, key_node_list=cesloc.key_node_ids_used, key_color=fig_highlight_color)
        # crosshairs to mark max point
        ax2.axvline(x=cesloc.max_C_pt[0], color=fig_highlight_color, lw=0.25, alpha=0.75)
        ax2.axhline(y=cesloc.max_C_pt[1], color=fig_highlight_color, lw=0.25, alpha=0.75)
        # annotate wide plot too
        ax1.axvline(x=cesloc.max_C_pt[0], color=fig_highlight_color, lw=0.25, alpha=0.75)
        ax1.axhline(y=cesloc.max_C_pt[1], color=fig_highlight_color, lw=0.25, alpha=0.75)

        # rough max error estimate
        err_thresh_factor = cfg.getfloat('CES', 'error_est_threshold')
        warning_string = ''
        if( err_thresh_factor > 0 ):
            err_est_thresh = err_thresh_factor * cesloc.max_C_val
            cidx = NP.column_stack(NP.where(cesloc.C > err_est_thresh))
            pt = cesloc.max_C_pt
            err_est_d = max( NP.sqrt((cesloc.x_vals[cidx[:,2]]-pt[0])**2 +
                                (cesloc.y_vals[cidx[:,1]]-pt[1])**2 +
                                (NP.array(cesloc.z_vals)[cidx[:,0]]-pt[2])**2 ) )
            err_est_d += NP.sqrt(refine_XY_resolution**2)
            if( pt[0]-err_est_d < cesloc.C_extent[0] or
                pt[0]+err_est_d > cesloc.C_extent[1] or
                pt[1]-err_est_d < cesloc.C_extent[2] or
                pt[1]+err_est_d > cesloc.C_extent[3] ):
                LOG.warn("Max point near edge of sum extent")
                warning_string += "Max point near edge of sum extent"
            # plot
            #ax2.scatter(cesloc.x_vals[cidx[:,2]], cesloc.y_vals[cidx[:,1]], marker='.')
            err_circ = mpl.patches.Circle((pt[0], pt[1]), radius=err_est_d,
                        facecolor='none', edgecolor=fig_highlight_color, alpha=0.75,
                        linewidth=0.25,
                        #path_effects=[matplotlib.patheffects.withStroke(linewidth=1, foreground="w", alpha=0.5)],
                        clip_on=False)
            ax2.add_artist(err_circ)
            # annotate wide plot too
            err_circ = mpl.patches.Circle((pt[0], pt[1]), radius=err_est_d,
                        facecolor='none', edgecolor=fig_highlight_color, alpha=0.75,
                        linewidth=0.25,
                        #path_effects=[matplotlib.patheffects.withStroke(linewidth=1, foreground="w", alpha=0.5)],
                        clip_on=False)
            ax1.add_artist(err_circ)
        else:
            err_est_d = None
            err_est_thresh = None

        # text (use unicode for special symbols)
        relpath_start = os.getcwd()
        if( data_path is not None and data_path != '' ):
            relpath_start = os.path.realpath(data_path)
        s = u""
        s += u"name: {0.name:s}\n".format(e)
        s += u"idx: {:d}\n".format(e.idx)
        s += u"key node: {0.node_id:s}\n".format(e)
        s += u"start: {0.start_time:.3f} s\n".format(e)
        s += u"duration: {0.duration:.3f} s\n".format(e)
        s += u"low  freq: {0:g} hz\nhigh freq: {1:g} hz\n".format(e.low_freq, e.high_freq)
        s += u"filter description: {0:s}\n".format(e.filter_desc)
        if( e.smooth_cenv_filter is not None and 
                e.smooth_cenv_filter != '' and
                e.smooth_cenv_filter != 'none' ):
            s += u"additional cenv smoothing filter: {}\n".format(e.smooth_cenv_filter)
        s += u"temperature / RH: {0.temperature:g}{1:s}C / {0.RH:g}%\n".format(e, u'\N{DEGREE SIGN}')
        s += u"key C val thresh: {0:g}\n".format(cfg.getfloat('CES','GCC_key_max_c_val_threshold'))
        s += u"all key nodes: {0:s}\n".format(",".join(sorted(cesloc.key_node_ids_used)))
        s += u"data dir: {0:s}\n".format(data_path)
        if( hasattr(e,'soundfile') ):
            s += u"soundfile: {0:s}\n".format(os.path.relpath(e.soundfile, relpath_start))
        if( hasattr(e,'metafile') ):
            s += u"metafile: {0:s}\n".format(os.path.relpath(e.metafile, relpath_start))
        if( hasattr(e,'soundfile_pattern') ):
            s += u"soundfile_pattern: {0:s}\n".format(e.soundfile_pattern)
        s += u"mics file: {0:s}\n".format(os.path.relpath(e.nodeloc_filename, relpath_start))
        s += u"\nMax point: [{0[0]:.3f}, {0[1]:.3f}, {0[2]:.3f}] = {1:.3f}\n\n".format(cesloc.max_C_pt, cesloc.max_C_val)
        s += u"Error indicator: {:.2f} m (thresh={:.3f}, c={:.3f})\n".format(err_est_d, err_thresh_factor, err_est_thresh)
        if( warning_string ):
            s += u'\nWARNING: '+warning_string+"\n"
        LOG.info("\n\t"+"\n\t".join(s.splitlines()))
        ax3 = fig.add_subplot(2,2,1)
        ax3.axis('off')
        tmp = ax3.text(0,1, s, va='top', ha='left', multialignment='left', fontsize='8', clip_on=False)

        # ave or show the main figure
        if( cfg.getboolean('Main','interactive') ):
                plt.show()
        else:
            plt.tight_layout()
            plt.savefig(os.path.join(outdir,fig_out_filename_base+'_loc.'+cfg.get('Fig','format')),
                        dpi=cfg.getint('Fig','dpi'))

        # anything else which needs to be prepared/set for csv output
        if( hasattr(e,'soundfile') ):
            sndfile = e.soundfile
        else:
            sndfile = e.metafile # Events should have a soundfile or a metafile set

        # output the CSV row
        outrow = [ e.idx,
                sndfile,
                e.node_id,
                e.start_time,
                e.duration,
                e.low_freq,
                e.high_freq,
                ]
        outrow.extend(cesloc.max_C_pt)
        outrow.append(cesloc.max_C_val)
        outrow.append(err_est_d)
        outrow.append(warning_string)
        csvwriter.writerow(outrow)

    # end for each event to localize
    csvfile.close()

    exit()

###################################################################################################################

def PlotNodeLocations(ax, sloc, plot_sym='wo', text_color='k', key_color='r', alpha=.65, text_alpha=.75, preserve_limits=True, key_node_list=None):
    """add plot of node locations to axis provided, preserves existing limits
      assumes 1 channel per node
    """
    if( text_alpha is None ):
        text_alpha = alpha
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.hold(True)
    # id
    for ni in range(len(sloc)):
        if( key_node_list is not None and
            sloc[ni][0] == key_node_list[0] ):
            tmp_color = key_color
            tmp_alpha = 1
        elif( key_node_list is not None and
              sloc[ni][0] in key_node_list ):
            tmp_color = text_color
            tmp_alpha = text_alpha
        else:
            tmp_color = '#666666'
            tmp_alpha = text_alpha*.5
        tmp_sym = plot_sym
        path_alpha = tmp_alpha*.75

        ax.plot([sloc[ni][2]], [sloc[ni][3]], tmp_sym, alpha=alpha)
        ax.text(sloc[ni][2], sloc[ni][3], str(sloc[ni][0]), \
        color=tmp_color,
        alpha=tmp_alpha,
        horizontalalignment='center', \
        verticalalignment='bottom', \
        transform=fig_misc.MPL_offset_transform(ax,0,2.5),\
        clip_on=True,
        path_effects=[matplotlib.patheffects.withStroke(linewidth=1.25, foreground="w", alpha=path_alpha)],
        fontsize=8,
        )
    # restore limits if asked for
    if( preserve_limits ):
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

def PlotKeySpectrogram(e, rs, fig, ax, padding, freq_padding, NFFT, overlap_fraction, vmin, selection_color='r'):
    # e is an event object
    # rs is rec_session
    freq_low = e.low_freq
    freq_high = e.high_freq
    noverlap = int(NFFT*overlap_fraction)
    view_freq_range=[max((0, freq_low-freq_padding)),
                     min((240000, freq_high+freq_padding))]
    FS = rs.nominal_FS
    filter = loc_misc.CreateFilter(FS, e.filter_desc)
    if( e.time_is_global ):
        start_time = e.start_time
    else:
        start_time = rs.FTime2GTime(e.node_id, e.soundfile, e.start_time)
    data = rs.ReadDataAtGTime(e.node_id, start_time-padding, e.duration+2*padding, 1) # @TCC 1 channel only
    data = lfilter(filter[1], filter[0], data)
    # spectrogram computation
    Pxx, freqs, times = matplotlib.mlab.specgram(data, NFFT=NFFT, Fs=FS, noverlap=noverlap, scale_by_freq=True)
    Pxx = 10. * NP.log10(Pxx) # convert to dB
    # plot
    times = [x-padding+start_time for x in times]
    extent = times[0], times[-1], freqs[0], freqs[-1]
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=NP.amax(Pxx), clip=False)
    cax = ax.imshow(Pxx,
                    extent=extent, origin='lower',
                    norm=norm,
                    interpolation='nearest')
    cbar = fig.colorbar(cax, drawedges=False)
    cbar.solids.set_edgecolor("face") # workaround common problem with pdf viewers
    ax.grid(False)
    ax.axis('tight')
    if( view_freq_range != None ):
        ax.set_ylim(view_freq_range)
    # draw a rectangle around the time+freq range of the event
    rect = mpl.patches.Rectangle( (start_time, freq_low),
            e.duration, freq_high-freq_low,
            edgecolor='w', facecolor='none', linewidth=1, alpha=0.5)
    ax.add_artist(rect)
    rect = mpl.patches.Rectangle( (start_time, freq_low),
            e.duration, freq_high-freq_low,
            edgecolor=selection_color, facecolor='none', linewidth=.25, alpha=1)
    ax.add_artist(rect)

def PlotPseudoLikelihoodSlice(cesloc, ax=None, z_idx=None, show_title=True):
    if( z_idx == None ):
        z_idx = cesloc.max_C_pt_idxs[0]
    if( ax is None ):
        fig = plt.figure()
        ax = plt.subplot(1,1,1)
    im = ax.imshow(cesloc.C[z_idx,:,:], extent=cesloc.C_extent[0:4], origin='lower', interpolation='nearest')
    ax.set_aspect('equal', 'datalim') # make the x and y scales the same
    ax.grid(True)
    cbar = plt.colorbar(im)
    cbar.solids.set_edgecolor("face") # workaround common problem with pdf viewers
    if( show_title ):
        title_str = "Z = {0:.3f}".format(cesloc.z_vals[z_idx])
        #title_str = "[{0[0]:.3f}, {0[1]:.3f}, {0[2]:.3f}] = {1:.3f}".format(cesloc.max_C_pt, cesloc.max_C_val)
        ax.set_title(title_str)
        ax.set_xlabel('X [m]')
        ax.set_ylabel('Y [m]')


#########################################################################
# Main loop hook... if run as script run main, else this is just a module
if __name__ == "__main__":
    sys.exit(Main())
