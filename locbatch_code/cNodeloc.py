#!/usr/bin/env python
"""The cNodeloc class and some constants for loading and dealing with self-survey data
    NOTE: nodeloc.mat files are assumed to number channels from 0, so we increment"""
import sys,os
import numpy as NP 
from numpy import sin, cos, pi
from scipy.io import loadmat 
from misc_utils import list_unique, str2bool
import loc_misc 
import ConfigParser

### Constants ###

## Sub-array geometry constants
# types = 'ENSBOXv1', 'VOX', or 'SINGLE'
# VoxNet (v2) sub-array geometry
d = .12; # meters
SUB_ARRAY_GEOMETRY = {}
SUB_ARRAY_GEOMETRY['VOX'] = [\
        [ +d/2.0, -d/2.0, -d/2.0, ],\
        [ -d/2.0, -d/2.0, +d/2.0, ],\
        [ -d/2.0, +d/2.0, -d/2.0, ],\
        [ +d/2.0, +d/2.0, +d/2.0, ] ]
# ENSBox (v1) sub-array geometry
d = .08; # meters
SUB_ARRAY_GEOMETRY['ENSBOXv1'] = [\
        [ +d/2.0, -d/2.0, -d/2.0, ],\
        [ +d/2.0, +d/2.0, +d/2.0, ],\
        [ -d/2.0, +d/2.0, -d/2.0, ],\
        [ -d/2.0, -d/2.0, +d/2.0, ] ]


### cNodeloc class ###

class cNodeloc:
    ids = []
    types = []
    pos = []
    sloc = []
    filename = None
    source_filename = None
    offset = None 
    force_zero_elevation = False
    force_zero_elevation_angle = False
    force_zero_roll_angle = False

    def __init__(self, filename, config, force_zero_elevation=False, force_zero_elevation_angle=False, force_zero_roll_angle=False):
        """@TCC to doc"""
        self.force_zero_elevation = str2bool(force_zero_elevation, False)
        self.force_zero_elevation_angle = str2bool(force_zero_elevation_angle, False)
        self.force_zero_roll_angle = str2bool(force_zero_roll_angle, False)
        mode = config.get('Main','mode')
        if( mode.lower() == 'wired' ):
            default_node_type = 'SINGLE'
        else:
            default_node_type = 'VOX'
        if( filename != None ):
            self.ReadSelfSurveyFile(filename, default_node_type)


    def Load(self, fn):
        self.ReadSelfSurveyFile(fn, default_node_type)


    def ReadSelfSurveyFile(self, fn, default_node_type):
        #print "ReadSelfSurveyFile:",fn
        base, ext = os.path.splitext(fn)
        if( ext.lower() == ".mat" ):
            # assume mat file files will be a nodeloc file
            nodes, pos, node_types, source_filename, offset, sloc = self.ReadNodelocFile(fn)
    
        elif( ext.lower() == ".loc" or ext.lower() == ".pl" ):
            # assume just a matrix of node locations... will need to compute sloc
            node_types = default_node_type 
            nodes, pos, node_types, source_filename, offset, sloc = self.ReadSelflocFile(fn, node_types)

        self.ids = nodes
        self.types = node_types
        self.pos = pos
        self.sloc = sloc
        self.filename = fn
        self.source_filename = source_filename
        self.offset = offset


    def ReadNodelocFile(self, fn):
        """ returns nodes, pos, node_types, source_filename, offset, sloc """
        ## read a nodeloc .mat file ##
        f = loadmat(fn, struct_as_record=True)
    
        if( not f.has_key('nodeloc') ):
            raise(RuntimeError,"'%s' does not appear to be a nodeloc file (does not define nodeloc)"%fn)
        if( not f.has_key('sloc') ):
            raise(RuntimeError,"'%s' does not appear to be a nodeloc file (does not define sloc)"%fn)
    
        # parse the nodeloc variable
        nl = f['nodeloc']
        varnames = nl.dtype.names
        if( True ):
            d = {}
            for v in varnames:
                d[v] = nl[v]
            
        # @TCC assumes node_ids are int numbers
        d['nodes'] = [ str(int(v)) for v in d['nodes'][0][0][0] ]
        #d['nodes'] = d['nodes'][0][0][0]
    
        d['pos'] = d['pos'][0][0]
        d['node_types'] = [ str(v[0][0]) for v in d['node_types'][0][0] ]
        if( not d['source_filename'][0][0] ):
            d['source_filename'] = None
        else:
            d['source_filename'] = d['source_filename'][0][0][0]
        if( not d['offset'][0][0] ):
            d['offset'] = None
        else:
            d['offset'] = d['offset'][0][0][0]

        if( self.force_zero_elevation ):
            for i in range(len(d['pos'])):
                d['pos'][i][2] = 0.

        # parse the sloc variable
        sloc_array = f['sloc'] # simple since it is just a matrix
        # convert sloc to python list
        sloc = []
        node_ids = []
        for ni in range(len(sloc_array)):
            sloc.append(list(sloc_array[ni]))
        for ni in range(len(sloc)):
            # assume node ids are ints (convert to str) @TCC!
            sloc[ni][0] = str(int(sloc[ni][0]))
            node_ids.append(sloc[ni][0])
            # channel numbers should start at 1
            sloc[ni][1] = int(sloc[ni][1])
        node_ids = list_unique(node_ids)

        if( True ): # recompute sloc (flags may differ from when nodeloc file was saved)
            sloc = self.Nodeloc2Sloc(node_ids, d['node_types'], d['pos'])

        if( sorted(d['nodes']) != sorted(node_ids) ):
            print >>sys.stderr, "WARNING: node IDs in nodeloc and sloc don't agree"
            print >>sys.stderr, "nodeloc nodes =", sorted(d['nodes'])
            print >>sys.stderr, "sloc nodes =   ", sorted(node_ids)

        return d['nodes'], \
                    d['pos'], \
                    d['node_types'], \
                    d['source_filename'], \
                    d['offset'], \
                    sloc


    def ReadSelflocFile(self, fn, node_types):
        """Read an older format self-loc / self-survey results file.
        returns node_ids, pos, node_types, source_filename, offset, sloc"""
        source_filename = fn
        offset = None
        node_ids = []
        pos = []
        infile = open(fn,'r')
        for l in infile:
            tmp = l.find("#")
            if( tmp == 0 ):
                continue # all of line is comment
            if( tmp > 0 ):
                l = l[0:tmp-1] # remove # and everything after
            l = l.strip()
            if( not l ):
                continue # empty line
            f = l.split()
            node_ids.append(f[0])
            pos.append([])
            pos[-1] = [ float(tmp) for tmp in f[1:] ]
    
        if( type(node_types) is str or not hasattr(node_types, '__iter__') ):
            node_types = [node_types]*len(node_ids)
        elif( len(node_types) == 1 ):
            node_types = node_types*len(node_ids)

        if( self.force_zero_elevation ):
            for i in range(len(pos)):
                pos[i][2] = 0.
    
        sloc = self.Nodeloc2Sloc(node_ids, node_types, pos)
        return node_ids, pos, node_types, source_filename, offset, sloc
    

    def Nodeloc2Sloc(self, node_ids, node_types, pos):
        """create sloc [[node, channel, x, y, z],...]
        NOTE: numbers channels from 1"""
        assert( len(node_ids) == len(node_types) )
        assert( len(node_ids) == len(pos) )
        sloc = []
        
        for ni in range(len(node_ids)):
    
            if( node_types[ni] == 'SINGLE' ):
                sloc.append([])
                sloc[-1].append(node_ids[ni])
                sloc[-1].append(1)
                sloc[-1].extend(pos[ni])
    
            elif( SUB_ARRAY_GEOMETRY.has_key(node_types[ni]) ):
                sl = NP.array(SUB_ARRAY_GEOMETRY[node_types[ni]])
                
                # rotation
                th_deg = pos[ni][3];
                if( len(pos[ni]) > 4 and not self.force_zero_elevation_angle ):
                    phi_deg = pos[ni][4]
                else:
                    phi_deg = 0.0
                if( len(pos[ni]) > 5 and not self.force_zero_roll_angle ):
                    roll_deg = pos[ni][5]
                else:
                    roll_deg = 0.0
    
                th = th_deg*pi/180;
                phi = phi_deg*pi/180;
                roll = roll_deg*pi/180;
    
                Rz = NP.array([ [cos(th), -sin(th), 0.], [sin(th), cos(th), 0.], [0., 0., 1.] ])
                Ry = NP.array([ [cos(phi), 0, sin(phi)], [0., 1., 0.], [-sin(phi), 0., cos(phi)] ])
                Rx = NP.array([ [1., 0., 0.], [0., cos(roll), -sin(roll)], [0., sin(roll), cos(roll)] ]) # roll
    
                sl = NP.dot(Rz, sl.T).T
                sl = NP.dot(Ry, sl.T).T
                sl = NP.dot(Rx, sl.T).T
    
                # translation
                for i in range(sl.shape[0]):
                    sl[i,:] = sl[i,:] + pos[ni][0:3]

                for i in range(sl.shape[0]):
                    sloc.append([])
                    sloc[-1].append(node_ids[ni])
                    sloc[-1].append(i+1)
                    sloc[-1].extend(list(sl[i,:]))

            else:
                raise(RuntimeError, "Sub-array type %s not in SUB_ARRAY_GEOMETRY constant",node_types[ni])

        return sloc


##### plotting using matplotlib #######

def Plot2DFigure(figure, nodeloc, title_str=None, show_z_values=True):
    """Actually plot the node locations in a matplotlib figure"""
    assert( len(nodeloc.pos) == len(nodeloc.ids) )
    pos = nodeloc.pos
    node_ids = nodeloc.ids

    ax = figure.add_axes([.1, .1, .85, .80])
    ax.hold(True)
    # node positions
    for ni in range(len(node_ids)):
        ax.plot([pos[ni][0]], [pos[ni][1]], 'bo')
    # node ids
    for ni in range(len(node_ids)):
        if( show_z_values ):
            ax.text(pos[ni][0], pos[ni][1], str(node_ids[ni]) + "\n(%.2f)"%(pos[ni][2]), \
                    horizontalalignment='center', \
                    verticalalignment='bottom', \
                    transform=loc_misc.MPL_offset_transform(ax,1,5))
        else:
            ax.text(pos[ni][0], pos[ni][1], str(node_ids[ni]), \
                    horizontalalignment='center', \
                    verticalalignment='bottom', \
                    transform=loc_misc.MPL_offset_transform(ax,1,5))
    # set the aspect to equal
    ax.set_aspect('equal', 'datalim') # make the x and y scales the same
    # node orientations (@TCC -- lots of variants of how to do this for future reference)
    xmin,xmax = ax.get_xlim()
    ymin,ymax = ax.get_ylim()
    dir_line_size = .015 * max([xmax-xmin, ymax-ymin])
    #tmpx1, tmpy1 = ax.transData.inverted().transform((0, 0))
    #tmpx2, tmpy2 = ax.transData.inverted().transform((5, 5))
    #dir_line_size = max([1/(tmpx2-tmpx1), 1/(tmpy2-tmpy1)]) * 2
    head_width = dir_line_size * 1
    alpha = .25
    for ni in range(len(node_ids)):
        if( len(pos[ni]) > 3 ):
            dx = dir_line_size*NP.cos((pos[ni][3])/180*NP.pi)
            dy = dir_line_size*NP.sin((pos[ni][3])/180*NP.pi)
            #z = dir_line_size*NP.sin((pos[ni][4])/180*NP.pi) # for 3D only
            #ax.plot([pos[ni][0], pos[ni][0]+dx], [pos[ni][1], pos[ni][1]+dy], 'r-') # just a line
            ax.arrow(pos[ni][0], pos[ni][1], dx, dy, head_width=head_width, alpha=alpha) # nice arrow
            #ax.annotate("", xy=(pos[ni][0], pos[ni][1]), xycoords='data',\
            #           xytext=(dx, dy), textcoords='offset points',\
            #           arrowprops=dict(arrowstyle="<|-") )
    # title, axis labels, grid, ect
    ax.tick_params(labelsize='x-small')
    if( title_str == None ): # default title (nodeloc filename)
        ax.set_title(nodeloc.filename, fontsize='small')
    elif( title_str == '' ): # no title
        pass
    else: # whatever title_str we passed in
        ax.set_title(title_str, fontsize='small')
    ax.set_xlabel('X [m]', fontsize='small')
    ax.set_ylabel('Y [m]', fontsize='small')
    ax.hold(False)
    ax.grid(True)



