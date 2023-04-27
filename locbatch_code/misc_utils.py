#!/usr/bin/python

__doc__ = """Miscellaneous utility functions and classes
@TCC--Some more usage documentation goes here
"""

import os,sys
from sys import stderr,stdout
import math
import configparser

### verbose helper funciton (would be macro in C) ###
VERBOSE_LEVEL_ = 8

def verbose(lvl,*mesg):
    global VERBOSE_LEVEL_
    if( lvl <= VERBOSE_LEVEL_ ): sys.stderr.write(' '.join([str(x) for x in mesg]) +"\n")


######### Programming/Python helpers ##################################################

def list_unique(hasDupes):
    """Return the sorted unique values from a list"""
    # order preserving
    from operator import itemgetter
    d = dict((x, i) for i, x in enumerate(hasDupes))
    return [k for k, _ in sorted(d.items(), key=itemgetter(1))]


def issequence(obj):
    """returns True if obj is non-string iterable (list, tuple, ect)"""
    return getattr(obj, '__iter__', False)


def StrOrNone(str):
    if( isinstance(str, basestring) and str.lower() == 'none' ):
        str = None
    return str

def CastOrNone(val, type_func):
    if( val == None ):
        return None
    # if val is a string, first check for "none"
    if( isinstance(val, basestring) ):
        if( not val or val.lower() == 'none' ):
            return None
    try:
        val = type_func(val)
    except ValueError:
        val = None
    return val


def str2bool(str, default_value=True):
    if( isinstance(str, basestring) ):
        if( default_value == True ):
            rv = True
            if( str.lower() == 'false' ):
                rv = False
        else: # default_value == False
            rv = False
            if( str.lower() == 'true' ):
                rv = True
    else:
        rv = str # assume what we have passed in is actually a bool to start with
    return rv


def minmax(a):
    """returns min_val, min_idx, max_val, max_idx for iterable a
         idexes will be first occurance"""
    min_val = a[0]
    min_idx = 0
    max_val = min_val
    max_idx = 0
    for i in range(len(a)):
        if( a[i] < min_val ):
            min_val = a[i]
            min_idx = i
        if( a[i] > max_val ):
            max_val = a[i]
            max_idx = i
    return min_val, min_idx, max_val, max_idx

def minmax_vals(a):
    """returns min_val, max_valfor iterable a"""
    min_val = a[0]
    max_val = min_val
    for i in range(len(a)):
        if( a[i] < min_val ):
            min_val = a[i]
        elif( a[i] > max_val ):
            max_val = a[i]
    return min_val, max_val


def ConfigGet(config, section, option):
    """Get a value (string) from a ConfigParser or return None if it isn't there
       the string 'none' (case insensitive) will return None too"""
    try:
        val = config.get(section, option)
        # 'none' should be translated to None
        if( val.lower() == 'none' ): val = None
    except (configparser.NoOptionError, configparser.NoSectionError):
        val = None
    return val

def ConfigGetBoolean(config, section, option):
    try:
        val = config.getboolean(section, option)
    except (configparser.NoOptionError, configparser.NoSectionError):
        val = None
    return val

def ConfigGetInt(config, section, option):
    try:
        val = config.getint(section, option)
    except (configparser.NoOptionError, configparser.NoSectionError):
        val = None
    return val

def ConfigGetFloat(config, section, option):
    try:
        val = config.getfloat(section, option)
    except (configparser.NoOptionError, configparser.NoSectionError):
        val = None
    return val



class TransplantFunc2Method:
    """Creates a class method in host (class) out of a function"""
    def __init__(self, method, host, method_name=None):
        self.host = host
        self.method = method
        setattr(host, method_name or method.__name__, self)
    def __call__(self, *args, **kwargs):
        nargs = [self.host]
        nargs.extend(args)
        return self.method(*nargs, **kwargs)



def ShortenPath(path):
    """takes path, converts it to absolute, and trims off cwd if that is the first part
       if path passed in is shorter, returns that instead"""
    in_path = path
    path = os.path.abspath(path)
    cwd = os.getcwd()+os.sep
    if( path[:len(cwd)] == cwd ):
        path = path[len(cwd):]
    if( len(in_path) < len(path) ): # if the input is shorter, just return it instead
        path = in_path
    return path


def range_string_to_indicies(string, list_len, inclusive=True, sort_unique=True, ones_based=False, zero_pad=0):
    """ converts a string like "-2,4,10-13,20-" into a list of indicies [1,2,10,11,12,13,20,21,22...]
        if string is None or 'all', will return the complete list of indicies (0 to list_len-1)
        ones_based=True makes it 1-based (start with 1 end with list_len)
        ranges (like "10-13") are inclusive of the end point if inclusive=True
        NOTE: doesn't check that values are actually in the proper range
    """
    if( string is None or string.lower() == 'all' ):
        idxs = range(list_len)
        if( ones_based ):
            idxs = [x+1 for x in idxs]
    else:
        idxs = []
        for tmp in string.split(','):
            if( '-' in tmp ):
                tmp2 = tmp.split('-')
                if( tmp2[0] == '' ):
                    tmp2[0] = 0
                    if( ones_based ):
                        tmp2[0] += 1
                if( tmp2[1] == '' ):
                    if( inclusive ):
                        tmp2[1] = list_len-1
                    else:
                        tmp2[1] = list_len
                    if( ones_based ):
                        tmp2[1] += 1
                if( inclusive ):
                    x = range(int(tmp2[0]),int(tmp2[1])+1)
                else:
                    x = range(int(tmp2[0]),int(tmp2[1]))
            else:
                x = [int(tmp)]
            idxs.extend(x)
        if( sort_unique ):
            idxs = sorted(set(idxs)) # remove duplicates and sort
        if( zero_pad > 0 ):
            for i in range(len(idxs)):
                idxs[i] = str(idxs[i]).zfill(zero_pad)
    return idxs


### Unused by maybe useful in future ####

def config_get_type(config, section, option):
    """parse config options including type like 'foo: :int: 3' """
    full_string = config.get(section, option)
    tmp = full_string.split(':',2)
    tmp = [x.strip() for x in tmp]
    if( len(tmp) == 3 and tmp[0] == '' ):
        if( tmp[1].lower() == 'bool' ):
            rv = tmp[2].lower() in ['true', '1', 't', 'y', 'yes']
        elif( tmp[1] in ['str', 'int', 'float', 'bool'] ):
            rv = eval(tmp[1]+'("'+tmp[2]+'")')
        else:
            err = "Unknown type '{}' for config in '[{}] {}: {}'".format(tmp[1], section, option, full_string)
            LOG.error(err)
            raise RuntimeError(err)
    else:
        # if no :type:, just return a string
        rv = full_string
    return rv

