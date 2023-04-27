#!/usr/bin/python

__doc__ = """Wrapper around mutiprocessing pool to simplify some usage
@TCC--Some more usage documentation goes here
"""

import sys
from sys import stderr,stdout

import multiprocessing as MP
from multiprocessing.pool import RUN as MP_RUN

import misc_utils
from misc_utils import verbose

import inspect

mppool = None

def mp_pool():
    global mppool
    #print mppool
    #print inspect.stack()[1]
    return mppool


def CreateProcessPool(num_processes=None, worker_globals_dict=None, quiet_flag=False):
    """Create a slightly nicer multiprocessing pool"""
    if( num_processes == None or ( type(num_processes) is str and num_processes.lower() == 'none' ) ): 
        return None
    num_processes = int(num_processes)
    if( num_processes <= 0 ): 
        num_processes = MP.cpu_count() 
    if not quiet_flag:
        print("Setting up multiprocessing pool with",num_processes,"worker processes", file=sys.stderr)
    # create the Pool
    mp_pool = MP.Pool(processes=num_processes, initializer=_WorkerProcessInit, initargs=[worker_globals_dict])
    # save the num_processes in mp_pool (should be built into class, but isn't)
    mp_pool.num_processes = num_processes # @TCC: len(mp_pool._pool) works too
    # and for convience, attach helper/utility methods
    misc_utils.TransplantFunc2Method(RespawnProcessPool, mp_pool, 'Respawn') 
    misc_utils.TransplantFunc2Method(GoodProcessPool, mp_pool, 'Good') 
    global mppool
    #globals().setattr('mppool', mp_pool)
    mppool = mp_pool


def RespawnProcessPool(mp_pool, num_processes=None, worker_globals_dict=None):
    """Close down all the pool worker processes and restart them.
             num_processes -- if None, use the same number of processes as initial pool
         Useful when the workers need to have new global variables... safer than the SetGlobals functions.
         Also good for freeing up resources / leaks.
         NOTE: use like: mp_pool = mp_pool.Respawn(...)"""
    global mppool
    print("Respawn multiprocessing pool", file=sys.stderr)
    if( num_processes == None ):
        num_processes = mppool.num_processes
    mppool.close()
    mppool.join()
    mppool.terminate()
    CreateProcessPool(num_processes, worker_globals_dict, True)
    

def GoodProcessPool(mp_pool):
    return mppool._state == MP_RUN

### Pool worker process init function ###
def _WorkerProcessInit(variables_dict=None):
    # the variables in variable_dict, (name, value) pairs, are assigned globals (well to this module)
    # they can be accessed in a function spawned to the workers via globals()['mppool'].name
    if( variables_dict ):
        for k,v in variables_dict.items():
            globals()[k] = v


