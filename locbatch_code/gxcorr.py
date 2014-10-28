#!/usr/bin/python

__doc__ = """Compute a generalized cross-correlation using FFTs
@TCC--Some more usage documentation goes here
"""

import sys
import math 
import numpy as NP

from numpy.fft import rfft, irfft, fft, ifft
#from scipy.fftpack import rfft, irfft, fft, ifft # a bit faster but doesn't seem to work quite right

## utility ##
def nextpow2(x):
    """Return 2**N for the first integer N such that 2**N >= abs(x)"""
    return 2**NP.ceil(NP.log2(NP.abs(x)))


## Main buisness #####

def GXCorr(x, y, maxlag, mode='coeff-real', comp_cenv=True):
    # assumes x is pattern and y is target (y should be >= x in length)
    NEXTPOW2_PADDING = True # both numpy and scipy.fftpack are much much faster with this set to True as of Feb2011
    SHIFT_OFFSET = -1
    CENTER_FLAG = True # if not True, have adjust lags at higher level
    maxlag = math.ceil(maxlag)

    # check inputs for the correct type
    assert( x.ndim == 1 )
    assert( y.ndim == 1 )

    ## prepad x if CENTER_FLAG and it is shorter than y
    if( CENTER_FLAG and len(x) < len(y) ):
        # prepad x & postpad x
        tmp = len(y)-len(x)
        if( tmp%2 == 1 ):
            print >>sys.stderr, 'WARNING: gxcorr trying to do symmetric padding with odd length difference'
        tmp = int(tmp/2)
        x = NP.concatenate( ( NP.zeros(tmp), x, NP.zeros(tmp) ), 1)
    else:
        x = NP.concatenate( ( x, NP.zeros(len(y)-len(x)) ), 1)

    # pad out signals to make them same length
    if( NEXTPOW2_PADDING ):
    # final size will be the next power of 2 of y (the longer signal)
        padded_size = nextpow2(len(y))
        y = NP.concatenate((y, NP.zeros(padded_size-len(y))), 1)
    # pad x to same length as y
    if( len(x) < len(y) ):
        x = NP.concatenate((x, NP.zeros(len(y)-len(x))), 1)
    if( len(y) < len(x) ):
        y = NP.concatenate((y, NP.zeros(len(x)-len(y))), 1)

    # convience
    N = len(x)

    ## blocks dependent on mode 
    # @TCC -- Just 'coeff' mode for now
    if( mode == 'coeff' ): # coeff normalization using FFT
        SHIFT_OFFSET = -1
        Fx = fft(x)
        Fy = fft(y)
        Pxy = NP.conj(Fx)*Fy 
        c = NP.real(ifft(Pxy))
        lag_idxs = NP.concatenate( (NP.arange(N-1-SHIFT_OFFSET-maxlag, N, dtype=NP.int), \
                                                            NP.arange(0, maxlag-SHIFT_OFFSET, dtype=NP.int)) )
        full_c = c # if needed below
        c = c[lag_idxs]
        tmp = NP.sqrt( NP.dot(x,x) * NP.dot(y,y) )
        if( tmp == 0 ):
            print >>sys.stderr, "Warning: GXCorr coeff normalization power == 0 ?!"
            c[:] = 0
        else:
            norm_factor = tmp
            c /= norm_factor
        maxlag_out = int(len(c)/2)

    elif( mode == 'coeff-real' ): # coeff normalization assuming real input
        SHIFT_OFFSET = 0
        Fx = rfft(x[::-1]) # flip x around
        Fy = rfft(y)
        Pxy = Fx * Fy # no need for conj
        c = irfft(Pxy)
        lag_idxs = NP.concatenate( (NP.arange(N-1-SHIFT_OFFSET-maxlag, N, dtype=NP.int), \
                                                            NP.arange(0, maxlag-SHIFT_OFFSET, dtype=NP.int)) )
        full_c = c 
        c = c[lag_idxs]
        tmp = NP.sqrt( NP.dot(x,x) * NP.dot(y,y) )
        if( tmp == 0 ):
            print >>sys.stderr, "Warning: GXCorr coeff-real normalization power == 0 ?!"
            c[:] = 0
        else:
            norm_factor = tmp
            c /= norm_factor
        maxlag_out = int(len(c)/2)

    # compute the hilbert amplitude envelope (optional but default)
    cenv = None
    if( comp_cenv ):

        if( True ): # use own hilbert and do clever padding (up to nextpow2 and from full_c), fastish [FAVORED METHOD]
            n = nextpow2(len(c))
            pre_pad = int((n-len(c))/2)
            pre_pad = min(pre_pad, lag_idxs[0]) 
            post_pad = int(n-len(c)-pre_pad) 
            tmp = full_c[NP.concatenate( (NP.arange(N-1-SHIFT_OFFSET-maxlag-pre_pad, N, dtype=NP.int), \
                                                        NP.arange(0, maxlag-SHIFT_OFFSET+post_pad, dtype=NP.int)) )] / norm_factor
            assert( len(tmp) == n )
            if( True ): # rfft (very tiny bit faster)
                h = rfft(tmp)
                h = NP.concatenate( ( h, NP.zeros(len(tmp)-len(h)) ), 1) 
                h[1:n/2] *= 2
            else: # fft (probably safer if using different fft libaray)  
                h = fft(tmp)
                h[1:n/2] *= 2
                h[n/2+1:] = 0.
            h = ifft(h)
            cenv = NP.sqrt(h*NP.conj(h))
            cenv = cenv[pre_pad:len(c)+pre_pad]

        if( False ): # use own hilbert method on existing Pxy, medium speed
            if( len(Pxy) < len(full_c) ): # rfft, so pad out
                tmp = NP.concatenate( ( Pxy, NP.zeros(len(full_c)-len(Pxy)) ), 1) 
            n = len(tmp)
            tmp[1:n/2] *= 2.
            tmp = ifft(tmp) # only one extra ifft needed to compute hilbert transform
            tmp = NP.sqrt(tmp*NP.conj(tmp))
            cenv = tmp[lag_idxs] / norm_factor

        if( False ): # use our own hilbert method, dumb padding, fastish 
            padding = 10  # extra padding to clean up the ends a bit 
            if( padding > 0 ):
                tmp = NP.concatenate( ( NP.zeros(padding), c, NP.zeros(padding) ), 1) 
            else:
                tmp = c
            n = nextpow2(len(tmp))
            if( len(tmp) < n ): 
                tmp = NP.concatenate( ( tmp, NP.zeros(n-len(tmp)) ), 1) 
            h = fft(tmp)
            h[1:n/2] *= 2
            h[n/2+1:] = 0.
            h = ifft(h)
            cenv = NP.sqrt(h*NP.conj(h))
            cenv = cenv[padding:len(c)+padding]

    return c, NP.real(cenv), maxlag_out


