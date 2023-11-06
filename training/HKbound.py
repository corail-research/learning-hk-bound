#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: augustinparjadis
"""

import sys
import ctypes
from numpy.ctypeslib import ndpointer

def thetas_C(file,N,UB):
    LP_c_char = ctypes.POINTER(ctypes.c_char)
    LP_LP_c_char = ctypes.POINTER(LP_c_char)
    
    N_NODES = N
    
    hkl = ctypes.CDLL("./hkfactors.so", mode = ctypes.RTLD_GLOBAL)
    hk = hkl.hkfactors
    hk.argtypes = [ctypes.c_int,LP_LP_c_char]
    hk.restype = ndpointer(dtype=ctypes.c_float, shape=(N_NODES+1,))
    argv = ["hkfactors","-filein",file,"-UB",str(UB)]
    argc = len(argv)
    
    p = (LP_c_char*len(argv))()
    for i, arg in enumerate(argv):
        enc_arg = arg.encode('utf-8')
        p[i] = ctypes.create_string_buffer(enc_arg)
    
    na = ctypes.cast(p, LP_LP_c_char)
    
    hk_res = hk(argc, na)
    
    print(f'\nHKbound = {hk_res[0]}')
    print(hk_res)
    
    return hk_res[0],hk_res[1:]
