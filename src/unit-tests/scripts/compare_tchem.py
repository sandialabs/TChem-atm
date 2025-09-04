#!/usr/bin/env python
# coding: utf-8


import numpy as np
import numpy.linalg as lin
import argparse,sys,os
import h5py

def norms(x_comp, x_ref) :
    """norms(x_comp, x_ref) - Returns L1, L2, and Linf norms for x_comp - x_ref."""
    diff = np.subtract(x_comp, x_ref)
    if 1 == diff.ndim :
      axis=None
    else :
      axis=0
    L1 = lin.norm(diff,1,axis)
    L2 = lin.norm(diff,2,axis)
    Linf = lin.norm(diff,np.inf,axis)
    return (L1, L2, Linf)

parser = argparse.ArgumentParser(prog="compare_tchem",
description="It takes two files and compares them using the relative root mean square error. If the differences are larger than the error threshold, this program produces an error.")
parser.add_argument('-r_file','--ref_file',type=str,required=True, default="", help="File with reference data, including its path. This file must be saved in HDF5 format if use_hdf5 is set to True.")
parser.add_argument('-t_file','--test_file',type=str,required=True, default="", help="File with test data, including its path")
parser.add_argument('-error','--error_threshold',type=float,required=False, default=1e-6, help="Threshold to pass a test: relative error")
parser.add_argument('-use_hdf5','--use_hdf5',type=bool,required=False, default=True, help="To use HDF5 format for the reference file. ")
parser.add_argument('-check_norms','--check_norms',type=bool,required=False, default=False, help="Check norms")
args = parser.parse_known_args(args=sys.argv)


check_norms=args[0].check_norms
error_threshold=args[0].error_threshold
ref_file = args[0].ref_file
test_file= args[0].test_file
use_hdf5= args[0].use_hdf5
## We consider any number smaller than small_number to be extremely small.
small_number=1e-23


if use_hdf5 :
    with h5py.File(ref_file, 'r') as hdf:
        ref = hdf['ref'][:]
else :
    data = np.genfromtxt(ref_file, dtype=str)
    #do not use header
    ref = (data[1:,:]).astype(float)

data = np.genfromtxt(test_file, dtype=str)
#do not use header
test = (data[1:,:]).astype(float)

# only works for 2D arrays.
d1, d2 = np.shape(ref)

L1, L2, Linf = norms(ref, test)
# let's use L2 if the ref array to compute relative error.
ref_L2 = lin.norm(ref,2,0)

# We assume that all test are passing.
pass_all_tests= np.full(d2, True)


if check_norms:
    # check for the bit-for-bit case because the relative errors will be NaN
    for i_out in range(d2):
        if L1[[i_out]] == 0 and L2[[i_out]] == 0 and Linf[[i_out]] == 0:
            pass_all_tests[[i_out]] = True
            continue
        max4norm = 1.0
        ## Avoid division by zero or very-small numbers.
        if ref_L2[[i_out]] > small_number:
            max4norm = ref_L2[[i_out]]
        rel_error = L1[i_out]/ max4norm
        print("L1 rel_error",rel_error)
        if rel_error > error_threshold: pass_all_tests[i_out] = False
        rel_error = L2[i_out]/ max4norm
        print("L2 rel_error",rel_error)
        if rel_error > error_threshold: pass_all_tests[i_out] = False
        rel_error = Linf[i_out]/ max4norm
        print("Linf rel_error",rel_error)
        if rel_error > error_threshold: pass_all_tests[i_out] = False
print(f'final pass array = {pass_all_tests}')

#
pass_test = np.all(pass_all_tests)
if (pass_test== False):
  print('test', 'ref', 'diff', 'rel')
  for i in range(d1):
    for j in range(d2):
      diff =test[i,j]-ref[i,j]
      rel = diff/test[i,j]
      if abs(rel) > error_threshold:
        print(test[i,j], ref[i,j], diff, rel)
assert(pass_test)
