import numpy as np
import h5py
import argparse, sys

parser = argparse.ArgumentParser(prog="make_hdf5_ref",
                                 description="It converts a TChem-atm output to HDF5 format.")

parser.add_argument('-r_file','--ref_file',type=str,required=True, default="",
 help="Name of the TChem-atm output to be converted to HDF5 format.")


args = parser.parse_known_args(args=sys.argv)
ref_file = args[0].ref_file

data = np.genfromtxt(ref_file, dtype=str)
ref = (data[1:,:]).astype(float)

new_filename= ref_file.replace('.dat', '')
file_name_hdf5=new_filename+".hdf5"
with h5py.File(file_name_hdf5, 'w') as hdf:
    hdf.create_dataset("ref", data=ref)
