# ----------------------------------------------------------------------
#   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
#   http://www.cs.sandia.gov/~sjplimp/spparks.html
#   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
#   Copyright (2008) Sandia Corporation.  Under the terms of Contract
#   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
#   certain rights in this software.  This software is distributed under 
#   the GNU General Public License.
#
#   See the README file in the top-level SPPARKS directory.
# -------------------------------------------------------------------------

# Python wrapper on SPPARKS library via ctypes

import sys,traceback
from ctypes import *

def _to_bytes(value):
  if isinstance(value,bytes): return value
  return str(value).encode('utf-8')

class spparks:
  def __init__(self,name="",cmdargs=None,comm=None):

    # load libspparks.so by default
    # if name = "g++", load libspparks_g++.so

    try:
      if not name: self.lib = CDLL("libspparks.so",RTLD_GLOBAL)
      else: self.lib = CDLL("libspparks_%s.so" % name,RTLD_GLOBAL)
    except:
      type,value,tb = sys.exc_info()
      traceback.print_exception(type,value,tb)
      raise Exception("Could not load SPPARKS dynamic library")

    # create an instance of SPPARKS
    # optionally pass an MPI communicator from mpi4py
    # no_mpi call lets SPPARKS use MPI_COMM_WORLD
    # cargs = array of C strings from args

    self.lib.spparks_open.argtypes = [c_int,POINTER(c_char_p),c_void_p,POINTER(c_void_p)]
    self.lib.spparks_open_no_mpi.argtypes = [c_int,POINTER(c_char_p),POINTER(c_void_p)]
    self.lib.spparks_close.argtypes = [c_void_p]
    self.lib.spparks_file.argtypes = [c_void_p,c_char_p]
    self.lib.spparks_command.argtypes = [c_void_p,c_char_p]
    self.lib.spparks_command.restype = c_char_p
    self.lib.spparks_extract.argtypes = [c_void_p,c_char_p]
    self.lib.spparks_energy.argtypes = [c_void_p]
    self.lib.spparks_energy.restype = c_double
    
    if cmdargs:
      cmdargs = list(cmdargs)
      cmdargs.insert(0,"spparks.py")
      cmdargs = [_to_bytes(arg) for arg in cmdargs]
      narg = len(cmdargs)
      cargs = (c_char_p*narg)(*cmdargs)
    else:
      narg = 0
      cargs = None

    self.spk = c_void_p()

    if comm is None:
      self.lib.spparks_open_no_mpi(narg,cargs,byref(self.spk))
    else:
      try:
        from mpi4py import MPI
      except:
        raise Exception("mpi4py is required to pass an MPI communicator")
      if hasattr(MPI,"Is_initialized") and not MPI.Is_initialized():
        MPI.Init()
      comm_handle = MPI._handleof(comm)
      try:
        comm_size = MPI._sizeof(comm)
      except:
        try:
          comm_size = MPI._sizeof(MPI.Comm)
        except:
          comm_size = sizeof(c_void_p)
      if comm_size == sizeof(c_int):
        self.lib.spparks_open.argtypes = [c_int,POINTER(c_char_p),c_int,POINTER(c_void_p)]
        comm_arg = c_int(comm_handle)
      else:
        self.lib.spparks_open.argtypes = [c_int,POINTER(c_char_p),c_void_p,POINTER(c_void_p)]
        comm_arg = c_void_p(comm_handle)
      self.lib.spparks_open(narg,cargs,comm_arg,byref(self.spk))
      self.comm = comm

  def __del__(self):
    if self.spk: self.lib.spparks_close(self.spk)

  def close(self):
    self.lib.spparks_close(self.spk)
    self.spk = None

  def file(self,file):
    file = _to_bytes(file)
    self.lib.spparks_file(self.spk,file)

  def command(self,cmd): 
    cmd = _to_bytes(cmd)
    self.lib.spparks_command(self.spk,cmd)

  def extract(self,name,type):
    name = _to_bytes(name)
    if type == 0:
      self.lib.spparks_extract.restype = POINTER(c_int)
      ptr = self.lib.spparks_extract(self.spk,name)
      return ptr[0]
    if type == 1:
      self.lib.spparks_extract.restype = POINTER(c_int)
      ptr = self.lib.spparks_extract(self.spk,name)
      return ptr
    if type == 2:
      self.lib.spparks_extract.restype = POINTER(POINTER(c_int))
      ptr = self.lib.spparks_extract(self.spk,name)
      return ptr
    if type == 3:
      self.lib.spparks_extract.restype = POINTER(c_double)
      ptr = self.lib.spparks_extract(self.spk,name)
      return ptr[0]
    if type == 4:
      self.lib.spparks_extract.restype = POINTER(c_double)
      ptr = self.lib.spparks_extract(self.spk,name)
      return ptr
    if type == 5:
      self.lib.spparks_extract.restype = POINTER(POINTER(c_double))
      ptr = self.lib.spparks_extract(self.spk,name)
      return ptr
    return None

  def extract_numpy(self,name,type,shape=None):
    try:
      import numpy as np
    except:
      raise Exception("NumPy is required for extract_numpy()")

    if isinstance(name,bytes): name_str = name.decode('utf-8')
    else: name_str = name

    if type == 0 or type == 3:
      return self.extract(name,type)

    ptr = self.extract(name,type)
    if ptr is None: return None

    if type == 1 or type == 4:
      if shape is None:
        shape = (self.extract("nlocal",0),)
      return np.ctypeslib.as_array(ptr,shape=shape)

    if type == 2 or type == 5:
      if shape is None:
        if name_str == "xyz":
          shape = (self.extract("nlocal",0),3)
        else:
          raise Exception("shape required for array extraction")
      base = ptr[0]
      flat = np.ctypeslib.as_array(base,shape=(shape[0]*shape[1],))
      return flat.reshape(shape)

    return None

  def energy(self):
    self.lib.spparks_energy.restype = c_double
    return self.lib.spparks_energy(self.spk)
