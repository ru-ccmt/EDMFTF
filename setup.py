#! /usr/bin/env python
# -*- coding: utf-8 -*-

###
#
# @file setup.py
#
#  DMFT is a software package provided by Rutgers Univiversity,
#  the State  University of New Jersey
#
# @version 1.0.0
# @author Kristjan Haule and Viktor Oudovenko
# @date 2016-02-15
#
###

__author__ = "Kristjan Haule and Viktor Oudovenko"
__version__ = "1.0.0"
__email__ = "haule@rutgers.edu"
__date__ = "Feb 28, 2016"

import sys
import os
import subprocess
import urllib

from script.dmft_install import Dmft_install
from script.blas         import Blas
from script.lapack       import Lapack
from script.fftw         import Fftw
from script.gsl          import Gsl

import configure

VERSION_MAJOR = 1
VERSION_MINOR = 1
VERSION_MICRO = 0


def main(argv):

  ### History of executed commands will be stored in  log.config
  logdir = 'log'
  if(not os.path.isdir(logdir)):
      print"Creating directory", logdir
      os.mkdir(logdir)
      # os.chdir(logdir)

  cmd = ""
  for arg in argv:
      cmd += arg+" "
  cmd += "\n"
  fp = open("log/log.config",'a')
  fp.write(cmd)
  fp.close()

  try:
    py_ver = sys.version_info
    print("\nDetected Python version %s" % ".".join(["%s" % i for i in py_ver]))
    if py_ver < (2, 7) or py_ver >= (2, 8):
        print("Python version 2.7+ required. Download and install the necessary "
              "python version from http://www.python.org/download/.")
        sys.exit(-1)
  except:
    print("\n Python version 2.7+ required. Download and install the necessary "
          "python version from http://www.python.org/download/.")
    sys.exit(-1)

  try:
    import numpy
    #from numpy.distutils.misc_util import get_numpy_include_dirs
    print("Detected numpy  version {}".format(numpy.__version__))
  except ImportError:
    print("numpy.distutils.misc_util cannot be imported. Please install ...")
    #subprocess.call(["pip", "install", "-q", "numpy>=1.8.0"])
    #from numpy.distutils.misc_util import get_numpy_include_dirs


  try:
    import scipy
    #from numpy.distutils.misc_util import get_numpy_include_dirs
    print("Detected scipy  version {}".format(scipy.__version__))
  except ImportError:
    print("scipy module cannot be imported. Please install ...")
    subprocess.call(["pip", "install", "-q", "scipy>=0.14.0"])
    #from numpy.distutils.misc_util import get_numpy_include_dirs

  print "\n"

  config = configure.Config((VERSION_MAJOR, VERSION_MINOR, VERSION_MICRO))
  dmft_install = Dmft_install(argv, config)

  #  if dmft_install.downblas :
  Blas(config, dmft_install)

  #  if dmft_install.downlapack :
  Lapack(config, dmft_install)

  #  if dmft_install.downfftw :
  Fftw(config, dmft_install)

  #  if dmft_install.downgsl :
  Gsl(config, dmft_install)

  dmft_install.resume()
  
  return 0


if "__main__" == __name__:
  sys.exit(main(sys.argv))  
