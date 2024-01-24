# Copyright (c) 2019-2023 Google LLC (USA).
# All rights reserved.
#
# This file is part of CGAL (www.cgal.org).
#
# $URL$
# $Id$
# SPDX-License-Identifier: GPL-3.0-or-later
#
#
# Author(s)     : Pierre Alliez
#                 Michael Hemmer
#                 Cedric Portaneri
#
#!/usr/bin/python

import os, sys, subprocess, datetime, time, getopt

def compute_quality_benchmark_data(execname, filename, alpha):

 output = ""
 cmd = (execname, "-i",
        filename, "-a", alpha)
 proc = subprocess.Popen(
     cmd,
     stdout=subprocess.PIPE,
     stderr=subprocess.PIPE,
     start_new_session=True)

 outs, errs = proc.communicate()
 output = outs.decode("utf-8") + errs.decode("utf-8")

 print(output)

def main(argv):
  execname=""
  filename=""
  alpha=""
  try:
    opts, args = getopt.getopt(sys.argv[1:], 'e:i:a:')
  except getopt.GetoptError:
    sys.exit(2)
  for opt, arg in opts:
    if opt == "-e":
      execname = arg
    elif opt == "-i":
      filename = arg
    elif opt == "-a":
      alpha = arg

  compute_quality_benchmark_data(execname, filename, alpha)

if __name__ == "__main__":
  main(sys.argv[1:])
