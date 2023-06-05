#!/usr/bin/env python
# _*_ coding=utf-8 _*_

# This script extracts depth information from the bam file name, and echos a min supp score for bash.
# Example file name: SKBR3_ONT_sort_25X 
#
# Example execution in bash
## for file in SKBR3_*.bam
## do
## name=${basename ${file} .bam}
## MinSup=$(python string_info.py ${name})
## bcftools.....${MinSup}...
## done
#
#
# Table of depth and min supp
# Depth		Min Supp
# 112X		15
#  61X		12
#  30X		7
#  25X		6
#  20X		5
#  15X		4
#  10X		3
#   5X		2

import glob
import subprocess
import os
import sys
import re
import math

# check that the correct number of arguments are provided
if len(sys.argv) != 2:
  print("Error: Incorrect number of arguments provided.")
  print("Usage: string_info.py <filename>")
  sys.exit()

filename = sys.argv[1]

Tech = f"{filename}".split('_')[1]	    
DepthX = f"{filename}".split('_')[3]
if DepthX != "MAX": 
  Depth = re.sub('[^0-9]', '', DepthX)
elif Tech == "ONT":
  Depth = 61
  MinSup = 12
elif Tech == "PBCLR":
  Depth = 112   
  MinSup = 15
else:
  Depth = "error"

Depth = int(Depth)
if Depth < 40:
  MinSup = 0.2*Depth + 1

MinSup = int(MinSup)
print(MinSup)
