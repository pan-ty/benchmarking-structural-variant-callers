#!/home/tpg8911/.conda/envs/SVcalling/bin/python
# _*_ coding=utf-8 _*_

# This script performs 
import glob
import subprocess
import os
import sys
import re
import pandas 

# check that the correct number of arguments are provided
if len(sys.argv) != 2:
  print("Error: Incorrect number of arguments provided.")
  print("Usage: string_info.py </path/to/benchmarking/folder/>")
  print("Example: /projects/b1171/tpg8911/SKBR3/benchmarking/")
  sys.exit()

# SV tool must be follow their path
pathname = sys.argv[1]
os.chdir(pathname)

# Empty list for pandas df
Caller_list = []
Tech_list = []
DepthX_list = []
Precision_list = []
Recall_list = []
F1_list = []

# Loop through benchmarking results folders
for path in glob.glob(f"{pathname}/*"):
  ID = os.path.basename(path)
  Caller = f"{ID}".split('_')[0]
  Caller_list.append(Caller)
  Tech = f"{ID}".split('_')[2]
  Tech_list.append(Tech)
  DepthX = f"{ID}".split('_')[4]
  DepthX_list.append(DepthX)
  os.chdir(pathname+ID)
  with open('summary.json', 'r') as file:
    for line in file:
      line.strip()
      if line.find('precision') != -1:
        precision = re.findall(r"[-+]?\d*\.\d+|\d+", line)
        Precision_list.append(precision[0])
      if line.find('recall') != -1:
        recall = re.findall(r"[-+]?\d*\.\d+|\d+", line)
        Recall_list.append(recall[0])
      if line.find('f1') != -1:
        f1 = re.findall(r"[-+]?\d*\.\d+|\d+", line)
        F1_list.append(f1[1])
	
dict = {'Caller': Caller_list, 'Tech': Tech_list, 'Depth': DepthX_list, 'Precision': Precision_list, 'Recall': Recall_list, 'F1': F1_list}
df = pandas.DataFrame(dict)
df = df.sort_values(by=['Caller','Tech','Depth'])
df.to_csv('/home/tpg8911/yanglab/ANALYSIS/benchmarking_results.csv', sep="\t")
