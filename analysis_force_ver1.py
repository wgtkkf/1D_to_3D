# Force results files analysis code
# Last modified: March 08 2019
# https://qiita.com/segavvy/items/51a515c19bcd29b13b7f
# Coded by Takuro TOKUNAGA

import math
import numpy as np
import cmath
import time
import sys
import pandas as pd
import glob
start = time.time()

# functions for begin & finish
def begin():
    print ("begin")

def finish():
    print ("finish")

# main
begin()

# total force
txt_files = glob.glob("../surface/force_total/*.txt")
list = []
for f in txt_files:
    list.append(pd.read_csv(f))
df = pd.concat(list, axis=1)
df.to_csv("../surface/force_total.txt")

# Lennard-Jones force
txt_files = glob.glob("../surface/force_lj/*.txt")
list = []
for f in txt_files:
    list.append(pd.read_csv(f))
df = pd.concat(list, axis=1)
df.to_csv("../surface/force_lj.txt")

# Casimir force
txt_files = glob.glob("../surface/force_casimir/*.txt")
list = []
for f in txt_files:
    list.append(pd.read_csv(f))
df = pd.concat(list, axis=1)
df.to_csv("../surface/force_casimir.txt")

# Capacitor force
txt_files = glob.glob("../surface/force_capacitor/*.txt")
list = []
for f in txt_files:
    list.append(pd.read_csv(f))
df = pd.concat(list, axis=1)
df.to_csv("../surface/force_capacitor.txt")

# Coulomb force
txt_files = glob.glob("../surface/force_coulomb/*.txt")
list = []
for f in txt_files:
    list.append(pd.read_csv(f))
df = pd.concat(list, axis=1)
df.to_csv("../surface/force_coulomb.txt")

# file close
finish()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
