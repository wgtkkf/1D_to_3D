# file rename program
# Last modified: July 03 2019
# Coded by Takuro TOKUNAGA
# Correctly working

import numpy as np
import os
import glob
import time
start = time.time()

# functions for begin & finish
def begin():
    print ("begin")

def finish():
    print ("finish")

# main
begin()

flist = glob.glob('../surface/bumpy/*.txt')
#flist = glob.glob('../surface/bumpy_tip/*.txt')

# rename
i=0
for file in flist:
    os.rename(file, '../surface/bumpy/rename/data' + str(i) + '.txt')
    #os.rename(file, '../surface/bumpy_tip/rename/data' + str(i) + '.txt')
    i=i+1

finish()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
