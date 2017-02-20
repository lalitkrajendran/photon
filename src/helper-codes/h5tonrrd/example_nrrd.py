# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 08:42:37 2016

@author: lalit
"""

import numpy as np
import nrrd

# some sample numpy data
data = np.zeros((5,4,3))
filename = 'testdata.nrrd'

# write to a nrrd file
nrrd.write(filename, data)

# read the data back from file
