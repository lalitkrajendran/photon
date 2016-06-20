import numpy as np
import math
#import bhmie_2
from bhmie import bhmie

a = 10
wavelength = 0.532


ref_med = 1.33
ref_part = np.complex(1.4,0.02)
refrel = ref_part/ref_med

nang = 10

x = 2*np.pi*a*ref_med/wavelength

s1 = np.zeros(2*nang-1)
s2 = np.zeros(2*nang-1)
qext = 0
qsca = 0
qback = 0
gsca = 0

#[s1,s2,qext,qsca,qback,gsca] = bhmie_2.bhmie(x,refrel,5.0,s1,s2,qext,qsca,qback,gsca)
[s1,s2,qext,qsca,qback,gsca,theta] = bhmie(x,refrel,nang)


s11 = 0.5*(abs(s1)**2 + abs(s2)**2)
s12 = 0.5*(abs(s2)**2 - abs(s1)**2)

import os
print os.path.realpath(__file__)

b = 1