import numpy as np
import math
from bhmie import bhmie
import matplotlib.pyplot as plt

a = 10
wavelength = 0.532


ref_med = 1.00
ref_part = np.complex(1.4,0.02)
refrel = ref_part/ref_med

nang = 180

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

dang = 0.5*np.pi/(nang-1)

# plot results
plt.polar(theta*np.pi, np.log10(s11), color='b')
plt.polar(-theta*np.pi, np.log10(s11), color='b')
#plt.title('mie scattering\n')

'''
# save plot to file
filename = 'mie_polar_' + str(a) + 'micron_ref_med_' + str(ref_med) + '.png'
plt.imsave(filename)
'''
plt.show()


