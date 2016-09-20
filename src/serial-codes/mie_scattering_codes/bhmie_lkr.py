import numpy as np

def bhmie(x,refrel,nang):
# This file is converted from mie.m, see http://atol.ucsd.edu/scatlib/index.htm
# Bohren and Huffman originally published the code in their book on light scattering

# Calculation based on Mie scattering theory
# input:
#      x      - size parameter = k*radius = 2pi/lambda * radius
#                   (lambda is the wavelength in the medium around the scatterers)
#      refrel - refraction index (n in complex form for example:  1.5+0.02*i;
#      nang   - number of angles for S1 and S2 function in range from 0 to pi/2
# output:
#        S1, S2 - funtion which correspond to the (complex) phase functions
#        Qext   - extinction efficiency
#        Qsca   - scattering efficiency
#        Qback  - backscatter efficiency
#        gsca   - asymmetry parameter


    max_nang = 100000
    nmxx = 150000

    if(nang>max_nang):
        print " error : nang > max_nang"

    if(nang<2):
        nang = 2

    dx = x
    drefrl = refrel
    y = x*drefrl
    ymod = abs(y)

    #    Series expansion terminated after NSTOP terms
    #    Logarithmic derivatives calculated from NMX on down
    xstop = x + 4*x**0.3333 + 2
    nmx = int(np.max([xstop,ymod])) + 15

    # BTD experiment 91/1/15: add one more term to series and compare resu<s
    #      NMX=AMAX1(XSTOP,YMOD)+16
    # test: compute 7001 wavelen>hs between .0001 and 1000 micron
    # for a=1.0micron SiC grain.  When NMX increased by 1, only a single
    # computed number changed (out of 4*7001) and it only changed by 1/8387
    # conclusion: we are indeed retaining enough terms in series!

    nstop = int(xstop)

    nmx = 300
    nstop = 300

    if(nmx > nmxx):
        print "Error : nmx>NMXX"

    dang = 0
    if(nang > 1):
        dang = 0.5*np.pi/(nang-1)

    theta = 0.5*np.pi*np.linspace(0.0,1.0,nang,endpoint=True)
    amu = np.cos(theta)

    pi0 = np.zeros(nang)
    pi1 = np.ones(nang)
    pi = np.zeros(nang)
    tau = np.zeros(nang)

    nn = 2*nang -1
    s1 = np.zeros(nn,dtype=complex)
    s2 = np.zeros(nn,dtype=complex)

    # Logarithmic derivative D(J) calculated by downward recurrence
    # beginning with initial value (0.,0.) at J=NMX

    d = np.zeros(nmx+1,dtype=complex)
    d[nmx] = 0.0
    nn = nmx-1

    for n in range(1,int(nmx+1)):
        en = nmx-n+1.0
        d[nmx-n] = en/y - 1.0/(d[nmx-n+1]+en/y)

    #*** Riccati-Bessel functions with real argument X
    #    calculated by upward recurrence
    psi0 = np.cos(dx)
    psi1 = np.sin(dx)
    chi0 = -np.sin(dx)
    chi1 = np.cos(dx)
    xi1 = psi1 -1j*chi1
    qsca = 0.0
    gsca = 0.0
    p = -1.0

    for n in range(1,int(nstop+1)):
        en = n
        fn = (2.0*en + 1)/(en*(en+1.0))

        psi = (2.0*en - 1.0)*psi1/dx - psi0
        chi = (2.0*en - 1.0)*chi1/dx - chi0
        xi = complex(psi,-chi)


        # Compute AN and BN
        an = (d[n]/drefrl + en/dx)*psi - psi1
        an /= (d[n]/drefrl + en/dx)*xi - xi1
        bn = (drefrl * d[n] + en/dx)*psi - psi1
        bn /= (drefrl*d[n] + en/dx)*xi - xi1

        if(n>1):
            an1 = an
            bn1 = bn

        # Augment sums for Qsca and g=<cos(theta)>
        qsca += (2.0*en + 1.0)*(abs(an)**2 + abs(en)**2)
        gsca += (2.0*en + 1.0)/(en*(en + 1.0)) * (an.real*bn.real + an.imag*bn.imag)

        if(n>1):
            gsca += ((en-1.0)*(en+1.0)/en) * (an1.real*an.real + an1.imag*an.imag + bn1.real*bn.real + bn1.imag*bn.imag)

        #*** Now calculate scattering intensity pattern
        #    First do angles from 0 to 90
        for j in range(0,nang):
            pi[j] = pi1[j]
            tau[j] = en*amu[j]*pi[j] - (en+1.0)*pi0[j]
            s1[j] += fn * (an*pi[j] + bn*tau[j])
            s2[j] += fn * (an*tau[j] + bn*pi[j])

        #*** Now do angles greater than 90 using PI and TAU from
        #    angles less than 90.
        #    P=1 for N=1,3,...; P=-1 for N=2,4,...
        p = -p
        for j in range(0,nang-1):
            jj = 2*nang - 1 - j - 1
            s1[jj] += fn*p*(an*pi[j] - bn*tau[j])
            s2[jj] += fn*p*(bn*pi[j] - an*tau[j])

        psi0 = psi1
        psi1 = psi
        chi0 = chi1
        chi1 = chi
        xi1 = complex(psi1,-chi1)

        # C*** Compute pi_n for next value of n
        # C    For each angle J, compute pi_n+1
        # C    from PI = pi_n , PI0 = pi_n-1
        for j in range(0,nang):
            pi1[j] = ((2*en + 1.0) * amu[j]*pi[j] - (en+1)*pi0[j])/en
            pi0[j] = pi[j]


    # C*** Have summed sufficient terms.
    # C    Now compute QSCA,QEXT,QBACK,and GSCA
    gsca = 2.*gsca/qsca
    qsca = 2.0/dx**2*qsca
    qext = 4.0/dx**2*s1[0].real
    qback = (abs(s1[2*nang-1-1])/dx)**2/np.pi

    return s1,s2,qext,qsca,qback,gsca,theta

