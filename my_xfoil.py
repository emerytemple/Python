import matplotlib.pyplot as plt
import subprocess as sp
import numpy as np
import os
from math import *
from scipy import interpolate

import atm
from Convert import *

def airfoil(naca, Re, M):
    # delete previous runs if they exist
    try:
        os.remove('pya')
    except OSError:
        pass

    try:
        os.remove('pyb')
    except OSError:
        pass

    # prepare input commands

    ifh = ""  
    ifh += "naca " + naca + "\n"
    ifh += "plop\n"
    ifh += "g\n"
    ifh += "\n"
    ifh += "oper\n"
    ifh += "re   " + str(Re) + "\n"
    ifh += "mach " + str(M) + "\n"
    ifh += "visc\n"
    ifh += "iter 500\n"
    ifh += "pacc\n"
    ifh += "pya\n"
    ifh += "\n"
    ifh += "aseq -4.0 10.5 0.1\n"
    ifh += "\n"
    ifh += "quit\n"

    # run xfoil
    proc = sp.Popen(['xfoil.exe'], stdin=sp.PIPE, stdout=sp.PIPE)
    out = proc.communicate(input=ifh.encode('ascii'))[0]

    # read from file
    lines = open('pya').readlines()
    open('pyb', 'w').writelines(lines[12:])
    table = np.loadtxt('pyb')

    alphat = table[:,0]
    clt = table[:,1]
    cdt = table[:,2]
    cmt = table[:,3]

    # print(alphat, len(alphat))

    # for i in range(len(alphat)):
    #     print(alphat[i], clt[i], cdt[i])

    # chop arrays
    a = []
    cl = []
    cd = []
    zero_lift = float('NaN')
    for i in range(len(alphat)-1):
        if clt[i] < 0 and clt[i+1] > 0: # first point starts at zero lift
            zero_lift = interpolate('y', 0.0, alphat[i], clt[i], alphat[i+1], clt[i+1])
            # print(zero_lift)
            a.append(zero_lift)
            cl.append(0.0)
            cd.append(interpolate('x', zero_lift, alphat[i], cdt[i], alphat[i+1], cdt[i+1]))            
        if isfinite(zero_lift):
            if alphat[i] > zero_lift: # continue array
                a.append(alphat[i])
                # print(a[i])
                cl.append(clt[i])
                cd.append(cdt[i])
        if clt[i+1] < clt[i]: # end at stall
            stall = alphat[i]
            break
        if -0.25 < alphat[i] and alphat[i] < 0.25:
            clo = clt[i]

    anp = np.array(a)
    clnp = np.array(cl)
    cdnp = np.array(cd)

    for i in range(len(anp)):
        print(anp[i], clnp[i], cdnp[i])

    # calculate polynomials
    cl = np.polyfit(anp, clnp, 6)
    cd = np.polyfit(anp, cdnp, 6)

    cl2 = np.polyfit(alphat, clt, 6)
    cd2 = np.polyfit(alphat, cdt, 6)

    ao = np.polyder(cl)

    '''
    a = cl[0]
    b = cl[1]
    c = cl[2]

    zl = (-b+sqrt((b*b)-(4.0*a*c)))/(2*a)
    sss = c

    err1 = abs((zl-zero_lift)/zero_lift)*100
    err2 = abs((sss-clo)/clo)*100
    '''
    # print(zero_lift, zl, clo, sss, err1, err2)

    # plot values
    pcl = np.poly1d(cl)
    pcd = np.poly1d(cd)
    pcl2 = np.poly1d(cl2)
    pcd2 = np.poly1d(cd2)

    clp2 = np.zeros(len(alphat))
    cdp2 = np.zeros(len(alphat))
    for i in range(len(alphat)):
        clp2[i] = pcl2(alphat[i])
        cdp2[i] = pcd2(alphat[i])

    n = 100
    da = (stall-zero_lift)/(n-1)
    a = np.zeros(n)
    clp = np.zeros(n)
    cdp = np.zeros(n)
    for i in range(n):
        a[i] = zero_lift + (i*da)
        # print(i, a[i])
        clp[i] = pcl(a[i])
        cdp[i] = pcd(a[i])

    plt.subplot(2,2,1)
    plt.plot(alphat, clt, label='clt')
    plt.plot(alphat, clp2, label='clp2')
    # plt.plot(anp, clnp, label='clnp')
    # plt.plot(a, clp, label='clp')
    plt.legend()

    plt.subplot(2,2,2)
    plt.plot(alphat, cdt, label='clt')
    plt.plot(alphat, cdp2, label='clp2')    
    # plt.plot(anp, cdnp, label='clnp')
    # plt.plot(a, cdp, label='clp')
    plt.legend()

    plt.subplot(2,2,3)
    plt.plot(cdt, clt, label='cldt')
    plt.plot(cdp2, clp2, label='cldp2')
    # plt.plot(cdnp, clnp, label='cldnp')
    # plt.plot(cdp, clp, label='cldp')
    plt.legend()

    plt.show()

    return zl, stall, cl, cd

def interpolate(var, n, x1, y1, x2, y2):
    # still need to check x2 > x1 and x in bounds

    m = (y2-y1)/(x2-x1)
    b = y1 - (m*x1)

    if var is 'x':
        retval = (m*n) + b
    if var is 'y':
        retval = (n-b)/m

    return retval

def wing(cl, cd, AR, e, M, lam, zero_lift, stall):
    pcl = np.poly1d(cl)
    pcd = np.poly1d(cd)

    ao = np.polyder(pcl)

    n = 100
    da = (stall-zero_lift)/(n-1)
    alpha = np.zeros(n)
    CL = np.zeros(n)
    CD = np.zeros(n)
    ncl = np.zeros(n)
    for i in range(n):
        alpha[i] = zero_lift + (i*da)
        acl = ao(alpha[i])*cos(lam)
        ncl[i] = pcl(alpha[i])
        f = acl/(pi*AR)
        a = acl/(sqrt(1.0-((M*cos(lam))**2.0)+(f*f))+f)
        CL[i] = a*(alpha[i]-zero_lift)
        CD[i] = pcd(alpha[i]) + ((CL[i]*CL[i])/(pi*e*AR))  

    PCL = np.polyfit(alpha, CL, 2)
    PCD = np.polyfit(alpha, CD, 2)

    # plt.plot(alpha, ncl, alpha, CL)
    # plt.show()

    return CL, CD

def plot_airfoil_full(zero_lift, stall, cl, cd):
    # plot values
    pcl = np.poly1d(cl)
    pcd = np.poly1d(cd)

    n = 100
    da = (stall-zero_lift)/(n-1)
    a = np.zeros(n)
    clp = np.zeros(n)
    cdp = np.zeros(n)
    for i in range(n):
        a[i] = zero_lift + (i*da)
        print(i, a[i])
        clp[i] = pcl(a[i])
        cdp[i] = pcd(a[i])

    plt.subplot(2,2,1)
    plt.plot(a, clp)
    plt.title('cl vs AoA')
    plt.xlabel('AoA (deg)')
    plt.ylabel('cl')
    plt.grid(True)

    plt.subplot(2,2,2)
    plt.plot(a, cdp)
    plt.title('cd vs AoA')
    plt.xlabel('AoA (deg)')
    plt.ylabel('cd')
    plt.grid(True)

    plt.subplot(2,2,4)
    plt.plot(cdp, clp)
    plt.title('cl vs cd (Drag Polar)')
    plt.xlabel('cd')
    plt.ylabel('cl')
    plt.grid(True)

    plt.show()

class Weight():
    empty = Convert(256.0, 'lb', 'kg')
    fuel = Convert(5.0*6.82, 'lb', 'kg')
    crew = Convert(200.0, 'lb', 'kg')
    payload = 0
    total = empty + fuel + crew + payload

def main():
    # input
    h = 0.0 # km
    v = 12.5 # m/s
    c = Convert(1.0, 'ft', 'm')
    naca = '2408'
    b = Convert(12.0, 'ft', 'm')
    lam = Convert(0.0, 'deg', 'rad')

    # constants
    gamma = 1.4
    Rgas = 287.15

    # atmosphere
    P, T = atm.atmosphere(h)
    rho = atm.density(P, T)
    mu = atm.viscosity(T)
    a = sqrt(gamma*Rgas*T)
    M = v/a

    # calculations
    Re = (rho*v*c)/mu

    S = b*c # fix to include lam later
    AR = (b*b)/S
    e = 0.95

    zero_lift, astall, cl, cd = airfoil(naca, Re, M)

    # CL, CD = wing(cl, cd, AR, e, M, lam, zero_lift, astall)
    # plot_airfoil_full(ao, astall, cl, cd)

    # print(h, P, T, rho, mu, Re, a, M)

    # w = Weight()
    # print(w.empty, w.crew, w.fuel, w.total)

if __name__ == "__main__":
    main()
