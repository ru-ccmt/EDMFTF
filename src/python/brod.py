#!/usr/bin/env python
# @Copyright 2007 Kristjan Haule
from scipy import *
from scipy import integrate, interpolate
# from pylab import *
import brd
import time

def Broad(width, kwidth, om, fw):
    " Broadens the data with gaussian of width=width"
    def MakeTanMesh(N, tanc, tanw, b0, b1):
        if not(b0<b1): print "Relation must hold: b0<b1!"
        if not(b0<tanw and tanw<b1): print "Relation mesu hold: b0<tanw<b1!"
        if not(b0>0): print "b0 must be positive!"
        du = arctan(((tanc-b0)/tanw))
        b1n = arctan((b1-tanc)/tanw)+du
        m0 = [tanc + tanw * tan(b1n*(i-1)/(N-2)-du) for i in range(1,N)]
        return hstack( (-array(m0[::-1]), array([0]+m0) ) )
    
    fwi = interpolate.interp1d(om, fw)
    fwn=[]
    for im in range(len(om)):
        w=width + kwidth*abs(om[im])
        if (om[im]-om[0]>w*4 and om[-1]-om[im]>w*4):  # Gaussian is fully within existing mesh
            x = brd.maketanmesh(200,0.0,w,w/50,w*20)
            x2,ni = brd.combinemesh(om[im],om,x)
            eps = x2[:ni]
            x3 = om[im]-eps
            tw = (2*w**2)
            gs = exp(-x3**2/tw)/(sqrt(2*pi)*w)
            norm = integrate.trapz(gs,x=x3)
            yn = integrate.trapz(fwi(eps) * gs, x=eps)/abs(norm)
        else:
            yn = fw[im]
        fwn.append(yn)
    return array(fwn)

if __name__ == '__main__':

    data = loadtxt('S2.dat').transpose()
    om = data[0]

    y=[]
    for i in range(1,(len(data)+1)/2):
        y.append( Broad(0.05, 0.1, om, data[2*i-1]+data[2*i]*1j) )
        
    plot(om, data[1])
    plot(om, real(y[0]))
    plot(om, data[2])
    plot(om, imag(y[0]))

    plot(om, data[15])
    plot(om, real(y[7]))
    plot(om, data[16])
    plot(om, imag(y[7]))

    show()
    
    
