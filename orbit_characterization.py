
import sys
from numpy import *
import spiceypy as sp

UMG = 1.989e+43           # g
ULG = 3.085678e+21      # cm
UVG = 1.0e+05               # cm/s
UTG = ULG/UVG               # s 

print UTG

GMKS = 6.672e-11             # m3 kg-1 s-2

GCGS = GMKS*(100**3)/1000.0 

print GCGS
print UTG/60.0/60.0/24.0/365.0

print GCGS*UMG*(UTG**2)/(ULG**3)

ULKM = ULG/100.0/1000.0

G = 43007.1
mHost = 7.83318368e+01*UMG
mSat = 3.18298641e+00*UMG

mu = GCGS*(mHost+mSat) # cgs , cm3 s-2

mu = mu/(100**3)/(1000**3) # km3 s-2

#print mu

indexIni = int(sys.argv[2])
indexFin = int(sys.argv[3])

print "\nSearching orbital elements from %lf to %lf\n" % (indexIni,indexFin)

time, x, y, z, vx, vy, vz = loadtxt(sys.argv[1], unpack=True)

x = x*ULG/100/1000
y = y*ULG/100/1000
z = z*ULG/100/1000

peri = zeros(len(time))
ecc =  zeros(len(time))
inclination = zeros(len(time))

for i in range(len(time)):
    #print time[i], x[i], y[i], z[i], vx[i], vy[i], vz[i]
    elements = sp.oscelt([x[i], y[i], z[i], vx[i], vy[i], vz[i]], time[i],mu)
    #print elements[0],elements[1],elements[2],elements[3],elements[4],elements[5],elements[6],elements[7]
    peri[i] = elements[0]
    ecc[i] = elements[1]
    inclination[i] = elements[2]

filter = logical_and( time>=indexIni , time<=indexFin )

print filter, peri[filter]*1000*100/ULG

print peri*1000*100/ULG
print mean(peri[filter])*1000*100/ULG, mean(ecc), mean(inclination)*180.0/pi
