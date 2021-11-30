#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 11:49:25 2021

@author: rogerslc
"""
import pylab as plt
import numpy as np
import pandas as pd
import proposal as pp  #installed with pip
import scipy


#to read file
pdmtn=pd.read_hdf('/home/lrogers/Proposal/MountainProfile.h5')
StepSize=100

eps=0.1   # This is a trick to stop divide by zero errors

m_to_cm=100
GeV=1000


#putting coordinates in format for contour plot
Ymax=int(np.round(pdmtn['Y'].max(),0))
Xmax=int(np.round(pdmtn['X'].max(),0))
Ymin=int(np.round(pdmtn['Y'].min(),0))
Xmin=int(np.round(pdmtn['X'].min(),0))



vals=[]
X1= range(Xmin,Xmax, StepSize)
Y1= range(Ymin,Ymax, StepSize)
for ys in Y1:
    for xs in X1:   
        pts=pdmtn[(pdmtn.X>=xs) &(pdmtn.X<=xs+StepSize)&(pdmtn.Y>=ys) &(pdmtn.Y<=ys+StepSize)].Z.unique()
        if len(pts)==0:
            ptsz=200
        else:
            ptsz=np.mean(pts)       
        vals.append(ptsz)

vals = np.array(vals)
#where_are_NaNs = np.isnan(vals)
#vals[where_are_NaNs] = 200
zz=vals.reshape(len(Y1), len(X1))
[xx,yy]=np.meshgrid(np.arange(Xmin,Xmax,StepSize),np.arange(Ymin,Ymax,StepSize))


cmap = plt.get_cmap('PiYG')
# Plot the height map
plt.figure(figsize=(6,5),dpi=150)
plt.pcolormesh(xx,yy,zz,cmap=cmap,shading='auto')
print("got here 0")

plt.colorbar(label='Vertical Distance for Lab to Surface (m)')
print("got here 1")
CS=plt.contour(xx,yy,zz,colors='black')

print("got here 2")
plt.clabel(CS, inline=1, fontsize=10,fmt='%1.0f')

print("got here 2.5")
plt.ylim(Ymin,Ymax)
plt.xlim(Xmin,Xmax)
plt.xlabel("X distance (m)")
plt.ylabel("Y distance (m)")
plt.title("Depth Profile Map")
plt.plot([0],[0],'x',label='Location of LSC',color='red')
plt.legend(loc='lower right')


print("got here 2.75")

#Find the thetas and phis and distance through rock for each grid square
phioffset=0                      # Orientation of detector relative to map - you need to figure this out.
rho=((xx+eps)**2+(yy+eps)**2)**0.5             # cylindrical rho coordinate
theta=np.arctan(rho/(zz+eps))           # spherical theta coordinate (0 = downgoing)
phi = np.arctan((yy+eps)/(xx+eps)) + phioffset # spherical phi coordinate
distancetodetector = np.sqrt((zz+eps)**2+(xx+eps)**2+(yy+eps)**2)


#Plot the distance from the surface to the lab at each place in XY

plt.figure(figsize=(6,5),dpi=150)
plt.scatter(xx,yy,StepSize,distancetodetector)

plt.colorbar(label='distance to detector (m)')
CS=plt.contour(xx,yy,distancetodetector,colors='black')
plt.clabel(CS, inline=1, fontsize=10,fmt='%1.0f')


plt.xlabel("X distance (m)")
plt.ylabel("Y distance (m)")
plt.title("Rock Distance Map")
plt.plot([0],[0],'x',label='Location of LSC',color='red')
plt.legend(loc='lower right')


# continuously interpolate rock length and depth functions so we can sample at any X,Y
rocklength=scipy.interpolate.RectBivariateSpline(xx[0],yy[:,1],np.transpose(distancetodetector),s=0,kx=3, ky=3)
depth     =scipy.interpolate.RectBivariateSpline(xx[0],yy[:,1],np.transpose(zz),s=0,kx=3, ky=3)


# a handy function to go from surface position to angles, depth, distance
def GetMuonInfo(startMuon):
    depthMuon   = np.round(depth(*startMuon)[0][0],2)
    distMuon    = np.round(rocklength(*startMuon)[0][0],2)
    thMuon      = np.round(np.arctan((startMuon[0]**2+startMuon[1]**2)**0.5/depthMuon),2)
    phiMuon     = np.round(np.arctan(startMuon[1]/startMuon[0]),2)
    return thMuon,phiMuon,depthMuon,distMuon


#This function does the business of calling PROPOSAL.
# It is configured in config.json to propagate the muon through a giant block of "standard rock"

def PropagateMuons(MuonEnergyToSimulate,DistToDetector,NumberToRun=1000):
    mu_def = pp.particle.MuMinusDef()
    prop = pp.Propagator(
        particle_def=mu_def,
        config_file=".config.json"   #in the PROPOSAL resources directory
    )
    print('A')
    mu = pp.particle.DynamicData(mu_def.particle_type)
    print('B')
    mu.energy = MuonEnergyToSimulate
    mu.direction = pp.Vector3D(0, 0, -1)
    print('C')
    mu_position = []
    mu_energy = []

    for i in range(NumberToRun):
        sec = prop.propagate(mu,DistToDetector)
        slop=100
        if(np.abs(sec.position[-1].magnitude()-DistToDetector)<slop):
            mu_energy.append(sec.energy[-1])
            mu_position.append(sec.position[-1].magnitude()-DistToDetector)
    return mu_energy,mu_position


#Decide on a place to start a muon from and get its angles and distances:
startMuon   = (50,-50)
thMuon,phiMuon,depthMuon,distMuon=GetMuonInfo(startMuon)
print("theta, phi, depth, distance to det: ", thMuon,phiMuon,depthMuon,distMuon)

thMuon=0.55 
phiMuon=1.55 
depthMuon=647.73 
distMuon=760.42


NumToRun=1000
# Run a few energies and make a plot
#plt.figure(figsize=(5,5),dpi=150)
print("got here 3")

E=1000*GeV

FinalMuons,FinalPos=np.array(PropagateMuons(E,distMuon*m_to_cm,NumberToRun=1000),dtype=object)
FracSurviving=round(len(FinalMuons)/NumToRun,1)

plt.hist(FinalMuons/E,bins=np.arange(0,1.02,0.02),label=str(E/1e6)+ " TeV muons, "+str(FracSurviving*100)+"% survive",alpha=0.5,color='DarkRed')

print("got here 4")
E=5000*GeV
FinalMuons,FinalPos=np.array(PropagateMuons(E,distMuon*m_to_cm,NumberToRun=1000),dtype=object)
FracSurviving=round(len(FinalMuons)/NumToRun,1)
plt.hist(FinalMuons/E,bins=np.arange(0,1.02,0.02),label=str(E/1e6)+ " TeV muons, "+str(FracSurviving*100)+"% survive",alpha=0.5,color='DarkBlue')

E=25000*GeV
FinalMuons,FinalPos=np.array(PropagateMuons(E,distMuon*m_to_cm,NumberToRun=1000),dtype=object)
FracSurviving=round(len(FinalMuons)/NumToRun,1)
plt.hist(FinalMuons/E,bins=np.arange(0,1.02,0.02),label=str(E/1e6)+ " TeV muons, "+str(FracSurviving*100)+"% survive",alpha=0.5,color='Black')


plt.legend(loc='upper right')
plt.xlabel("Fraction of original energy")
plt.title(r"$\theta,\phi$ = "+str(thMuon)+"," +str(phiMuon))
plt.ylim(0,plt.gca().get_ylim()[1]*1.2)
plt.savefig('Survival.png',dpi=250,bbox_inches='tight')
plt.show()
