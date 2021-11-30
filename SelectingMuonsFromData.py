#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 00:03:35 2021

@author: rogerslc
"""
import sys
import glob

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d as m3d
import math
from matplotlib import image
from matplotlib import pyplot

import os.path
from os import path

from scipy.signal import convolve as scipy_convolve

def randrange(n, vmin, vmax):
    return (vmax-vmin)*np.random.rand(n) + vmin


def Ry (R, phi, theta):
    return (R) * np.cos(phi)

def Rx (R, phi, theta):
    return (R) * np.sin(phi) * np.cos(theta)

def Rz (R, phi, theta):
    return (R) * np.sin(phi) * np.sin(theta)

def Distance (P1,P2,P3,U1,U2,U3):
    S_1=P2*U3-P3*U2
    S_2=P3*U1-P1*U3
    S_3=P1*U2-P2*U1
    Dist=np.sqrt(S_1**2+S_2**2+S_3**2)
    return Dist



low_lim=2
MaxTimefromMuon=.01

# path to save the file                                                                                                                                                           
savelocat=sys.argv[1]


#run numbers
start_run = int(sys.argv[2])
end_run   = int(sys.argv[3])


FidCut=178*2**-.5 #fiducial cut for x and y in mm                                                                                                                            
depth=510

numb=0
vars=np.arange(-205,205,100)
R=np.arange(-3500,3500,400)
width=430
height=430

angle_step=1

radius =180
allowedspread=35 #mm on either side of the track                                                                                                                             
allowedT=1
tracklenmin=120 #mm minimum accepted length of track within detector                                                                                                         
dedxmin=.017 #MeV/mm maximum dE/dx, determined from MC                                                                                                                       


diag_lenxy = int(np.ceil(np.sqrt(width * width + height * height)))   # max_dist                                                                                          
diag_lenzy = int(np.ceil(np.sqrt(depth * depth + height * height)))   # max_dist                                                                                          

rhos = np.linspace(-diag_lenxy, diag_lenxy, diag_lenxy * 2)
rhozed = np.linspace(-diag_lenzy, diag_lenzy, diag_lenzy * 2)


thetas = np.deg2rad(np.arange(-90, 90, angle_step))
num_thetas = len(thetas)


#For convolution with gausians in Fourier space                                                                                                                              
#First a 1-D  Gaussian                                                                                                                                                       
t = np.linspace(-10, 10, 30)
bump = np.exp(-0.1*t**2)
bump /= np.trapz(bump) # normalize the integral to 1                                                                                                                         

# make a 2-D kernel out of it                                                                                                                                                
kernel = bump[:, np.newaxis] * bump[np.newaxis, :]

#locatformat='/analysis/{}/hdf5/prod/v1.2.0/20191122/cdst/trigger2'


#actual data from Next-White                                                                                                                                                                                                                                                                
locatformat='/analysis/{}/hdf5/prod/v1.2.0/20191122/cdst/trigger2/'



RunV = [6971,6972,6973,6976,6977,
             6988, 6990,6991,6992,6993,6994,6995,6996,6997,7006,7007,
             7010,7011,7012,7013,7014,7015,7016,7023,
             7036,7037,7038,7039,
             7042,7044,7046,7047,7048,7049,
             7050,7065,7066,7067,7069,
             7070,7071,7072,7073,7074,7075,7079,
             7084, 7085,7106,7107,7108,
             7112,7113,7114,7115,7116,7117,7118,
             7121,7122,7123,7124,7142,7143,7144,7145,
             7175,7204,7205,7210,7211,7212,7214,7234,7247,7517,7518,7519,
            7520,7540,7541,7542,7543,7544,7545,7546,7547,7548,7549,
            7550,7551,7555,7556,7557,
            7563,7564,7565,7566,7567,7568,7569,
            7591,7592,7593,7594,7595,7596,7597,7598,7599,
            7600,7601,7602,7603,7604,7605,7608,7609,
            7610,7611,7612, 7751,7752,7753,7755,7757,7758,7759,
            7761,7762,7763,7764,7765,7766,7767,
            7770,7773,7774,7775,7776,
            7782,7801,7813,7814,7817,7818,
            7821,7824,7827,7828,
            7830,7836,7837,7838,7839,
            7840,7841,7842,7849,
            7850,7851,7865,7866,7867,7869,
            7870,7871,7872,7877,7878,7879,
            7880,7881,7882,7883,7884,7885,7886, # March 10: PMT issue                                                  
            7911,7913,7914,7915,7916,7917,7919,
            7920,7921,7922,7923,7924,
            7931,7932,7933,7934,7935,7936,7937,7938,7939,
            7940,7941,7942,7943,7944,
            7945,7946,7947,7948,7949,
            7951,7952,7953,7954,7955,
            7973,7974,7975,7976,7977,7978,7979,
            7981,7982,7983,7984,7985,7986,7987,7988,7989,
            7990,7991,7992,7993,7994,7995,7996,7997,7998,7999,
            8003,8004,8006,8009,
            8010,8011,8012,8013,8014,8015,8016,8017,8018,
            8020,8021,8022,8023,8025,8026,8027,8028,8029,
            8030,8031,8032,8033,8034,
            8053,8054,8055,8056]






num =0
probnum=0
ttlevts=0

#for i in range(start_run, end_run):
for i in RunV[start_run:end_run]:
    run_num=str(i).zfill(4)
    locat=locatformat.format(str(i))
    fileformat= locat+'cdst_{}_v1.2.0_trigger2_bg.h5'''
    files = fileformat.format(run_num)
    print('Starting on file ',files)
    
    try:
        events_data=pd.read_hdf(files,key='/CHITS/highTh')
        line_data=pd.read_hdf(files,key='/Tracking/Tracks')
        

    except FileNotFoundError:
        print('file not there')
        continue    
    except KeyError:
        print('KeyError')
        continue
    

    outside_fid=line_data[(line_data.r_max > 197) | (line_data['z_min'] < 20) | (line_data['z_max'] > 490) & line_data.z_max < depth]

    ttlevts+=len(outside_fid['event'].unique())
    
    
    dedxall=[]
    pcline=[]
    linlen=[]
    munrg=[]
    muvts=[]
    mu_theta=[]
    mu_phi=[]
    xint=[]
    yint=[]
    zint=[]
    Pass=[]
    FidCross=[]
    truth_theta=[]
    truth_phi=[]
    truth_per=[]
    truth_dedx=[]
    truth_len=[]
    Weight=[]
    zvar=[]

    costhetas=[]
    costhetasall=[]

    for evts in events_data.event.unique():
        #print(run_num,evts,MC_particles[(MC_particles.event_id==evts)&(MC_particles.particle_id==1)])

        #print("it should have passed",run_num)

        #fiducial cuts                                                                                                                                                       
        events_dataplt=events_data[events_data['event']==evts][(events_data['Ec']>0) & (events_data.Z<=depth) & (events_data.Z>=20) ]
        events_dataplt=events_dataplt[(events_dataplt.X<FidCut) & (events_dataplt.X>-FidCut)]
        events_dataplt=events_dataplt[(events_dataplt.Y<FidCut) & (events_dataplt.Y>-FidCut)]


        ttl_nrg=sum(events_dataplt['Ec'])
        #if ttl_nrg==0:
            #print ("no energy deposited",ttl_nrg, run_num,evts)
            #continue
        max_nrg=np.max(events_data[events_data['event']==evts]['Ec'])


        #use this if doing full muon cuts 
        if ttl_nrg < low_lim:
            continue
        EmxEtl=max_nrg/ttl_nrg


        #use this to remove sparks
        #sparkcheck=(np.var(events_dataplt.X)-np.var(events_dataplt.Y))/((np.var(events_dataplt.X)+np.var(events_dataplt.Y))*.5)
        #if (np.var(events_dataplt.Z) <10) & (sparkcheck<1):
        #    print(evts,np.var(events_dataplt.Z),"spark spark")
        #    continue
        
        #now for accumlators to find best rhos and thetas
        accumulatorXY = np.zeros((2 * diag_lenxy, num_thetas), dtype=np.uint64)
        accumulatorZY = np.zeros((2 * diag_lenzy, num_thetas), dtype=np.uint64)
        accumulatorXZ = np.zeros((2 * diag_lenzy, num_thetas), dtype=np.uint64)


        for x,y,z in events_dataplt[['X','Y','Z']].values:
            posx=x
            posy=y
            posz=z*.97-depth/2 #to change from microsec to mm, the minus is just for a mid axis for hough trans                                                              
            tval=0
            for t in thetas:
                rhoXY=round(posx * np.cos(t) + posy * np.sin(t),0) + diag_lenxy #XY plane                                                                                    
                rhoZY=round(posz * np.cos(t) + posy * np.sin(t),0) + diag_lenzy #ZY plane   
                rhoXZ=round(posx * np.cos(t) + posz * np.sin(t),0) + diag_lenzy #XZ plane          
                accumulatorXY[int(rhoXY), tval] += 1
                accumulatorZY[int(rhoZY), tval] += 1
                accumulatorXZ[int(rhoXZ), tval] += 1
                tval+=1


        # Convolution: scipy's direct convolution mode spreads out NaNs                                                               
        conv_accuXY = scipy_convolve(accumulatorXY, kernel, mode='same', method='direct')
        conv_accuZY = scipy_convolve(accumulatorZY, kernel, mode='same', method='direct')
        conv_accuXZ = scipy_convolve(accumulatorXZ, kernel, mode='same', method='direct')

        locatXY = np.unravel_index(np.argmax(conv_accuXY, axis=None), accumulatorXY.shape)
        locatZY = np.unravel_index(np.argmax(conv_accuZY, axis=None), accumulatorZY.shape)
        locatXZ = np.unravel_index(np.argmax(conv_accuXZ, axis=None), accumulatorXZ.shape)


        rhoXY = rhos[locatXY[0]]
        thetaXY = thetas[locatXY[1]]
        rhoZY=rhozed[locatZY[0]]
        thetaZY=thetas[locatZY[1]]
        rhoXZ=rhozed[locatXZ[0]]
        thetaXZ=thetas[locatXZ[1]]

        
        if (thetaZY ==0) or (thetaXY==0) or (thetaXZ==0):
            print ("need to implement catch for thetas =0",run_num,evts)
            continue

        
        alpha1=np.arctan2(np.sqrt(np.tan(thetaXY)**2+np.tan(thetaZY)**2),1)
        beta1=np.arctan2(-np.tan(thetaZY),-np.tan(thetaXY))

        alpha2=np.arctan2(np.sqrt(np.tan(thetaZY)**2+np.tan(thetaXZ)**2*np.tan(thetaZY)**2),1)
        alpha3=np.arctan2(np.sqrt(np.tan(thetaXY)**2+np.tan(thetaXY)**2/np.tan(thetaXZ)**2),1)     



        if thetaXY*thetaXZ >0:
            beta2=thetaXZ+np.pi/2
            beta3=thetaXZ+np.pi/2

        if thetaXY*thetaXZ <0:
            beta2=thetaXZ-np.pi/2
            beta3=thetaXZ-np.pi/2    


        #there has to be a cleaner way to do this        
        if alpha1>np.pi/2:
            beta1+=np.pi
            alpha1=np.pi-alpha1
        if beta1<0:
            beta1+=2*np.pi
        if beta1>2*np.pi:
            beta1-=2*np.pi    

        if alpha2>np.pi/2:
            beta2+=np.pi
            alpha2=np.pi-alpha2
        if beta2<0:
            beta2+=2*np.pi
        if beta2>2*np.pi:
            beta2-=2*np.pi 

        if alpha3>np.pi/2:
            beta3+=np.pi
            alpha3=np.pi-alpha3
        if beta3<0:
            beta3+=2*np.pi
        if beta3>2*np.pi:
            beta3-=2*np.pi     

        
        #intersections for each set where they share a common axis
        x0=rhoXY/np.cos(thetaXY)
        z0=rhoZY/np.cos(thetaZY) 

        x1=rhoXZ/np.cos(thetaXZ)
        y1=rhoZY/np.sin(thetaZY)

        y2=rhoXY/np.sin(thetaXY)
        z2=rhoXZ/np.sin(thetaXZ)

        weightXY=conv_accuXY[int(diag_lenxy+rhoXY), int(90+np.rad2deg(thetaXY))]
        weightXZ=conv_accuXZ[int(diag_lenzy+rhoXZ), int(90+np.rad2deg(thetaXZ))]
        weightZY=conv_accuZY[int(diag_lenzy+rhoZY), int(90+np.rad2deg(thetaZY))]
        
        SHAPE0=np.shape(conv_accuXY)[0]
        SHAPE1=np.shape(conv_accuXZ)[0]
        SHAPETHET=np.shape(conv_accuXY)[1]
        
        #finding set of points from each pair that projects onto 3rd plane's hough space for best 3d line definition

        if (np.abs(np.rad2deg(beta1)==180)) or (alpha1==0):
            weightXZ1=0
        else:
            Rpz=-z0/(np.sin(alpha1)*np.sin(beta1))
            Rpx=-x0/(np.sin(alpha1)*np.cos(beta1))
            xp0=Rx(Rpz,alpha1,beta1)+x0
            zp0=Rz(Rpx,alpha1,beta1)+z0
            thetxz=np.arctan(xp0/(zp0))
            roxz=xp0*np.cos(thetxz)

            if (int(diag_lenzy+roxz) > SHAPE1) or (int(90+np.rad2deg(thetxz))>SHAPETHET) or (int(diag_lenzy+roxz) <0) or (int(90+np.rad2deg(thetxz))<0):
                weightXZ1=0
            else:
                weightXZ1=conv_accuXZ[int(diag_lenzy+roxz), int(90+np.rad2deg(thetxz))] #put in roxz                                                                                         



        if (np.abs(np.rad2deg(beta2)==180)) or (alpha2==0):
            weightXY2=0
        else:
            Rpy=-y1/(np.cos(alpha2))
            Rpx=-x1/(np.sin(alpha2)*np.cos(beta2))
            xp1=Rx(Rpy,alpha1,beta1)+x1
            yp1=Ry(Rpx,alpha1,beta1)+y1
            thetxy=np.arctan(xp1/yp1)
            roxy=xp1*np.cos(thetxy)

            if (int(diag_lenxy+roxy) > SHAPE0) or (int(90+np.rad2deg(thetxy))>SHAPETHET) or (int(diag_lenxy+roxy) <0) or (int(90+np.rad2deg(thetxy))<0):
                weightXY2=0
            else:
                weightXY2=conv_accuXY[int(diag_lenxy+roxy), int(90+np.rad2deg(thetxy))]


        if (np.abs(np.rad2deg(beta3)==180)) or (alpha3==0):
            weightZY3=0
        else:
            Rpz=-z2/(np.sin(alpha3)*np.sin(beta3))
            Rpy=-y2/(np.cos(alpha3))
            yp2=Ry(Rpz,alpha3,beta3)+y2
            zp2=Rz(Rpy,alpha3,beta3)+z2
            thetzy=np.arctan((zp2)/yp2)
            rozy=zp2*np.cos(thetzy)
            if (int(diag_lenzy+rozy) > SHAPE1) or (int(90+np.rad2deg(thetzy))>SHAPETHET) or (int(diag_lenzy+rozy) <0) or (int(90+np.rad2deg(thetzy))<0):
                weightZY3=0
            else:
                weightZY3=conv_accuZY[int(diag_lenzy+rozy), int(90+np.rad2deg(thetzy))]

        ProdW1=weightXY*weightXZ1*weightZY    
        ProdW2=weightXY2*weightXZ*weightZY    
        ProdW3=weightXY*weightXZ*weightZY3

        WTs=[ProdW1,ProdW2,ProdW3]
        As=[alpha1,alpha2,alpha3]
        Bs=[beta1,beta2,beta3]
        xs=[x0,x1,0]
        ys=[0,y1,y2]
        zs=[z0,0,z2]


        #selecting the best set
        LOC=np.argmax(WTs)
        alpha=As[LOC]
        beta=Bs[LOC]


        Xint=xs[LOC]
        Yint=ys[LOC]    
        Zint=zs[LOC]+depth/2

        bestx=Rx(R,alpha,beta)+Xint
        bestz=Rz(R,alpha,beta)+Zint
        besty=Ry(R,alpha,beta)+Yint

        
        events_dataplt['cross']=Distance((events_dataplt.X-Xint),(events_dataplt.Z-Zint),(events_dataplt.Y-Yint),Rx(1,alpha,beta),Rz(1,alpha,beta),Ry(1,alpha,beta))



  
        linelike=events_dataplt[(np.abs(events_dataplt.cross)<=allowedspread)]
        linelen=np.sqrt((np.min(linelike.X)-np.max(linelike.X))**2+(np.min(linelike.Y)-np.max(linelike.Y))**2+(np.min(linelike.Z)-np.max(linelike.Z))**2)
        nrgonline=sum(linelike['Ec'])/ttl_nrg
        nrgovrlen=sum(linelike['Ec'])/linelen        
        

        #collecting all the data for saving at end
        if (nrgovrlen<.015) & (linelen>73.5) & (nrgonline>.79):
            dedxall.append(nrgovrlen)  
            pcline.append(nrgonline)
            linlen.append(linelen)
            munrg.append(ttl_nrg)
            muvts.append(evts)
            mu_theta.append(beta)
            mu_phi.append(alpha)
            xint.append(Xint)
            yint.append(Yint)
            zint.append(Zint)
            zvar.append(np.var(events_dataplt.Z))

        #2d plots
        '''xdataplot=events_dataplt.X
        ydataplot=events_dataplt.Z
        fig = plt.figure(figsize=(10,8))
        ax = fig.add_subplot()
        ax.set_xlabel('X')
        ax.plot(bestx,bestz,color='magenta',linewidth=3,label='best')
        ax.plot([FidCut,FidCut],[-depth,depth],color='k')
        ax.plot([-FidCut,-FidCut],[-depth,depth],color='k')
        ax.plot([-width/2,width/2],[490,490],color='k')
        ax.plot([-width/2,width/2],[20,20],color='k')
        ax.set_ylabel('Z')
        ax.legend(fontsize=20)
        ax.set_xlim(-215, 215)
        ax.set_ylim(0, 500)
        colmap = cm.ScalarMappable(cmap=cm.jet)
        ax.scatter(xdataplot,ydataplot,c=cm.jet(events_dataplt['Ec']/np.max(events_dataplt['Ec'])))


        xdataplot=events_dataplt.X
        ydataplot=events_dataplt.Y
        fig = plt.figure(figsize=(10,8))
        ax = fig.add_subplot()
        ax.set_xlabel('X')
        ax.plot(bestx,besty,color='magenta',linewidth=3,label='best')
        ax.plot([FidCut,FidCut],[-depth,depth],color='k')
        ax.plot([-FidCut,-FidCut],[-depth,depth],color='k')
        ax.plot([-width/2,width/2],[-FidCut,-FidCut],color='k')
        ax.plot([-width/2,width/2],[FidCut,FidCut],color='k')
        ax.set_ylabel('Y')
        ax.legend(fontsize=20)
        ax.set_xlim(-215, 215)
        ax.set_ylim(-215, 215)
        colmap = cm.ScalarMappable(cmap=cm.jet)
        ax.scatter(xdataplot,ydataplot,c=cm.jet(events_dataplt['Ec']/np.max(events_dataplt['Ec'])))


        xdataplot=events_dataplt.Z
        ydataplot=events_dataplt.Y
        fig = plt.figure(figsize=(10,8))
        ax = fig.add_subplot()
        ax.set_xlabel('Z')
        ax.plot(bestz,besty,color='magenta',linewidth=3,label='best')
        ax.plot([0,depth],[-FidCut,-FidCut],color='k')
        ax.plot([0,depth],[FidCut,FidCut],color='k')
        ax.plot([20,20],[-FidCut,FidCut],color='k')
        ax.plot([490,490],[-FidCut,+FidCut],color='k')
        ax.set_ylabel('Y')
        ax.set_xlim(0, 500)
        ax.set_ylim(-215, 215)
        colmap = cm.ScalarMappable(cmap=cm.jet)
        ax.legend(fontsize=20)
        ax.scatter(xdataplot,ydataplot,c=cm.jet(events_dataplt['Ec']/np.max(events_dataplt['Ec'])))'''
        
        #hough transform convolved plots                                                                                                                                       
        '''fig = plt.figure(figsize=(9,9))
        ax = fig.add_subplot(111)
        ax.imshow(np.log(1+conv_accuXY), cmap='jet',extent=[np.rad2deg(thetas[0]), np.rad2deg(thetas[-1]), diag_lenxy, -diag_lenxy])
        ax.plot([-90,90],[rhoXY,rhoXY],color='white')
        ax.plot([np.rad2deg(thetaXY),np.rad2deg(thetaXY)],[-200,200],color='white')
        ax.scatter(np.rad2deg(thetxy),roxy,color='white')
        ax.set_aspect('equal', adjustable='box')
        ax.set_title('Hough trans_conv XY')
        ax.set_xlabel('theta (degrees)')
        ax.set_ylabel('Distance')
        plt.ylim(200,-200)
        #fig.savefig(savelocat+'hough_convolved'+run_num+'_'+str(evts)+'XY.png')
                                                                                                                                       \

        fig = plt.figure(figsize=(9,9))
        ax = fig.add_subplot(111)
        ax.imshow(np.log(1+conv_accuZY), cmap='jet',extent=[np.rad2deg(thetas[0]), np.rad2deg(thetas[-1]), diag_lenzy, -diag_lenzy])
        ax.plot([-90,90],[rhoZY,rhoZY],color='white')
        ax.scatter(np.rad2deg(thetzy),rozy,color='white')
        ax.plot([np.rad2deg(thetaZY),np.rad2deg(thetaZY)],[-width/2,width/2],color='white')
        ax.set_aspect('equal', adjustable='box')
        ax.set_title('Hough trans_conv ZY')
        ax.set_xlabel('theta (degrees)')
        #plt.xlim(-5,5)
        plt.ylim(200,-200)
        ax.set_ylabel('Distance')
        #fig.savefig(savelocat+'hough_convolved'+run_num+'_'+str(evts)+'ZY.png')

        fig = plt.figure(figsize=(9,9))
        ax = fig.add_subplot(111)
        ax.imshow(np.log(1+conv_accuXZ), cmap='jet',extent=[np.rad2deg(thetas[0]), np.rad2deg(thetas[-1]), diag_lenzy, -diag_lenzy])
        ax.plot([-90,90],[rhoXZ,rhoXZ],color='white')
        ax.scatter(np.rad2deg(thetxz),roxz,color='white')
        ax.plot([np.rad2deg(thetaXZ),np.rad2deg(thetaXZ)],[-width/2,width/2],color='white')
        ax.set_aspect('equal', adjustable='box')
        ax.set_title('Hough trans_conv XZ')
        ax.set_xlabel('theta (degrees)')
        #plt.xlim(-5,5)
        plt.ylim(200,-200)
        ax.set_ylabel('Distance')
        #fig.savefig(savelocat+'hough_convolved'+run_num+'_'+str(evts)+'ZY.png')'''
        
        #3d plot                                                                                                                                                             
        '''xdataplot=events_dataplt.X
        ydataplot=events_dataplt.Z
        zdataplot=events_dataplt.Y

        fig = plt.figure(figsize=(10,8))
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlabel('X')
        ax.plot(bestx,bestz,besty,color='magenta',linewidth=3,label='best')
        ax.set_ylabel('Z')
        ax.set_zlabel('Y')
        ax.set_xlim(-215, 215)
        ax.set_zlim(-215, 215)
        ax.set_ylim(0, 500)
        colmap = cm.ScalarMappable(cmap=cm.jet)
        ax.scatter(xdataplot,ydataplot,zdataplot,c=cm.jet(events_dataplt['Ec']/np.max(events_dataplt['Ec'])))
        #cb = fig.colorbar(colmap)                                                                                                                                           
        #fig.savefig(savelocat+'3Dplot'+run_num+'_'+str(evts)+'.png')'''



    #save for each run of events, things tend to run out of memory and crash
    data_out4 = savelocat+'/run'+run_num+'alleventscutinfo_highTH.h5'
    pd.DataFrame({'dEdx':dedxall,
                  'perconline':pcline,
                  'linelength':linlen,
                  'muenergy':munrg,
                  'eventnum':muvts,
                  'beta':mu_theta,
                  'alpha':mu_phi,
                  'zvariance':zvar,
                  'xintercept':xint,
                  'yintercept':yint,
                  'zintercept':zint}).to_hdf(data_out4,'dEdx')
        
    
    


