#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 16:50:56 2021

@author: leohoinaski
"""
#import netCDF4 as nc4
import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import MultiLineString
from shapely.ops import polygonize
import datetime
import numpy.matlib

#%% Gridding and populatingGrid functions

def gridding(lon,lat):
    xv, yv = np.meshgrid(lon, lat)
    hlines = [((x1, yi), (x2, yi)) for x1, x2 in zip(lon[:-1], lon[1:]) for yi in lat]
    vlines = [((xi, y1), (xi, y2)) for y1, y2 in zip(lat[:-1], lat[1:]) for xi in lon]
    grids = list(polygonize(MultiLineString(hlines + vlines)))
    grids = gpd.GeoDataFrame(grids) 
    grids.columns =['geometry'] 
    grids['geometry'] = grids['geometry']
    grids.crs = "EPSG:4326"  
    grids['X'] = grids.geometry.centroid.x
    grids['Y'] = grids.geometry.centroid.y
    xX = np.array(grids['X']).reshape((lon.shape[0]-1,lat.shape[0]-1)).transpose()
    yY = np.array(grids['Y']).reshape((lon.shape[0]-1,lat.shape[0]-1)).transpose()
    return grids,xv,yv,xX,yY

def populatingGrid(dataEmiss,center,xX,yY,xv,yv):   
    data = np.zeros([1,dataEmiss.shape[1],np.size(yv,0)-1, np.size(xv,1)-1])
    xcenter = center.geometry.centroid.x
    ycenter = center.geometry.centroid.y
   
    for ii in range(0,dataEmiss.shape[0]):
        dist = ((xcenter[ii]-xX)**2 + (ycenter[ii]-yY)**2)**(1/2)
        mindist = np.where(dist == np.amin(dist))
        print('Fire number '+str(ii)+' from '+str(dataEmiss.shape[0]))
        for kk in range (0,dataEmiss.shape[1]):
            data[0,kk,mindist[0][0],mindist[1][0]]= data[0,kk,mindist[0][0],mindist[1][0]]+dataEmiss.iloc[ii,kk]     
    return data

def populatingGridMatBB(bb,dataEmissBB,center,xX,yY,year,month):
    
    # Extracting year, month and day frmo FINN file
    bb['datetimeBB']=''
    for ii in range(0,bb.shape[0]):
        dateT= [str(year)+'-'+str(bb.iloc[ii,1])]
        datetimebb = datetime.datetime.strptime(dateT[0],'%Y-%j')
        bb['datetimeBB'][ii] = datetimebb.strftime('%Y-%m-%d %H:%M:%S')   
    bb['year'] = pd.DatetimeIndex(bb['datetimeBB']).year
    bb['month'] = pd.DatetimeIndex(bb['datetimeBB']).month  
    bb['day'] = pd.DatetimeIndex(bb['datetimeBB']).day
    
    # Creating a perfect day array
    startDate = datetime.datetime(year, month, 1, 0, 0)
    endDate = datetime.datetime(year, month+1, 1, 0, 0)
    datePfct = np.arange(np.datetime64(startDate),np.datetime64(endDate),3600000000)
    datePfct = pd.DataFrame(datePfct)
    datePfct['year'] = pd.DatetimeIndex(datePfct.iloc[:,0]).year
    datePfct['month'] = pd.DatetimeIndex(datePfct.iloc[:,0]).month  
    datePfct['day'] = pd.DatetimeIndex(datePfct.iloc[:,0]).day
    # Creating a new array with perfect date, var and cells
    #dataTempo = np.zeros([np.size(datePfct,0),dataEmissBB.shape[1],np.size(yY,0), np.size(xX,1)])
    dataTempo1 = np.zeros([384,dataEmissBB.shape[1],np.size(yY,0), np.size(xX,1)])
    dataTempo2 = np.zeros([np.size(datePfct,0)-384,dataEmissBB.shape[1],np.size(yY,0), np.size(xX,1)])
    # extracting centroids from BB cases
    xcenter = center.geometry.centroid.x
    ycenter = center.geometry.centroid.y
    
    for jj in range(0,datePfct.shape[0]):
        bbUsed = bb[(bb['month']==datePfct['month'][jj]) & (bb['day']==datePfct['day'][jj])]
        xcenterUsed = xcenter[(bb['month']==datePfct['month'][jj]) & (bb['day']==datePfct['day'][jj])]
        ycenterUsed = ycenter[(bb['month']==datePfct['month'][jj]) & (bb['day']==datePfct['day'][jj])]      
        if jj<384:
            if bbUsed.shape[0]>0:
                for ii in range(0,bbUsed.shape[0]):
                    dist = ((xcenterUsed.iloc[ii]-xX)**2 + (ycenterUsed.iloc[ii]-yY)**2)**(1/2)
                    mindist = np.where(dist == np.amin(dist))
                    print('cell number = ' + str(ii))
                    for kk in range (0,dataEmissBB.shape[1]):
                        dataTempo1[jj,kk,mindist[0][0],mindist[1][0]]= dataTempo1[jj,kk,mindist[0][0],mindist[1][0]]+dataEmissBB.iloc[ii,kk]
        else:
            if bbUsed.shape[0]>0:
                for ii in range(0,bbUsed.shape[0]):
                    dist = ((xcenterUsed.iloc[ii]-xX)**2 + (ycenterUsed.iloc[ii]-yY)**2)**(1/2)
                    mindist = np.where(dist == np.amin(dist))
                    print('cell number = ' + str(ii))
                    for kk in range (0,dataEmissBB.shape[1]):
                        dataTempo2[jj-384,kk,mindist[0][0],mindist[1][0]]= dataTempo2[jj-384,kk,mindist[0][0],mindist[1][0]]+dataEmissBB.iloc[ii,kk]
    datePfct1 = datePfct.iloc[0:384,:]
    datePfct2 = datePfct.iloc[384:,:]
    datePfct2 = datePfct2.reset_index(drop=True)
    
    return dataTempo1,datePfct1,dataTempo2,datePfct2 