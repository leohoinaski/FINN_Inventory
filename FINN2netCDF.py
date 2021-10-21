#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------
                             FINN2netCDF.py
                             
                             
This is the main script to convert daily fire emissions from 
Fire INventory from NCAR (FINN) (https://www.acom.ucar.edu/Data/fire/).

Inputs: 
    
    rootPath: Path to functions and BRAVESdatabase_main.py
    
    outPath: Path to folder with roadLenght....csv
    
    lati: Initial latitude (lower-left)
    
    latf: Final latitude (upper-right)
    
    loni: Initial longitude (lower-left)
    
    lonf: Final longitude (upper-right)
    
    deltaX: Grid resolution/spacing in x direction
    
    deltaY: Grig resolution/spacing in y direction
    
    year: Base-year for your simulation
            
    fileId = identification of your output files
    
    cutORNotFINN = cut or not FINN origial file using Brazilian shapefile. If you
        have already the intermediate file, you do not need to run (0)
    
    conv = conversion factor from mol/day to mol/s
    
    speciation = speciation profiles. You could use GEOS-CHEM, MOZART or SAPRC99
        you can find more details at:https://www.acom.ucar.edu/Data/fire/
    
Outputs:
    
    'FINNannualEmiss_'+speciation+'_'+fileId+'_'+str(deltaX)+'x'+str(deltaY)+'_'+str(year)+'.nc'
    
    
External functions:
    cutFINNstatesUtil, FINNgridding, netCDFcreator
    
Last update = 21/10/2021

Author: Leonardo Hoinaski - leonardo.hoinaski@ufsc.br

---------------------------------------------------------------
"""

#import netCDF4 as nc4
import geopandas as gpd
import pandas as pd
import os 
import numpy as np
import datetime
from netCDFcreator import createNETCDFtemporalFINN
import numpy.matlib
from shapely.geometry import Polygon
from FINNgridding import gridding, populatingGrid, populatingGridMatBB
from cutFINNstatesUtil import cutFINN

#%%============================= INPUTS =========================================

rootPath = '/media/leohoinaski/HDD/FINN_inventory'

outPath = rootPath +"/Outputs"

#-------------------------Setting grid resolution------------------------------

# Users can change the domain and resolution here.
lati =-30 #lati = int(round(bound.miny)) # Initial latitude (lower-left)

latf = -24 #latf = int(round(bound.maxy)) # Final latitude (upper-right)

loni = -54 #loni = int(round(bound.minx)) # Initial longitude (lower-left)

lonf = -47 #lonf = int(round(bound.maxx)) # Final longitude (upper-right)

deltaX = 0.01 # Grid resolution/spacing in x direction

deltaY = 0.01 # Grig resolution/spacing in y direction

fileId = 'SC' # Code to identify your output files

prefix = str(deltaX)+'x'+str(deltaY) # grid definition identification

year = 2019

month = 1

cutORnotFINN = 0

convt = 3600

speciation  = 'GEOS-CHEM' 


#%%-------------------------- PROCESSING --------------------------------------

print('Reading biomass burning emissions from ' + rootPath+'/Inputs')

if cutORnotFINN == 1:    
    bb = cutFINN(rootPath,year)
    
else:   
    file_path = [filename for filename in os.listdir(rootPath+'/Inputs') if 
                 filename.startswith("BR_FINN_"+str(year))]
    df_fire = pd.read_csv(rootPath+'/Inputs/'+file_path[0])
    bb = gpd.GeoDataFrame(
        df_fire, geometry=gpd.points_from_xy(df_fire.LONGI, df_fire.LATI))
    bb.crs = "EPSG:4326"
 
domain = Polygon(zip([loni,loni,lonf,lonf],[lati,latf,latf,lati]))
domain = gpd.GeoDataFrame(index=[0],geometry=[domain])
domain.crs = "EPSG:4326"
pip_mask = bb.within(domain.iloc[0,0])
bb = bb.loc[pip_mask]
bb=bb.reset_index(drop=True)
    
# Collecting emissions' centroids 
centerBB = bb.geometry.centroid
centerBB.to_crs("EPSG:4326")

# Condition for each speciation profile
if speciation =='GEOS-CHEM':
    bb.columns = bb.columns.str.replace(' ', '')
    bb['VOC'] =  0.04782609*bb.C3H8 + 0.03260870*bb.C2H6 + 0.03043478*bb.C2H4  +\
            0.06304348*bb.ACET + 0.0782609*bb.MEK + 0.04782609*bb.ALD2 +\
                0.03260870*bb.CH2O + 0.0848*bb.BENZ + 0.1*bb.TOLU + 0.1152*bb.XYLE
    dataEmissBB = bb[['ACET','ALD2','BENZ','CH4','CO','CO2','C2H4',
                      'C2H6','C2H2','CH2O','MEK','NH3','NO','NO2','BC',
                      'TPC','OC','SO2','TOLU','VOC','XYLE']].copy() 

# Converting to the right unit of measurement
dataEmissBB = dataEmissBB/convt 

# cd to the main folder
os.chdir(rootPath)

# Creating output directory
if os.path.isdir(outPath)==0:
    os.mkdir(outPath)
    
    
print('Setting domain borders')
x = np.linspace(loni, lonf, int((lonf-loni)/deltaX))
y = np.linspace(lati, latf, int((latf-lati)/deltaY))

#Loop over each cel in x direction
polygons=[]
for ii in range(1,x.shape[0]):
    #Loop over each cel in y direction
    for jj in range(1,y.shape[0]):
        #roadClip=[]
        lat_point_list = [y[jj-1], y[jj], y[jj], y[jj-1]]
        lon_point_list = [x[ii-1], x[ii-1], x[ii], x[ii]]
        cel = Polygon(zip(lon_point_list, lat_point_list))
        polygons.append(cel)

# Creating basegridfile
baseGrid = gpd.GeoDataFrame({'geometry':polygons})
baseGrid.to_csv(outPath+'/baseGrid_'+prefix+'.csv')
baseGrid.crs = "EPSG:4326" 
print('baseGrid_'+prefix+'.csv was created at ' + outPath )


# Calling gridding function 
grids,xv,yv,xX,yY = gridding(x,y)

# Calling populatingGrid function 
dataBB = populatingGrid(dataEmissBB,centerBB,xX,yY,xv,yv)

# Name of you output file
name = 'FINNannualEmiss_'+speciation+'_'+fileId+'_'+str(deltaX)+'x'+str(deltaY)+'_'+str(year)+'.nc'

# Calling createNETCDFtemporalFINN - ANNUAL EMISSIONS
startDate = datetime.datetime(year, month, 1, 0, 0)
endDate = datetime.datetime(year, month, 1, 1, 0)
datePfct = np.arange(np.datetime64(startDate),np.datetime64(endDate),3600000000)
dates = pd.DataFrame(datePfct)   
dates['year'] = year
dates['month'] = 1
dates['day'] = 1
dates['hour'] = 00
createNETCDFtemporalFINN(outPath,name,dataBB,xv,yv,y,x,dates,month,speciation)

#%% Calling netcdf creator

# # Creating annual basis inventory - results in grams per year
# dataTempoBB1,datePfct1,dataTempoBB2,datePfct2  = populatingGridMatBB(bb,dataEmissBB,centerBB,xX,yY,year, month)
# name = 'FINNtemporalEmiss_'+str(year)+'_'+str(month)+'_0to16.nc'
# createNETCDFtemporalFINN(folder,name,dataTempoBB1,xv,yv,lat,lon,datePfct1,month)

# name = 'FINNtemporalEmiss_'+str(year)+'_'+str(month)+'_16toEnd.nc'
# createNETCDFtemporalFINN(folder,name,dataTempoBB2,xv,yv,lat,lon,datePfct2,month)


