# -*- coding: utf-8 -*-
"""
This script cuts the FINN biomass burningh emission inventory by states and 
regular areas. 

Inputs: 
    years: respective years of emission inventories for cutting.
    limits: coordinates in degrees to cut the emission data
    staticFolder: Foldes with FINN data and shapefile with Brazilian states

Last update: 27/08/2020

Author: Leonardo Hoinaski - leonardo.hoinaski@ufsc.br

"""
import geopandas as gpd
import pandas as pd
import os

#================================ INpUTS ======================================

def cutFINN(rootPath,year):
    
    #Open shapefile with Brazilian states
    brSHP = gpd.read_file(rootPath+'/Inputs/Brasil.shp')
    brSHP.crs = "EPSG:4326"
       
    #Reading FINN data
    file_path = [filename for filename in os.listdir(rootPath+'/Inputs') if 
                 filename.startswith("FINNv2.4_MOD_GEOSCHEM_"+str(year))]
    
    df_fire = pd.read_csv(rootPath+"/Inputs/"+file_path[0])
    gdf_fire = gpd.GeoDataFrame(
        df_fire, geometry=gpd.points_from_xy(df_fire.LONGI, df_fire.LATI))
    gdf_fire.crs = "EPSG:4326"
    
    brSHP['BR'] = 1
    brALL = brSHP.dissolve(by='BR')
    brALL.reset_index(drop=True, inplace=True)
    # Simplifying the shapefile
    brsimple = brALL.simplify(1, preserve_topology=True)
    brBuffer = brsimple.buffer(0.5, resolution=16)
    brBuffer.crs = "EPSG:4326"
    brBuffer.reset_index(drop=True, inplace=True)
    
    # Fire events inside brBuffer
    pip_mask = gdf_fire.within(brBuffer.loc[0])
    brFire = gdf_fire.loc[pip_mask]
    brFire.to_csv(rootPath+"/Inputs/BR_FINN_"+str(year)+".txt")
    #brFire.plot()
    return brFire
    
    
  