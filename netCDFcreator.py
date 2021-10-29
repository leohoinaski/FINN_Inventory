#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------
                            netCDFcreator.py
                            
This function creates the netCDF files ready to use in CMAQ from FINN inventory.

Inputs:
    
    folder: folter to output files
    
    name: output names
    
    data: matrix with data ready to convert in netCDF
    
    xv, yv: meshgrid outputs - grid definition
    
    lat, lon = grid latitude and longitudes
    
    year: respective years of emission inventories
    
    month: respective month of emission inventories
    
    dates: dates from emission file
    
    specation = speciation profile used. It could be GEOS-CHEM, 
            MOZART or SAPRC99
    
    outType = "Annual" for annual basis and "Temporal" for g/s 

Outputs:
        
    netdCDF files

        

Last update = 21/10/2021

Author: Leonardo Hoinaski - leonardo.hoinaski@ufsc.br
---------------------------------------------------------------
"""
import netCDF4 as nc4
import numpy as np
import datetime


def createNETCDFtemporalFINN(folder,name,data,xv,yv,lat,lon,dates,month,
                             speciation,outType):
    cdate = datetime.datetime.now()
    cdateStr = int(str(cdate.year)+str(cdate.timetuple().tm_yday))
    ctime = int(str(cdate.hour)+str(cdate.minute)+str(cdate.second))
    tflag = np.empty([dates.shape[0],62,2],dtype='i4')
    
    for ii in range(0,dates.shape[0]):
        tflag[ii,:,0]=int(dates['year'][0]*1000 + dates.iloc[ii,0].timetuple().tm_yday)
        tflag[ii,:,1]=int(str(dates.iloc[ii,0].hour)+'0000')
    
    sdate =  dates['year'][0]*1000 + dates.iloc[0,0].timetuple().tm_yday  
    
    f2 = nc4.Dataset(folder+'/'+name,'w', format='NETCDF3_CLASSIC') #'w' stands for write    
    #Add global attributes
    f2.IOAPI_VERSION ='$Id: @(#) ioapi library version 3.1 $'
    f2.EXEC_ID = '???????????????'
    f2.FTYPE =  1
    f2.CDATE= cdateStr
    f2.CTIME= ctime
    f2.WDATE= cdateStr
    f2.WTIME= ctime
    f2.SDATE= sdate
    f2.STIME= 0
    f2.TSTEP= 10000
    f2.NTHIK= 1
    f2.NCOLS= np.size(xv,1)-1
    f2.NROWS= np.size(yv,0)-1
    f2.NLAYS= 1
    f2.NVARS= 62 #dataEmiss.shape[1]
    f2.GDTYP= 1
    f2.P_ALP= -10
    f2.P_BET= 0
    f2.P_GAM= np.mean(xv)
    f2.XCENT= np.mean(xv)
    f2.YCENT= np.mean(yv)
    f2.XORIG= xv.min()
    f2.YORIG= yv.min()
    f2.XCELL= xv[0,1] - xv[0,0]
    f2.YCELL= yv[1,0] - yv[0,0]
    f2.VGTYP= -1
    f2.VGTOP= 0.0
    f2.VGLVLS= [0,0]
    
    f2.GDNAM= 'SE53BENCH'       
    f2.UPNAM= 'M3WNDW'   
    strVAR = ' ACET            ACROLEIN        ALD2           \n\
    ALD2_PRIMARY    ALDX            BENZ            BUTADIENE13\n\
    CH4             CH4_INV         CL2             CO\n\
    CO2_INV         ETH             ETHA            ETHY\n\
    ETOH            FORM            FORM_PRIMARY    HCL\n\
    HONO            IOLE            ISOP            KET\n\
    MEOH            N2O_INV         NAPH            NH3\n\
    NH3_FERT        NO              NO2             NVOL\n\
    OLE             PAL             PAR             PCA\n\
    PCL             PEC             PFE             PH2O\n\
    PK              PMC             PMG             PMN\n\
    PMOTHR          PNA             PNCOM           PNH4\n\
    PNO3            POC             PRPA            PSI\n\
    PSO4            PTI             SO2             SOAALK\n\
    SULF            TERP            TOL             UNK\n\
    UNR             VOC_INV         XYLMN'       
    f2.VAR_LIST=strVAR
    f2.FILEDESC= 'Merged emissions output file from Mrggrid'
    f2.HISTORY ='' 
       
    # # Specifying dimensions
    #tempgrp = f.createGroup('vehicularEmissions_data')
    f2.createDimension('TSTEP', dates.shape[0])
    f2.createDimension('DATE-TIME', 2)
    f2.createDimension('LAY', 1)
    f2.createDimension('VAR', 62)
    f2.createDimension('ROW', len(lat)-1)
    f2.createDimension('COL', len(lon)-1)
    
    # Building variables
    TFLAG = f2.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
    print(f2['TFLAG'])
    ACET = f2.createVariable('ACET', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    ACROLEIN = f2.createVariable('ACROLEIN', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    ALD2 = f2.createVariable('ALD2', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    ALD2_PRIMARY = f2.createVariable('ALD2_PRIMARY', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    ALDX = f2.createVariable('ALDX', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    BENZ = f2.createVariable('BENZ', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    BUTADIENE13 = f2.createVariable('BUTADIENE13', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    CH4 = f2.createVariable('CH4', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    CH4_INV = f2.createVariable('CH4_INV', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    CL2 = f2.createVariable('CL2', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    CO = f2.createVariable('CO', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    CO2_INV = f2.createVariable('CO2_INV', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    ETH = f2.createVariable('ETH', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    ETHA = f2.createVariable('ETHA', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    ETHY = f2.createVariable('ETHY', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    ETOH = f2.createVariable('ETOH', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    FORM = f2.createVariable('FORM', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    FORM_PRIMARY = f2.createVariable('FORM_PRIMARY', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    HCL = f2.createVariable('HCL', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    HONO = f2.createVariable('HONO', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    IOLE = f2.createVariable('IOLE', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    ISOP = f2.createVariable('ISOP', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    KET = f2.createVariable('KET', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    MEOH = f2.createVariable('MEOH', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    N2O_INV = f2.createVariable('N2O_INV', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    NAPH = f2.createVariable('NAPH', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    NH3 = f2.createVariable('NH3', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    NH3_FERT = f2.createVariable('NH3_FERT', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    NO = f2.createVariable('NO', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    NO2 = f2.createVariable('NO2', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    NVOL = f2.createVariable('NVOL', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    OLE = f2.createVariable('OLE', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PAL = f2.createVariable('PAL', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PAR = f2.createVariable('PAR', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PCA = f2.createVariable('PCA', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PCL = f2.createVariable('PCL', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PEC = f2.createVariable('PEC', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PFE = f2.createVariable('PFE', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PH2O = f2.createVariable('PH2O', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PK = f2.createVariable('PK', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PMC = f2.createVariable('PMC', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PMG = f2.createVariable('PMG', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PMN = f2.createVariable('PMN', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PMOTHR = f2.createVariable('PMOTHR', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PNA = f2.createVariable('PNA', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PNCOM = f2.createVariable('PNCOM', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PNH4 = f2.createVariable('PNH4', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PNO3 = f2.createVariable('PNO3', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    POC = f2.createVariable('POC', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PRPA = f2.createVariable('PRPA', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PSI = f2.createVariable('PSI', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PSO4 = f2.createVariable('PSO4', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    PTI = f2.createVariable('PTI', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    SO2 = f2.createVariable('SO2', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    SOAALK = f2.createVariable('SOAALK', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    SULF = f2.createVariable('SULF', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    TERP = f2.createVariable('TERP', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    TOL = f2.createVariable('TOL', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    UNK = f2.createVariable('UNK', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    UNR = f2.createVariable('UNR', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    VOC_INV = f2.createVariable('VOC_INV', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    XYLMN = f2.createVariable('XYLMN', 'f4', ('TSTEP', 'LAY', 'ROW','COL'))
    
    
    # Passing data into variables
    if speciation=='GEOS-CHEM':
        print(tflag.shape)
        TFLAG[:,:,:] = tflag
        #TFLAG[:] = dates.iloc[:,0].dt.strftime('%Y-%m-%d %H:%M:%S')
        ACET[:,:,:,:] =  data[:,0,:,:]
        ACROLEIN[:,:,:,:] =np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        ALD2[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        ALD2_PRIMARY[:,:,:,:] = data[:,1,:,:]
        ALDX[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        BENZ[:,:,:,:] = data[:,2,:,:]
        BUTADIENE13[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        CH4[:,:,:,:] = data[:,3,:,:]
        CH4_INV[:,:,:,:] = data[:,3,:,:]
        CL2[:,:,:,:] =  np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        CO[:,:,:,:] = data[:,4,:,:]
        CO2_INV[:,:,:,:] =  data[:,5,:,:]
        ETH[:,:,:,:] =  data[:,6,:,:]
        ETHA[:,:,:,:] =  data[:,7,:,:]
        ETHY[:,:,:,:] =  data[:,8,:,:]
        ETOH[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        FORM[:,:,:,:] = data[:,9,:,:]
        FORM_PRIMARY[:,:,:,:] = data[:,9,:,:]
        HCL[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        HONO[:,:,:,:] =np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        IOLE[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        ISOP[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])  
        KET[:,:,:,:] = data[:,10,:]  
        MEOH[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        N2O_INV[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        NAPH[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])    
        NH3[:,:,:,:] = data[:,11,:,:]
        NH3_FERT[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        NO[:,:,:,:] = data[:,12,:,:]
        NO2[:,:,:,:] = data[:,13,:,:]
        NVOL[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        OLE[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        PAL[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        PAR[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        PCA[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        PCL[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        PEC[:,:,:,:] = data[:,14,:,:]*1000
        PFE[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        PH2O[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        PK[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])    
        PMC[:,:,:,:] = data[:,15,:,:]    
        PMG[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        PMN[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        PMOTHR[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)]) #data[:,26,:,:]*1000 # está errado!!! 
        PNA[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        PNCOM[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])    
        PNH4[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])   
        PNO3[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        POC[:,:,:,:] = data[:,16,:,:]*1000
        PRPA[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        PSI[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        PSO4[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        PTI[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        SO2[:,:,:,:] = data[:,17,:,:]
        SOAALK[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        SULF[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        TERP[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        TOL[:,:,:,:] = data[:,18,:,:]
        UNK[:,:,:,:]= np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        UNR[:,:,:,:] = np.zeros([np.size(data,0),1,np.size(data,2),np.size(data,3)])
        VOC_INV[:,:,:,:] = data[:,19,:,:]
        XYLMN[:,:,:,:] = data[:,20,:,:]
   


    
    #Add local attributes to variable instances
    TFLAG.units = '<YYYYDDD,HHMMSS>'
    
    if outType == 'Annual':
    
        ACET.units = 'moles/year '
        ACROLEIN.units = 'moles/year '
        ALD2.units = 'moles/year '
        ALD2_PRIMARY.units = 'moles/year '
        ALDX.units = 'moles/year '
        BENZ.units = 'moles/year '
        BUTADIENE13.units = 'moles/year '
        CH4.units = 'moles/year '
        CH4_INV.units = 'g/year '
        CL2.units = 'moles/year ' 
        CO.units = 'moles/year '
        CO2_INV.units = 'g/year '
        ETH.units = 'moles/year '
        ETHA.units = 'moles/year '
        ETHY.units = 'moles/year '
        ETOH.units = 'moles/year '
        FORM.units = 'moles/year '
        FORM_PRIMARY.units = 'moles/year '
        HCL.units = 'moles/year '
        HONO.units = 'moles/year '
        IOLE.units = 'moles/year '
        ISOP.units = 'moles/year '
        KET.units = 'moles/year '
        MEOH.units = 'moles/year '
        N2O_INV.units = 'g/year '
        NAPH.units = 'moles/year '
        NH3.units = 'moles/year '
        NH3_FERT.units = 'moles/year '
        NO.units = 'moles/year '
        NO2.units = 'moles/year '
        NVOL.units = 'moles/year '
        OLE.units = 'moles/year '
        PAL.units = 'moles/year '
        PAR.units = 'moles/year '
        PCA.units = 'g/year '
        PCL.units = 'g/year '
        PEC.units = 'g/year '
        PFE.units = 'g/year '
        PH2O.units = 'g/year '
        PK.units = 'g/year '
        PMC.units = 'g/year ' 
        PMG.units = 'g/year '
        PMN.units = 'g/year ' 
        PMOTHR.units = 'g/year ' 
        PNA.units = 'g/year ' 
        PNCOM.units = 'g/year ' 
        PNH4.units = 'g/year ' 
        PNO3.units = 'g/year ' 
        POC.units = 'g/year ' 
        PRPA.units = 'moles/year '
        PSI.units = 'g/year ' 
        PSO4.units = 'g/year ' 
        PTI.units = 'g/year ' 
        SO2.units = 'moles/year ' 
        SOAALK.units = 'moles/year ' 
        SULF.units = 'moles/year ' 
        TERP.units = 'moles/year ' 
        TOL.units = 'moles/year ' 
        UNK.units = 'moles/year ' 
        UNR.units = 'moles/year ' 
        VOC_INV.units = 'g/year ' 
        XYLMN.units = 'moles/year '
        
    else:
        
        ACET.units = 'moles/s '
        ACROLEIN.units = 'moles/s '
        ALD2.units = 'moles/s '
        ALD2_PRIMARY.units = 'moles/s '
        ALDX.units = 'moles/s '
        BENZ.units = 'moles/s '
        BUTADIENE13.units = 'moles/s '
        CH4.units = 'moles/s '
        CH4_INV.units = 'g/s '
        CL2.units = 'moles/s ' 
        CO.units = 'moles/s '
        CO2_INV.units = 'g/s '
        ETH.units = 'moles/s '
        ETHA.units = 'moles/s '
        ETHY.units = 'moles/s '
        ETOH.units = 'moles/s '
        FORM.units = 'moles/s '
        FORM_PRIMARY.units = 'moles/s '
        HCL.units = 'moles/s '
        HONO.units = 'moles/s '
        IOLE.units = 'moles/s '
        ISOP.units = 'moles/s '
        KET.units = 'moles/s '
        MEOH.units = 'moles/s '
        N2O_INV.units = 'g/s '
        NAPH.units = 'moles/s '
        NH3.units = 'moles/s '
        NH3_FERT.units = 'moles/s '
        NO.units = 'moles/s '
        NO2.units = 'moles/s '
        NVOL.units = 'moles/s '
        OLE.units = 'moles/s '
        PAL.units = 'moles/s '
        PAR.units = 'moles/s '
        PCA.units = 'g/s '
        PCL.units = 'g/s '
        PEC.units = 'g/s '
        PFE.units = 'g/s '
        PH2O.units = 'g/s '
        PK.units = 'g/s '
        PMC.units = 'g/s ' 
        PMG.units = 'g/s '
        PMN.units = 'g/s ' 
        PMOTHR.units = 'g/s ' 
        PNA.units = 'g/s ' 
        PNCOM.units = 'g/s ' 
        PNH4.units = 'g/s ' 
        PNO3.units = 'g/s ' 
        POC.units = 'g/s ' 
        PRPA.units = 'moles/s '
        PSI.units = 'g/s ' 
        PSO4.units = 'g/s ' 
        PTI.units = 'g/s ' 
        SO2.units = 'moles/s ' 
        SOAALK.units = 'moles/s ' 
        SULF.units = 'moles/s ' 
        TERP.units = 'moles/s ' 
        TOL.units = 'moles/s ' 
        UNK.units = 'moles/s ' 
        UNR.units = 'moles/s ' 
        VOC_INV.units = 'g/s ' 
        XYLMN.units = 'moles/s '
    
    f2.close()

