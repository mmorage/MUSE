import xml.etree.ElementTree as ET
import os,sys,string
from io import StringIO
import numpy as np
from lxml.etree import parse as lparse
import os
import numpy as np



astrometry='astrometry_wcs_wfm_2015-09-10.fits'





os.system("dfits  M*.fz |fitsort  OBS.ID  dpr.type OBS.PROG.ID ORIGFILE OBS.NAME OBS.DID  >science_exposures.txt")

os.system ("dfits  *.fits  |fitsort  pro.catg >M_files.txt" )

MFILE=open('M_files.txt','r')
Mlines=MFILE.readlines()
MFILE.close()

Marchivot=[]
Mcatalogt=[]

for Mline in Mlines:
    if Mline.startswith('FILE'):
        a=1
    else:
        Mp=Mline.split()
        Marchivot.append(str(Mp[0]))   
        Mcatalogt.append(str(Mp[1]))    

Marchivo=np.array(Marchivot)      
Mcatalog=np.array(Mcatalogt)    

for j in range(len(Marchivo)):
    if Mcatalog[j]=='LINE_CATALOG':        	
        LINE_CATALOG_v=Mcatalog[j]
        FILE_LINE_CATALOG=Marchivo[j]
    if Mcatalog[j]=='FILTER_LIST':        	         	
        FILTER_LIST_v=Mcatalog[j]        	         	
        FILE_FILTER_LIST=Marchivo[j]
    if Mcatalog[j]=='EXTINCT_TABLE':
        EXTINCT_TABLE_v=Mcatalog[j]
        FILE_EXTINCT_TABLE=Marchivo[j]
    if Mcatalog[j]=='STD_FLUX_TABLE':      	
        STD_FLUX_TABLE_v=Mcatalog[j]      	
        FILE_STD_FLUX_TABLE=Marchivo[j]      	
    try:
        if Mcatalog[j]=='ASTROMETRY_REFERENCE':	
            ASTROMETRY_REFERENCE_v=Mcatalog[j]	
            FILE_ASTROMETRY_REFERENCE=Marchivo[j]
    except:
        pass
    if Mcatalog[j]=='BADPIX_TABLE':        	
        BADPIX_TABLE_v=Mcatalog[j]        	
        FILE_BADPIX_TABLE=Marchivo[j]      	
    if Mcatalog[j]=='VIGNETTING_MASK':     	
        VIGNETTING_MASK_v=Mcatalog[j]     	
        FILE_VIGNETTING_MASK=Marchivo[j]   	
    if Mcatalog[j]=='LSF_PROFILE':         	
        LSF_PROFILE_v=Mcatalog[j]         	
        FILE_LSF_PROFILE=Marchivo[j]       	
    if Mcatalog[j]=='ASTROMETRY_WCS':      	
        ASTROMETRY_WCS_v=Mcatalog[j]      	
        FILE_ASTROMETRY_WCS=Marchivo[j]    	
    if Mcatalog[j]=='GEOMETRY_TABLE':      	
        GEOMETRY_TABLE_v=Mcatalog[j]      	
        FILE_GEOMETRY_TABLE =Marchivo[j]    	

    if Mcatalog[j]=='SKY_LINES':      	
        SKY_LINES_v=Mcatalog[j]      	
        FILE_SKY_LINES =Marchivo[j]    	



print FILTER_LIST_v           
print EXTINCT_TABLE_v         
print STD_FLUX_TABLE_v        
print ASTROMETRY_REFERENCE_v  
print BADPIX_TABLE_v          
print VIGNETTING_MASK_v       
print LSF_PROFILE_v           
print ASTROMETRY_WCS_v        
print GEOMETRY_TABLE_v        
print SKY_LINES_v     	



print FILE_FILTER_LIST           
print FILE_EXTINCT_TABLE         
print FILE_STD_FLUX_TABLE        
print FILE_ASTROMETRY_REFERENCE  
print FILE_BADPIX_TABLE          
print FILE_VIGNETTING_MASK       
print FILE_LSF_PROFILE           
print FILE_ASTROMETRY_WCS        
print FILE_GEOMETRY_TABLE       
print FILE_SKY_LINES  	




FILE=open('science_exposures.txt','r')
lines=FILE.readlines()
FILE.close()

archivot=[]
OBS_IDt=[]
DTYPEt=[]
OBS_Nt=[]
ESO_NAMEt=[]
ODB_NAMEt=[]
ESO_ODB_IDt=[] 

for line in lines:
    if line.startswith('FILE'):
        a=1
    else:
        p=line.split()
        archivot.append(str(p[0]))   
        OBS_IDt.append(str(p[1]))   
        DTYPEt.append(str(p[2]))     
        OBS_Nt.append(str(p[3]))     
        ESO_NAMEt.append(str(p[4]))  
        ODB_NAMEt.append(str(p[5]))  
        ESO_ODB_IDt.append(str(p[6]))



archivo=np.array(archivot)      
OBS_ID=np.array(OBS_IDt)    
DTYPE=np.array(DTYPEt)     
OBS_N=np.array(OBS_Nt)     
ESO_NAME=np.array(ESO_NAMEt)  
ODB_NAME=np.array(ODB_NAMEt)  
ESO_ODB_ID=np.array(ESO_ODB_IDt)

MUSE_BIAS=open("muse_bias.sof","w")
MUSE_FLAT=open("muse_flat.sof","w")
MUSE_ARC=open("muse_wavecal.sof","w")
MUSE_LSF=open("muse_lsf.sof","w")
MUSE_TWL=open("muse_twilight.sof","w")
MUSE_SCB_OB=open("muse_scibasic_obj.sof","w")
MUSE_SCB_STD=open("muse_scibasic_std.sof","w")
MUSE_STD=open("muse_standard.sof","w")
MUSE_SKY=open("muse_create_sky.sof","w")

N_objects=0
N_std=0
for i in range(len(archivo)):
    if DTYPE[i]=='BIAS':
        print archivo[i], DTYPE[i]
        MUSE_BIAS.write(str(archivo[i])+'  '+str(DTYPE[i])+'\n')

    if DTYPE[i]=='FLAT,LAMP':
        print archivo[i], str(DTYPE[i])[:-5]
        #MUSE_FLAT.write(str(archivo[i])+'  '+str(DTYPE[i])+'\n')
        MUSE_FLAT.write(str(archivo[i])+'   FLAT\n')

    if DTYPE[i]=='WAVE':
        print archivo[i]+ '    ARC'
        MUSE_ARC.write(str(archivo[i])+' ARC\n')
        MUSE_LSF.write(str(archivo[i])+'  ARC\n')
    if DTYPE[i]=='FLAT,SKY':
        print archivo[i], str(DTYPE[i])[:-5]
        #MUSE_TWL.write(str(archivo[i])+'  '+str(DTYPE[i])+'\n')
        MUSE_TWL.write(str(archivo[i])+'   SKYFLAT \n')

    if DTYPE[i]=='OBJECT':
        print archivo[i], str(DTYPE[i])[:-5]
        #MUSE_TWL.write(str(archivo[i])+'  '+str(DTYPE[i])+'\n')
        MUSE_SCB_OB.write(str(archivo[i])+'   OBJECT \n')
        N_objects=N_objects +1
        
    if DTYPE[i]=='STD':
        print archivo[i], str(DTYPE[i])[:-5]
        #MUSE_TWL.write(str(archivo[i])+'  '+str(DTYPE[i])+'\n')
        MUSE_SCB_STD.write(str(archivo[i])+'   STD \n')
        N_std=N_std+1
    if DTYPE[i]=='FLAT,LAMP,ILUM':
        print archivo[i], str(DTYPE[i])[:-5]
        #MUSE_FLAT.write(str(archivo[i])+'  '+str(DTYPE[i])+'\n')
        MUSE_FLAT.write(str(archivo[i])+'   ILUM\n')



 
#BIAS

MUSE_BIAS.write(str(FILE_BADPIX_TABLE)+'  '+ str(BADPIX_TABLE_v)+' \n') 

# FLAT

MUSE_FLAT.write('MASTER_BIAS.fits'+'  '+'MASTER_BIAS'+'\n')
MUSE_FLAT.write(str(FILE_BADPIX_TABLE)+'  '+ str(BADPIX_TABLE_v)+' \n') 


# Arc muse_wavecal.sof


#MUSE_ARC.write(str(archivo[i])+'  '+str(DTYPE[i])+'\n')
MUSE_ARC.write('MASTER_BIAS.fits'+'  '+'MASTER_BIAS'+'\n')
MUSE_ARC.write('MASTER_FLAT.fits'+'  '+'MASTER_FLAT'+'\n')
#TRACE_TABLE created from FLATS
MUSE_ARC.write('TRACE_TABLE.fits'+'  '+'TRACE_TABLE'+'\n')
#line catalog.fits must be copied from the uper directiry
MUSE_ARC.write('line_catalog.fits'+'  '+'LINE_CATALOG'+'\n')


#MUSE_LSF.write(str(archivo[i])+'  ARC\n') 
MUSE_LSF.write('MASTER_BIAS.fits'+'  '+'MASTER_BIAS'+'\n')
MUSE_LSF.write('MASTER_FLAT.fits'+'  '+'MASTER_FLAT'+'\n') 
MUSE_LSF.write('TRACE_TABLE.fits'+'  '+'TRACE_TABLE'+'\n')
MUSE_LSF.write('WAVECAL_TABLE.fits'+'  '+'WAVECAL_TABLE'+'\n')   
MUSE_LSF.write(str(FILE_BADPIX_TABLE)+'  '+ str(BADPIX_TABLE_v)+' \n') 
MUSE_LSF.write('line_catalog.fits'+'  '+'LINE_CATALOG'+'\n')

#
#
# MUSE GEOMETRY IS IGNORED!
#
#

MUSE_TWL.write('MASTER_BIAS.fits'+'  '+'MASTER_BIAS'+'\n')
MUSE_TWL.write('MASTER_FLAT.fits'+'  '+'MASTER_FLAT'+'\n')
MUSE_TWL.write(str(FILE_BADPIX_TABLE)+'  '+ str(BADPIX_TABLE_v)+' \n') 
MUSE_TWL.write('TRACE_TABLE.fits  '+ 'TRACE_TABLE  \n') 
MUSE_TWL.write('WAVECAL_TABLE.fits  '+ 'WAVECAL_TABLE  \n') 
MUSE_TWL.write(str(FILE_GEOMETRY_TABLE)+'  '+str(GEOMETRY_TABLE_v)+' \n')             
MUSE_TWL.write(str(FILE_VIGNETTING_MASK)+'  '+str(VIGNETTING_MASK_v)+' \n')                           



MUSE_SCB_OB.write('MASTER_BIAS.fits'+'  '+'MASTER_BIAS'+'\n')
MUSE_SCB_STD.write('MASTER_BIAS.fits'+'  '+'MASTER_BIAS'+'\n')
MUSE_SCB_OB.write('MASTER_FLAT.fits'+'  '+'MASTER_FLAT'+'\n')
MUSE_SCB_STD.write('MASTER_FLAT.fits'+'  '+'MASTER_FLAT'+'\n')
MUSE_SCB_OB.write('TRACE_TABLE.fits  '+ 'TRACE_TABLE  \n') 
MUSE_SCB_STD.write('TRACE_TABLE.fits  '+ 'TRACE_TABLE  \n')
MUSE_SCB_OB.write('WAVECAL_TABLE.fits  '+ 'WAVECAL_TABLE  \n') 
MUSE_SCB_STD.write('WAVECAL_TABLE.fits  '+ 'WAVECAL_TABLE  \n')
 
MUSE_SCB_OB.write(str(FILE_GEOMETRY_TABLE)+'  '+str(GEOMETRY_TABLE_v)+' \n')             
MUSE_SCB_STD.write(str(FILE_GEOMETRY_TABLE)+'  '+str(GEOMETRY_TABLE_v)+' \n')             
#MUSE_SCB_OB.write('GEOMETRY_TABLE  '+str(GEOMETRY_TABLE_v)+' \n')             
#MUSE_SCB_STD.write('GEOMETRY_TABLE  '+str(GEOMETRY_TABLE_v)+' \n')             
MUSE_SCB_OB.write('TWILIGHT_CUBE.fits TWILIGHT_CUBE  \n')
MUSE_SCB_STD.write('TWILIGHT_CUBE.fits TWILIGHT_CUBE  \n')


MUSE_STD.write('PIXTABLE_STD_0001-01.fits  PIXTABLE_STD \n')
MUSE_STD.write('PIXTABLE_STD_0001-02.fits  PIXTABLE_STD \n')
MUSE_STD.write('PIXTABLE_STD_0001-03.fits  PIXTABLE_STD \n')
MUSE_STD.write('PIXTABLE_STD_0001-04.fits  PIXTABLE_STD \n')
MUSE_STD.write('PIXTABLE_STD_0001-05.fits  PIXTABLE_STD \n')
MUSE_STD.write('PIXTABLE_STD_0001-06.fits  PIXTABLE_STD \n')
MUSE_STD.write('PIXTABLE_STD_0001-07.fits  PIXTABLE_STD \n')
MUSE_STD.write('PIXTABLE_STD_0001-08.fits  PIXTABLE_STD \n')
MUSE_STD.write('PIXTABLE_STD_0001-09.fits  PIXTABLE_STD \n')
MUSE_STD.write('PIXTABLE_STD_0001-10.fits  PIXTABLE_STD \n')
MUSE_STD.write('PIXTABLE_STD_0001-11.fits  PIXTABLE_STD \n')
MUSE_STD.write('PIXTABLE_STD_0001-12.fits  PIXTABLE_STD \n')
MUSE_STD.write('PIXTABLE_STD_0001-13.fits  PIXTABLE_STD \n')
MUSE_STD.write('PIXTABLE_STD_0001-14.fits  PIXTABLE_STD \n')
MUSE_STD.write('PIXTABLE_STD_0001-15.fits  PIXTABLE_STD \n')
MUSE_STD.write('PIXTABLE_STD_0001-16.fits  PIXTABLE_STD \n')
MUSE_STD.write('PIXTABLE_STD_0001-17.fits  PIXTABLE_STD \n')
MUSE_STD.write('PIXTABLE_STD_0001-18.fits  PIXTABLE_STD \n')
MUSE_STD.write('PIXTABLE_STD_0001-19.fits  PIXTABLE_STD \n')
MUSE_STD.write('PIXTABLE_STD_0001-20.fits  PIXTABLE_STD \n')
MUSE_STD.write('PIXTABLE_STD_0001-21.fits  PIXTABLE_STD \n')
MUSE_STD.write('PIXTABLE_STD_0001-22.fits  PIXTABLE_STD \n')
MUSE_STD.write('PIXTABLE_STD_0001-23.fits  PIXTABLE_STD \n')
MUSE_STD.write('PIXTABLE_STD_0001-24.fits  PIXTABLE_STD \n')
#MUSE_STD.write('PIXTABLE_STD.fits PIXTABLE_STD.fits \n')
#MUSE_STD.write('PIXTABLE_STD.fits PIXTABLE_STD.fits \n')
MUSE_STD.write(str(FILE_EXTINCT_TABLE)+' '+str(EXTINCT_TABLE_v)+' \n')           
MUSE_STD.write(str(FILE_STD_FLUX_TABLE)+'  '+str(STD_FLUX_TABLE_v)+' \n')



MUSE_SKY.write('PIXTABLE_SKY.fits  PIXTABLE_SKY \n')           
MUSE_SKY.write(str(FILE_EXTINCT_TABLE)+' '+str(EXTINCT_TABLE_v)+' \n')           
MUSE_SKY.write('STD_RESPONSE.fits STD_RESPONSE  \n')
MUSE_SKY.write('SKY_LINES.fits SKY_LINES  \n')
#MUSE_AST.write('PIXTABLE_ASTROMETRY.fits  PIXTABLE_ASTROMETRY  \n')


           
###
#
#
# INDIVIDUAL CUBE CREATIION. DO NOT FORGET TO MV THE OUTPUTFILE!!!!
#
#
###


print N_objects
print type(N_objects)
for l in range(N_objects):
    #MUSE_SCIPOST=[]
    print("muse_scipost_"+str(l+1).zfill(4)+".sof",'w' ) 
    
    MUSE_SCiPOST=open("muse_scipost_"+str(l+1).zfill(4)+".sof",'w' )
     
    for k in range(24):
#        MUSE_SCiPOST_str(l+1).write('PIXTABLE_OBJECT_'+str(N_Object+1).zfill(4)+'-'+str(k+1).zfill(2)+'.fits  \n')
        MUSE_SCiPOST.write('PIXTABLE_OBJECT_'+str(l+1).zfill(4)+'-'+str(k+1).zfill(2)+'.fits  PIXTABLE_OBJECT\n')
#   MUSE_SCiPOST_str(l+1).write('LSF_PROFILE  LSF_PROFILE\n')        
    MUSE_SCiPOST.write('LSF_PROFILE.fits  LSF_PROFILE\n')        
    for m in range(N_std):
#        MUSE_SCiPOST_str(l+1).write('STD_RESPONSE'+str(m+1).zfill(4)+' STD_RESPONSE\n')        
        MUSE_SCiPOST.write('STD_RESPONSE_'+str(m+1).zfill(4)+'.fits  STD_RESPONSE\n')        
#        MUSE_SCiPOST_str(l+1).write('STD_TELLURIC'+str(m+1).zfill(4)+' STD_TELLURIC\n')        
        MUSE_SCiPOST.write('STD_TELLURIC_'+str(m+1).zfill(4)+'.fits  STD_TELLURIC\n')        
#    MUSE_SCiPOST_str(l+1).write(astrometry+'  ASTROMETRY_WCS \n')
#    MUSE_SCiPOST_str(l+1).write('sky_lines.fits SKY_LINES \n')
#    MUSE_SCiPOST_str(l+1).write('extinct_table.fits EXTINCT_TABLE \n')
#    MUSE_SCiPOST_str(l+1).write('filter_list.fits FILTER_LIST  \n')

    MUSE_SCiPOST.write(astrometry+'  ASTROMETRY_WCS \n')
#    MUSE_SCiPOST.write(str(FILE_ASTROMETRY_WCS)+' ASTROMETRY_WCS \n')

    MUSE_SCiPOST.write('sky_lines.fits SKY_LINES \n')
    MUSE_SCiPOST.write('extinct_table.fits EXTINCT_TABLE \n')
    MUSE_SCiPOST.write('filter_list.fits FILTER_LIST  \n')
    MUSE_SCiPOST.close()


#for l in range(len(N_objects)):
#    MUSE_SCiPOST_str(l+1).close()


####


MUSE_BIAS.close()
MUSE_FLAT.close()
MUSE_ARC.close()
MUSE_LSF.close()
MUSE_TWL.close()
MUSE_SCB_STD.close()
MUSE_SCB_OB.close()
MUSE_STD.close()
MUSE_SKY.close()




#MUSE.2015-01-10T00:41:36.458 STD

### Creating execution script
print("Creating execution script\n")
script=open('REDUCE_MUSE.sh','w')
combine=open('muse_exp_combine.sof','w')

script.write('OMP_NUM_THREADS=24 esorex --log-file=bias.log muse_bias --nifu=-1 --merge muse_bias.sof \n')
#script.write('OMP_NUM_THREADS=24 esorex --log-file=dark.log muse_dark --nifu=-1 --merge DARK_test.sof \n')
script.write('OMP_NUM_THREADS=24 esorex --log-file=flat.log muse_flat --nifu=-1 --merge muse_flat.sof \n')
script.write('OMP_NUM_THREADS=24 esorex --log-file=wavecal.log muse_wavecal --nifu=-1 --resample --residuals --merge muse_wavecal.sof \n')
script.write('OMP_NUM_THREADS=24 esorex --log-file=lsf.log muse_lsf --nifu=-1 --merge muse_lsf.sof \n')
script.write('OMP_NUM_THREADS=24 esorex --log-file=twilight.log muse_twilight muse_twilight.sof \n')
script.write('OMP_NUM_THREADS=24 esorex --log-file=sci_basic_OBJ.log muse_scibasic --nifu=-1 --merge muse_scibasic_obj.sof \n')
script.write('OMP_NUM_THREADS=24 esorex --log-file=sci_basic_STD.log muse_scibasic --nifu=-1 --merge muse_scibasic_std.sof \n')
script.write('OMP_NUM_THREADS=24 esorex --log-file=std.log muse_standard --filter=white muse_standard.sof \n')
for l in range (N_objects): 
    script.write(' OMP_NUM_THREADS=24 esorex --log-file=scipost_'+str(l+1).zfill(4)+'.log muse_scipost  --skymethod=subtract-model --save=individual,cube,skymodel --crtype=median --filter=white  muse_scipost_'+str(l+1).zfill(4)+'.sof\n')
    script.write(' mv DATACUBE_FINAL.fits DATACUBE_FINAL_'+str(l+1).zfill(4)+'.fits\n')
    script.write(' mv IMAGE_FOV_0001.fits IMAGE_FOV_'+str(l+1).zfill(4)+'END.fits\n')
    script.write(' mv PIXTABLE_REDUCED_0001.fits  PIXTABLE_REDUCED_'+str(l+1).zfill(4)+'_EXP'+str(l+1).zfill(2)+'.fits\n')
    combine.write('PIXTABLE_REDUCED_'+str(l+1).zfill(4)+'_EXP'+str(l+1).zfill(2)+'.fits  PIXTABLE_REDUCED\n')

combine.write('filter_list.fits  FILTER_LIST')
script.write('OMP_NUM_THREADS=24 esorex --log-file=combine.log muse_exp_combine --save=cube,combined muse_exp_combine.sof \n')  


combine.close()
script.close()

#mv DATACUBE_FINAL
#IMAGEN_FOV
#PIXTABLE_REDUCED

exit()


