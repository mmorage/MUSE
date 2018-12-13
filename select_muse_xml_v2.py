#!/usr/bin/env python


import xml.etree.ElementTree as ET
import os,sys,string
from io import StringIO
import numpy as np
from lxml.etree import parse as lparse
import os
import numpy as np
import astropy
#from astropy import io.fits as fits
from astropy.io import fits as pyfits
################
#
# usando dfits y fitsorot
#
#
##################


astrometry_wcs='astrometry_wcs_wfm_2015-09-10.fits'


FINAL=open('copia_todo.sh','w')

def find_OBJ_MUSE_xml():
    #lee xml from ESO in the current directory.
    #return File Name and OBS number
    OBSIDt=[]
    NAMEt=[]
    TYPE=[]
    #out=open(file_out,'w')
    #reading the .xml files that are located in the current directory
    files_xml= [f for f in os.listdir('.') if os.path.isfile(f)]
    for f in files_xml:
        if f.endswith(('.xml')):
            #print "reading "+f+" \n"
            ltree= lparse(f) 
            tree = ET.parse(f) #setup ETread
            root = tree.getroot() #define root
            arreglo=[]
            for associatedFiles  in root.findall('associatiatedFiles'):
                category=associatedFiles.find('association').text
                arr=str(category).split(' ')
            for node in tree.iter('file'):
                name=node.attrib.get('category')
                archivo=node.attrib.get('name')
                #selecting only the xml related to OBJECTS
                if name=="OBJECT": 
                    #reading the header and selecting OBS.ID
                    header=pyfits.getheader(str(f)[:-4]+'.fits.fz')
                    #print header['HIERARCH ESO OBS ID']
                    NAMEt.append((str(f)[:-4]+'.fits.fz'))
                    OBSIDt.append(str(header['HIERARCH ESO OBS ID']))
    NAME=np.array(NAMEt)
    OB=np.array(OBSIDt)
    return(NAME,OB)

# Executing the function 
print "Reading the .xml files on the directory "
OBJ,OB=find_OBJ_MUSE_xml()


def repetidos(arreglo):
    #look for repeated elements in an array
    import collections
    return([item for item, count in collections.Counter(arreglo).items() if count >= 1])

print "Cleaning the repeated names from the xml files "
ODB=repetidos(OB)
#print ODB[0]

print "writting in final file directories and science to be copied "

ARCHIVOxmlt=[]
Numero_OBSt=[]
TEMP=90
for j in range(len(ODB)):
    for i in range(len(OB)):
        if str(OB[i])==str(ODB[j]) and ODB[j]!=TEMP:
            TEMP=ODB[j]
            #print ' mkdir  '+str(ODB[j])
            FINAL.write(' mkdir  '+str(ODB[j])+' \n ')
            #            print str(ODB[j])
            ARCHIVOxmlt.append(str(OBJ[i])[:-8]+'.xml')
            Numero_OBSt.append(str(ODB[j]))
                  #print OBJ[i]+' '+ ODB[j] 
        if str(OB[i])==str(ODB[j]):
            #print 'cp '+OBJ[i]+' '+ODB[j]+'\.' 
            FINAL.write('cp '+OBJ[i]+' '+ODB[j]+'/. \n ') 
archivo_xml=np.array(ARCHIVOxmlt)
Numero_OBS=np.array(Numero_OBSt)

#print" aCACAACCAC\n\n\n\n\n" 
#print archivo_xml







def lee_ESO_MUSE_xml(test,out):
    ltree= lparse(test)
    tree = ET.parse(test) #setup ETread
    root = tree.getroot() #define root 
    #print root.tag
    #print root.attrib
    #for child in root:  
    #     print child.tag
    #    print child.attrib


    f=open(out,'w')
    
    ##
    #
    #
    arreglo=[]
    arreglo_HD=[]
    for associatedFiles  in root.findall('associatiatedFiles'):
        category=associatedFiles.find('association').text
        arr=str(category).split(' ')
#        print arr[1]
    for node in tree.iter('file'):
        name=node.attrib.get('category')
        archivo=node.attrib.get('name')
        if name=="OBJECT": 
            arreglo.append(str(archivo))
            arreglo_HD.append(str(name))

            #print   archivo,name
            f.write(str(archivo)+' '+str(name)+' \n' )
        if name=="BIAS": 
            arreglo.append(str(archivo))
            arreglo_HD.append(str(name))

            #print   archivo,name
            f.write(str(archivo)+' '+str(name)+' \n' )
        if name=="FLAT": 
            arreglo.append(str(archivo))
            arreglo_HD.append(str(name))
            #print    archivo,name
            f.write(str(archivo)+' '+str(name)+' \n' )
        if name=="BADPIX_TABLE":
            arreglo.append(str(archivo)) 
            arreglo_HD.append(str(name))
            #print     archivo,name
            f.write(str(archivo)+' '+str(name)+' \n' )
        if name=="LINE_CATALOG":
            arreglo.append(str(archivo))
            arreglo_HD.append(str(name))
            #print    archivo,name
            f.write(str(archivo)+' '+str(name)+' \n' )  
#        if name=="SKY_FLAT":
	if name=="SKYFLAT":
            arreglo.append(str(archivo))    
            arreglo_HD.append(str(name))
            #print    archivo,name
            f.write(str(archivo)+' '+str(name)+' \n' )
        if name=="LAMP_FLAT":
            arreglo.append(str(archivo)) 
            arreglo_HD.append(str(name))
            #print    archivo,name
            f.write(str(archivo)+' '+str(name)+' \n' ) 
#        if name=="ARC ILLUM":
#            arreglo.append(str(archivo)) 
#            print   archivo,name 
        if name=="ARC":
            arreglo_HD.append(str(name))
            arreglo.append(str(archivo)) 
            #print   archivo,name
            f.write(str(archivo)+' '+str(name)+' \n' ) 
        if name=="SCI_SINGLE":
            arreglo_HD.append(str(name))
            arreglo.append(str(archivo))
            #print   archivo,name 
            f.write(str(archivo)+' '+str(name)+' \n' )
        if name=="VIGNETTING_MASK":
            arreglo_HD.append(str(name))
            arreglo.append(str(archivo))
            #print    archivo,name
            f.write(str(archivo)+' '+str(name)+' \n' )  
        if name=="GEOMETRY_TABLE":
            arreglo_HD.append(str(name))
            arreglo.append(str(archivo))
            #print    archivo,name
            f.write(str(archivo)+' '+str(name)+' \n' )  
        if name=="EXTINCT_TABLE":
            arreglo_HD.append(str(name))
            arreglo.append(str(archivo))
            #print   archivo,name 
            f.write(str(archivo)+' '+str(name)+' \n' ) 
        if name=="FILTER_LIST":
            arreglo_HD.append(str(name))
            arreglo.append(str(archivo))
            #print   archivo,name 
            f.write(str(archivo)+' '+str(name)+' \n' ) 
        if name=="LSF_PROFILE":
            arreglo_HD.append(str(name))
            arreglo.append(str(archivo))
            #print   archivo,name 
            f.write(str(archivo)+' '+str(name)+' \n' ) 
        if name=="SKY_LINES":
            arreglo_HD.append(str(name))
            arreglo.append(str(archivo)) 
            #print   archivo,name 
            f.write(str(archivo)+' '+str(name)+' \n' )
        if name=="STD_FLUX_TABLE":
            arreglo_HD.append(str(name))
            arreglo.append(str(archivo)) 
            #print   archivo,name
            f.write(str(archivo)+' '+str(name)+' \n' )
        if name=="ASTROMETRY_REFERENCE":
            arreglo_HD.append(str(name))
            arreglo.append(str(archivo))
            #print   archivo,name
            f.write(str(archivo)+' '+str(name)+' \n' )
        if name=="ASTROMETRY_WCS":
            arreglo_HD.append(str(name))
            arreglo.append(str(archivo)) 
            #print   archivo,name
            f.write(str(archivo)+' '+str(name)+' \n' )
        if name=="BADPIX_TABLE":
            arreglo_HD.append(str(name))
            arreglo.append(str(archivo)) 
            #print   archivo,name
            f.write(str(archivo)+' '+str(name)+' \n' )
        if name=="ILLUM": 
            arreglo_HD.append(str(name))
            arreglo.append(str(archivo))
            #print   archivo,name
            f.write(str(archivo)+' '+str(name)+' \n' )
        if name=="STD": 
            arreglo_HD.append(str(name))
            arreglo.append(str(archivo))
            #print   archivo,name
            f.write(str(archivo)+' '+str(name)+' \n' )
            
        if name=="FILTER_LIST": 
            arreglo_HD.append(str(name))
            arreglo.append(str(archivo))
            #print   archivo,name
            f.write(str(archivo)+' '+str(name)+' \n' )

        if name=="DARK": 
            arreglo_HD.append(str(name))
            arreglo.append(str(archivo))
            #print archivo,name
            f.write(str(archivo)+' '+str(name)+' \n' )

    arr_archivo=np.array(arreglo)
    arr_nombre=np.array(arreglo_HD)
    f.close()    
    #return(arreglo,arreglo_HD)
    return(arr_archivo,arr_nombre)


print "Reading inside the .xml files on the directory" 
for i in range(len(archivo_xml)):
    a,b=lee_ESO_MUSE_xml(archivo_xml[i],str(archivo_xml[i])[:-4]+'.dat')   
    c=repetidos(a)
    for j in range(len(c)):
        print c[j]+' '+Numero_OBS[i]+' '+archivo_xml[i] #+ b[i]
        #for k in range(len(OBJ)):
        #    
        #print c[j]+'.fits.fz '+str(archivo_xml[i])[:-4]+'/.'+archivo_xml[i]+' '+Numero_OBS[i]
        #print 'cp '+c[j]+'.fits.fz '+Numero_OBS[i]+'/. ' 
        
        if str(c[j]).startswith("M.M"):
            FINAL.write('cp -u '+c[j]+'.fits '+Numero_OBS[i]+'/.  \n' )
        else:    
            FINAL.write('cp -u '+c[j]+'.fits.fz '+Numero_OBS[i]+'/.  \n' )

    #copy of the files from the pipelines
    FINAL.write('cp -u line_catalog.fits  '+Numero_OBS[i]+'/.  \n' ) 
    FINAL.write('cp -u sky_lines.fits  '+Numero_OBS[i]+'/.  \n' )  
    FINAL.write('cp -u extinct_table.fits  '+Numero_OBS[i]+'/.  \n' )  
    FINAL.write('cp -u filter_list.fits  '+Numero_OBS[i]+'/.  \n' )  
    FINAL.write('cp -u ordena_v2.py  '+Numero_OBS[i]+'/.  \n' )  
    FINAL.write('cp '+astrometry_wcs+' '+Numero_OBS[i]+'/.  \n' )  


FINAL.close()

print "files created copia_todo.sh and several files.dat" 

#os.system('cat copia_todo.sh ')

#os.system('csh copia_todo.sh ')

print "\n\nMarcelo Mora Nov-30-2018 v2.0"
