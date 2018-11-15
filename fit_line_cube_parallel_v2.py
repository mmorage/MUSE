#!/User/bin/env python
#importing system commands


import sys,os,string

#scientific packages
import astropy.io.fits as pyfits
import scipy
from scipy.optimize import curve_fit, leastsq
import numpy as np
from numpy import random, exp,sqrt
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting, polynomial
import multiprocessing
from multiprocessing import Pool,cpu_count,Queue,Manager

#time
import time
from astropy import units as u


import warnings
warnings.filterwarnings("ignore")


start_time1 = time.time()
start_time2 = time.time()


sys.tracebacklimit = 0


imagen_in=sys.argv[1]               #MUSE datacube ONLY science not variance 
l_min_izq=float(sys.argv[2])        #  minimum lambda selected to estimate the continuum at the left side of the frature
l_max_izq=float(sys.argv[3])        #  maximum lambda selecte  to estimate the continuum at the left side of the frature

l_min_der=float(sys.argv[4])        # minimum lambda selected to estimate the continuum at the right side of the frature
l_max_der=float(sys.argv[5])        # maximum lambda selected to estimate the continuum at the right side of the frature

guess_line=float(sys.argv[6])       # guess lamnda of the absorption or emission line in angstroms
guess_FWHM=float(sys.argv[7])       # estimated FWHM of the feature

orden_pol=float(sys.argv[8])        # Polinomia used to fit the continuos (recomended 1)

FILE_OUT=sys.argv[9]                # Output file

REST_FRAME=float(sys.argv[10])      # Rest frame line or feature

# output file: 
# x,y,Maximum_value,lambda_measured,FWHM,gaussian_integration,velocity





def ajusta(x_cube,y_cube,imagen_in,l_min_izq,l_max_izq,l_min_der,l_max_der,guess_line,guess_FWHM,orden_pol):
    sys.tracebacklimit = 0
    #scientific packages
    import astropy.io.fits as pyfits
    import scipy
    from scipy.optimize import curve_fit, leastsq
    import numpy as np
    from numpy import random, exp,sqrt
    import matplotlib.pyplot as plt
    from astropy.modeling import models, fitting, polynomial

    def Lee_cubo(spectra,XX,YY):
        global imagen
        imagen=pyfits.getdata(spectra,header=False)
        header=pyfits.getheader(spectra)
        
        #print len(imagen)
        #empty array
        Lambda_t=[]
        Flux_t=[]
    
    
        for i in range (len(imagen)):
            y=imagen[i][XX][YY]
            #        x=i*header['CDELT1']+header['CRVAL1']
            x=i*header['CD3_3']+header['CRVAL3']
            Lambda_t.append(float(imagen[i][XX][YY]))
            #Flux_t.append(float(i*header['CDELT1']+header['CRVAL1']))
            Flux_t.append(float(i*header['CD3_3']+header['CRVAL3']))
            #print x,y

        Flux=np.array(Lambda_t)
        Lambda=np.array(Flux_t)
        x=Lambda
        y=Flux
        return x,y


    ##########################
    ##
    ## Funcion Region 
    ## Toma una region de un espectro entre
    ## un minimo lambda y un maximo lambda
    ## x e y corresponden a lamba y cuentas o flujo
    ##
    ##
    ## region in between lambda_int and lambda_end
    ##
    ##
    ##
    ##
    ###########################
    #
    def region(minimo,maximo,x,y):
        xar=[]
        yar=[]
        for i in range(len(x)):
            if (x[i] > minimo) and (x[i] <maximo):
                xar.append(float(x[i]))
                yar.append(float(y[i]))
            
        xar=np.array(xar)
        yar=np.array(yar)
        return xar,yar


    #########################
    #
    # Funcion Region_discontinuo
    # Toma dos regiones de un espectro separadoas entre
    # un minimo lambda y un maximo lambda a la izquierda y un 
    # un minimo lambda y un maximo lambda a la deracha de la emission o obsorption 
    # x e y corresponden a lamba y cuentas (o flujo)
    #
    #
    # Takes 2 regions: one on at the left of the feature, the other one at the right
    #
    #
    #
    #
    #
    #
    ##########################

    def region_discontinua(minimo1,maximo1,minimo2,maximo2,x,y):
        try:
            xar=[]
            yar=[]
            for i in range(len(x)):
                if ((x[i] > minimo1) and (x[i] <maximo1)) or ((x[i] > minimo2) and (x[i] <maximo2)):
                    xar.append(float(x[i]))
                    yar.append(float(y[i]))
            
            xar=np.array(xar)
            yar=np.array(yar)
            return xar,yar
        except:
            pass

        
    #######
    # poly_fit, fitea un polinomio a datos
    # xp e yp correspondend a x e y a ser fiteado
    #
    # en este caso corresponden al x e y del output de la
    # region discontinua
    #
    #
    #
    #  Polinomial fit to the continuum
    #
    #######

    def poly_fit(xp,yp,grado_pol):
        try:
            t_init = polynomial.Polynomial1D(degree=int(grado_pol))
            fit_t = fitting.LevMarLSQFitter()
            t = fit_t(t_init, xp, yp)
            return t
        except:
            pass

        
    

    #calclulo original
    #First reading    

    
    x_sci,y_sci=Lee_cubo(imagen_in,x_cube,y_cube)
    x=x_sci 
    y=y_sci 
    xspec_o,yspec_o=region(l_min_izq,l_max_der,x,y)

    #################



    ###Fitting regions with a polynomio


    x_cont_o,y_cont_o=region_discontinua(l_min_izq,l_max_izq,l_min_der,l_max_der,x,y)


    cont3_o=poly_fit(x_cont_o,y_cont_o,orden_pol)
    
    res3_o= -cont3_o(x_cont_o)+ y_cont_o
    
    #se aplica el polinomio al espectro en la zona de interes
    # polinom to the continuum
    res4_o= -cont3_o(xspec_o)+yspec_o
    
    #####################
    #
    # Normalization!!!
    #
    ######################
    res4_oN= yspec_o/cont3_o(xspec_o)
    
    
    ######################
    
    
    ################
    
    #iteracion 1
    t_init4_o = models.Gaussian1D(amplitude=1, mean=guess_line, stddev=guess_FWHM)
    fit_t4_o = fitting.LevMarLSQFitter()
    t4_o = fit_t4_o(t_init4_o, xspec_o, res4_o)
    
    a_science=t4_o.mean.value
    b_science=t4_o.stddev.value
    Amplitud=t4_o.amplitude.value
    

    import numpy
    from scipy.optimize import curve_fit
    import matplotlib.pyplot as plt
    import math
        
    CWt= t4_o.mean.value
    FWHMt=2*sqrt(2*math.log(2))*t4_o.stddev.value
    At=t4_o.amplitude.value


    xspec_o,yspec_o=region(CWt-5*FWHMt,CWt+5*FWHMt,x,y)    
    x_cont_o,y_cont_o=region_discontinua(CWt-5*FWHMt,CWt-3*FWHMt,CWt+3*FWHMt,CWt+5*FWHMt,x,y)
    cont3_o=poly_fit(x_cont_o,y_cont_o,orden_pol)
    res4_o= -cont3_o(xspec_o)+yspec_o


    #iteracion 2
    t_init4_o = models.Gaussian1D(amplitude=Amplitud, mean=a_science, stddev=b_science)
    fit_t4_o = fitting.LevMarLSQFitter()
    t4_o = fit_t4_o(t_init4_o, xspec_o, res4_o)

    a_science=t4_o.mean.value
    b_science=t4_o.stddev.value
    Amplitud=t4_o.amplitude.value
    
    
    #iteracion 2 para gaussiana
    t_init4_o = models.Gaussian1D(amplitude=Amplitud, mean=a_science, stddev=b_science)
    fit_t4_o = fitting.LevMarLSQFitter()
    t4_o = fit_t4_o(t_init4_o, xspec_o, res4_o)
    
    residuo_o=-t4_o(xspec_o)+res4_o
    
    #print "resultados",t_init4,t4
    #print "resultados",t4_o
    a_science=t4_o.mean.value
    b_science=t4_o.stddev.value
    c_science_amplitude=t4_o.amplitude.value
    
    #Aplicando FWHM guess
    guess_FWHM_gauss=-float(-b_science)
    #print guess_FWHM_gauss
    #exit(0)
    
    
    ##print t4.mean
    Lambda_gauss_fit_sci="{:10.3f}".format(a_science)
    Sigma_gauss_fit_sci="{:10.3f}".format(b_science)
    #print  Lambda_gauss_fit_sci, Sigma_gauss_fit_sci
    #print a ,b
    
    central_wavelenght = t4_o.mean.value
    FWHM=t4_o.stddev.value
    Amplitude=t4_o.amplitude.value
    
        
    
    import numpy
    from scipy.optimize import curve_fit
    import matplotlib.pyplot as plt
    import math
    
    # Define model function to be used to fit to the data above:
    def gauss(x, *p):
        A, mu, sigma = p
        return A*numpy.exp(-(x-mu)**2/(2.*sigma**2))
    
        
    CW= t4_o.mean.value
    FWHM=2*sqrt(2*math.log(2))*t4_o.stddev.value
    A=t4_o.amplitude.value
    
    
    #INTEGRAL
    from scipy import integrate
    args = A,CW,FWHM
    #print args
    # 
    ## Integrate myfunc() from 0.5 to 1.5
    #print  CW,CW-4*FWHM, CW+4*FWHM,FWHM
    #results = integrate.quad(gauss, min(x_sci), max(x_sci), args)
    results = integrate.romberg(gauss, CW-5*FWHM, CW+5*FWHM, args)
    
    
    
 
    #print x_cube,y_cube,A,CW,FWHM,results
    #print  Lambda_gauss_fit_sci, Sigma_gauss_fit_sci

    

    ############
    #
    #
    # Uncomment the lines for plot fitted visualization
    #
    #
    ############
    
    #plt.plot(xspec_o,res4_o, 'c-',lw=1,label='gaus')
    #plt.plot(xspec_o,t4_o(xspec_o),'r-', lw=2,label='gauss')
    #plt.pause(0.005)
    #plt.clf()
    #
    # until here
    #########################################

    RESULTADO=int(x_cube),int(y_cube),float(A),float(CW),float(FWHM),float(results)

    linea=REST_FRAME * u.AA 
    linea_rest_frame=linea
    doppler=u.doppler_optical(linea_rest_frame)

    vel=(float(CW)*u.AA).to(u.km/u.s, equivalencies=doppler).value
    
    f=open(FILE_OUT,'a')    
    
    f.write(str(RESULTADO[0])+' '+str(RESULTADO[1])+' '+str(RESULTADO[2])+' '+str(RESULTADO[3])+' '+str(RESULTADO[4])+' '+str(RESULTADO[5])+' ' +str(vel) +'\n')
    f.close()
    #return


# x,y,Maximum_value,lambda_measured,FWHM,gaussian_integration,velocity


########################
#
#
# Cambio de pp por multiprocessing
# parallel work
#
#########################




def main():
    sys.tracebacklimit = 0
    subprocs = []
    for i in range (0,315): #315  #lenght y of the data cube
        for j in range(0,314): #314 lenght x of the data cube
            archivo=file
            #print archivo
            try:
                p=multiprocessing.Process(target=ajusta,args=(i,j,imagen_in,l_min_izq,l_max_izq,l_min_der,l_max_der,guess_line,guess_FWHM,orden_pol))
                
            except: 
                pass
            subprocs.append(p)
    print("Starting all processes now")

    # now start them all
    for p in subprocs: p.start()
    
    # wait on them all
    for p in subprocs: p.join()#

    return True

    
main()


#f.close()
print "Total execution time: %.2f min" % ((time.time() - start_time1)/60)

