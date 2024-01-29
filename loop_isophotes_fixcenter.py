import numpy as np           # to define our table
import math                  # for computing the maximun ellipticity and fbar
import decimal               # for computing the maximun ellipticity and fbar
#import fitsio                # for fits image
from astropy.io import fits  # isteade of fitsio
import bz2
import numpy.ma as ma
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from pylab import *
from photutils.isophote import EllipseGeometry, Ellipse         # main package for isophotes and ellipticity
from photutils import EllipticalAperture
import sys

scale=0.396

pi=22/7                      # used in fbar equation     #fbar= 2/pi [arctan(1-emax)^(-1/2) - arctan(1-emax)^(1/2)]




# First define your table for geometric parameters

d= "combined_pairs_output.dat"
data = np.loadtxt(d,  dtype={'names':('Index','ID','RA','Dec','x0', 'y0', 'a','b', 'eps', 'pa','mag','magerr'),'formats': ('i','O','O','O','f','f', 'f', 'f', 'f', 'f', 'f','f')}, skiprows=1)

# then define what each column has

Index = data['Index']
ID=data['ID']
RA= data['RA']
Dec=data['Dec']
x0 = data['x0']
y0 = data['y0']
sma= data['a']
b = data['b']
eps = data['eps']
pa = data['pa']   # we have to convert the PA into radians
mag=data['mag']
mag_error=data['magerr']



for row in data:

    [Index,ID, RA, Dec, x0, y0, sma, b, eps, pa, mag, mag_error] = row
    if Index > 226:
        pa = pa * np.pi/180
        print ("Iter %s..." % Index)  #print the object name before its output isophote table
        g = EllipseGeometry(x0, y0, sma, eps, pa, fix_center=True)
        hdu=fits.open(ID)
        dat = hdu[0].data

        ell = Ellipse(dat, geometry=g)

        isolist = ell.fit_image()
        #print(isolist.to_table())
        print("E_all_len", len(isolist.eps))
        E_all=isolist.eps
        sma_all= isolist.sma
        print("sma_all_len", len(sma_all))
        sma_max=max(sma_all)
        print('max_sma=', sma_max)
        pa_li=isolist.pa
        print("pa_len:",len(pa_li))


        E1_li= []       # corresponding values of E in the first third of sma
        sma1_li=[]      #for the first third of the sma
        smabar1_li= []  # sma region after the first Emax
        Ebar1_li= []    # E region after Enmax1 to identify the drop (Emin)
        sma2_li=[]      #for the second third of the sma
        E2_li= []       # corresponding values of E in the second third of sma
        smabar2_li= []  # sma region after the second Emax
        Ebar2_li= []    # E region after Emax2
        #pa1_li=[]
        #pa2_li=[]
        #pa3_li=[]
        #pa4_li=[]


        for i in range(len(sma_all)):
            if 3 < sma_all[i] <=(max(sma_all)*0.33):  ## 0.33 is 1/3 the length of the sma
                sma_1=sma_all[i]
                #print('sma1', sma_1)
                sma1_li.append(sma_1)
                #print('sma1_li', sma1_li)
                E1= E_all[i]
                #print('E1:', E1)
                E1_li.append(E1)
                #print('E1_li', E1_li)
                E1_array=np.array(E1_li)
                #pa1= pa_all[i]
                #pa1_li.append(pa1)
                #pa1_array=np.array(pa1_li)



            if (max(sma_all)*0.33) <=  sma_all[i] <= (max(sma_all)*0.67):   ## 0.67 is 2/3 the length of the sma
                sma_2=sma_all[i]
                #print('sma2', sma_2)
                sma2_li.append(sma_2)
                #print('sma2_li', sma2_li)
                E2= E_all[i]
                #print('E2:', E2)
                E2_li.append(E2)
                #print('E2_li', E2_li)
                E2_array=np.array(E2_li)
                #pa2= pa_all[i]
                #pa2_li.append(pa2)
                #pa2_array=np.array(pa2_li)



                E1=Emax1= np.max(E1_array)
                print('Emax1:', Emax1)
                sma1_x = sma1_li[np.argmax(E1_array)]
                #print('sma1_x',sma1_x)
                sma_emax1=np.interp(E1,E1_li,sma1_li)        # to estimate the value of sma corresponding to Emax
                print('sma_emax1:',sma_emax1)
                #pa_max1=np.interp(sma_emax1,sma1_li,pa1_li)
                #print('pa_max1:', pa_max1)
                #print('E1_len', len(E1_li))
                E2=Emax2= np.max(E2_array)
                print('Emax2:', Emax2)
                sma2_x = sma2_li[np.argmax(E2_array)]  # to estimate the number of array (e.g [13]) that has sma corresponding to Emax
                sma_emax2=np.interp(E2,E2_li,sma2_li)        # to estimate the value of sma corresponding to Emax
                print('sma_emax2:',sma_emax2)
                #print('sma2_x',sma2_x)
                #pa_max2=np.interp(sma_emax2,sma2_li,pa2_li)
                #P1=pa_max1+(10*np.pi/180)
                #P2=pa_max2+(10*np.pi/180)
                #sma_P1=np.interp(P1,pa1_li,sma1_li)   #####
                #sma_P2=np.interp(P2,pa2_li,sma2_li)   #####


        for i in range(len(sma_all)):
            if sma1_x <= sma_all[i]  <= (sma_max*0.33):
                sma1_bar=sma_all[i]
                smabar1_li.append(sma1_bar)
                #print('sma1_bar',smabar1_li)
                Ebar1=E_all[i]
                Ebar1_li.append(Ebar1)
                E3_array=np.array(Ebar1_li)
                #print('Ebar1_li', Ebar1_li)
                #print(type (Ebar1_li))
                #Emin1=min(Ebar1_li)
                #pa3= pa_all[i]
                #pa3_li.append(pa3)
                #pa3_array=np.array(pa3_li)



            if sma2_x <= sma_all[i] <= (sma_max*0.67):
                sma2_bar=sma_all[i]
                smabar2_li.append(sma2_bar)
                #print('sma2_bar',smabar2_li)
                Ebar2=E_all[i]
                Ebar2_li.append(Ebar2)
                #print('Ebar2_li', Ebar2_li)
                #pa4= pa_all[i]
                #pa4_li.append(pa4)
                #pa4_array=np.array(pa4_li)

                E3=Emin1=np.min(E3_array)
                print('Emin1',E3)
                sma_emin1=np.interp(E3,Ebar1_li, smabar1_li)
                print('sma_emin1:',sma_emin1)
                #pa_min1=np.interp(sma_emin1,smabar1_li,pa3_li)

                E4=Emin2=np.min(Ebar2_li)
                sma_emin2=np.interp(E4,Ebar2_li,smabar2_li)
                print('sma_emin2:',sma_emin2)
                #pa_min2=np.interp(sma_emin2,smabar2_li,pa4_li)


        Edelta1=Emax1-Emin1
        Edelta2=Emax2-Emin2
        print('ED1',Edelta1)
        print('ED2', Edelta2)
        smadelta1=sma_emax1-sma_emin1
        smadelta2=sma_emax2-sma_emin2
        print('SD1', smadelta1)
        print('SD2', smadelta2)
        lbar1= smadelta1*scale
        print("Bar length_1", lbar1)
        lbar2= smadelta2*scale
        print("Bar length_2", lbar2)
        Sbar1= (2/pi)*(np.arctan(1-(Emax1))**(-1/2) - np.arctan(1-(Emax1))**(1/2))
        Sbar2= (2/pi)*(np.arctan(1-(Emax2))**(-1/2) - np.arctan(1-(Emax2))**(1/2))




        ####Equation of bar length###
        ####Lbar,0 =AVG(Lpeak,MIN(Lmin,LPA))###

        if Edelta1 >= 0.1:
            flag1=round(Sbar1,2)
            print('Sbar1:', Sbar1)

        elif Edelta1 < 0.1:
             flag1='NF'

        if Edelta2 >= 0.1:
            flag2=round(Sbar2,2)
            print('Sbar2:', Sbar2)

        elif Edelta2 < 0.1:
            flag2='NF'

        if flag1 == flag2 == 'NF':
            Bar='F'
        else:
            Bar='T'

        f_out = open('iso_%s.dat'%Index, 'w')
        print ("#Index ID RA Dec  x0  y0  sma  eps  pa edelta1 edelta2  Bar  ", file=f_out)
        print(Index,ID, RA, Dec, x0,  y0,  sma,  eps, round(pa/np.pi*180.,1), round(Edelta1,2),round(Edelta2,2), Bar,   file=f_out)



        plt.figure(2)
        fig, axs = plt.subplots(2)
        axs[0].errorbar(isolist.sma, isolist.eps, yerr=isolist.ellip_err, fmt='o', markersize=4)
        axs[0].set( xlabel='Semimajor axis ($pixels$)', ylabel='Ellipticity')
        axs[1].errorbar(isolist.sma, isolist.pa/np.pi*180., yerr=isolist.pa_err/np.pi* 80., fmt='o', markersize=4)
        axs[1].set(xlabel='Semimajor axis ($pixels$)', ylabel='PA (deg)')
        savefig('eps_%s.png'%Index)
        plt.close()
