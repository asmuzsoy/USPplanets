#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 15 21:06:54 2018

@author: anasofiauzsoy
"""

import numpy as np
import csv

earthrad = 6371000
sunrad = 695508000

kepid = []
koi = []
per = []
adivr = []
delta = []
derror = []

def constrain(teff, logg, mkep):
    if ((row[teff] != "") and (4100 < float(row[teff]) < 6100)):
        if (((row[logg]) != "") and (4.0  < float(row[logg]) < 4.9)):
            if ((row[mkep] != "") and (float(row[mkep]) < 16)):
                return True
            

with open('USPCandidates.csv', encoding='latin-1') as File:
    reader = csv.reader(File, delimiter=',')
    rownum = 0
    for row in reader:
        if rownum > 0:
            if constrain(3,6,2):
                if ((row[12] != "") and (4 < (float(row[12])*24) < 24)):
                    if ((row[23] != "") and (.84 < (float(row[23])) < 4)):
                        kepid.append(int(row[0]))
                        koi.append(int(round(float(row[1]))))
                        per.append(float(row[12]))
                        adivr.append(float(row[20]))
                        delta.append(float(row[14]))
                        derror .append(float(row[15]))
        rownum += 1 
        
print(len(kepid))

stars1 = []
snr = np.zeros(len(kepid))
teff= np.zeros(len(kepid))
logg = np.zeros(len(kepid))
feh = np.zeros(len(kepid))

with open('keplerstellar.csv') as File:
    reader = csv.reader(File, delimiter=',')
    rownum = 0
    for row in reader:
        if rownum > 0:
            if constrain(6,8,5):
                if "1" in (str(row[19])[:16]):
                    #if (row[1] != "") and ((int(row[1]) in stars1) == False):
                    if ((row[1]) != "") and "q1_q16" in str(row[0]):
                        stars1.append(int(row[1])) 
                        if row[21] != "" and row[22] != "":
                            time = float(row[21])*float(row[22])
                            if row[23] != "":
                                for i in range(len(kepid)):
                                    if int(row[1]) == kepid[i]:
                                        t = (per[i]*24)*(1/adivr[i])*(1/np.pi) # transit duration
                                        snr[i] = ((delta[i]/(float(row[23])))*np.sqrt((t*time)/(per[i]*6)))
                                        teff[i] = (float(row[6]))
                                        logg[i] = (float(row[8]))
                                        feh[i] = (float(row[10]))
                
        rownum += 1
print("Total number of stars searched: ",len(stars1))

with open('ajaa7b7ct1_ascii.csv') as File:
    reader = csv.reader(File, delimiter=',')
    rownum = 0
    for row in reader:
        if rownum > 2:
            #print(int(row[0][3:]))
            for i in range(len(kepid)):
                if int(row[0][3:]) == koi[i]:
                    teff[i] = (float(row[1]))
                    logg[i] = (float(row[4]))
                    feh[i] = (float(row[7]))        
        rownum += 1

#print(snr)
from astropy.table import Table

data = Table.read('kepler_dr2_1arcsec.fits', format='fits')
#print(data.info)

kepid1 = []
parallax = []
perror = []

for i in range(len(data['kepid'])):
    kepid1.append(data['kepid'][i])
    parallax.append(data['parallax'][i])
    perror.append(data['parallax_error'][i])

parallax2 = []
perror2 = []
for i in range(len(kepid)):
    for j in range(len(kepid1)):
        if kepid[i] == kepid1[j]:
            parallax2.append(parallax[j])
            perror2.append(perror[j])
            
# =============================================================================
# print(parallax2[kepid.index(3112129)])
# print(perror2[kepid.index(5972334)])
# =============================================================================


koi = [koi[i] for i in range(0,len(kepid)) if kepid[i] != 11030475]
adivr = [adivr[i] for i in range(0,len(kepid)) if kepid[i] != 11030475]
per = [per[i] for i in range(0,len(kepid)) if kepid[i] != 11030475]
teff = [teff[i] for i in range(0,len(kepid)) if kepid[i] != 11030475]
logg = [logg[i] for i in range(0,len(kepid)) if kepid[i] != 11030475]
feh = [feh[i] for i in range(0,len(kepid)) if kepid[i] != 11030475]
delta = [delta[i] for i in range(0,len(kepid)) if kepid[i] != 11030475]
derror = [derror[i] for i in range(0,len(kepid)) if kepid[i] != 11030475]
kepid = [i for i in kepid if i != 11030475] # there is no Gaia parallax for this KepID.

print(len(kepid))


from isochrones import StarModel
from isochrones.mist import MIST_Isochrone

mist = MIST_Isochrone()

meanrad = []


# =============================================================================
# print(teff[2], logg[2], feh[2], parallax[2], perror[2])
# model = StarModel(mist, Teff=(teff[2],200), logg=(logg[2],.03), feh=(feh[2],.04), parallax = (parallax[2],perror[2]))
# StarModel.fit_multinest(model,n_live_points=1000, basename=None, verbose=True, refit=True, overwrite=True,
# test=False)
# print(model.samples.radius_0_0.quantile(0.5))
# 
# =============================================================================
for i in range(0,len(kepid)):
    model = StarModel(mist, Teff=(teff[i],200), logg=(logg[i],.03), feh=(feh[i],.04), parallax = (parallax2[i],perror2[i]))
    model.fit_multinest(n_live_points=1000, basename=None, verbose=False, refit=True, overwrite=True,test=False)
    meanrad.append(model.samples.radius_0_0.quantile(0.5))
    print(i)

pradii = []
for i in range(len(delta)):
    pradii.append((meanrad[i]*sunrad*np.sqrt(delta[i]/1000000))/earthrad)
        
# =============================================================================
# file = open('USPHostProperties.csv', 'w')
# with file:
#     writer = csv.writer(file)
#     header = ['KepID', 'KOI', 'a/R','Period','Teff', 'logg', 'feh','delta(ppm)','delta error', 'parallax', 'parallax error', '50% star radius', 'planet radius']
#     writer.writerow(header)
#     for k in range(0,len(kepid)):
#         row = [kepid[k], koi[k], adivr[k], per[k],teff[k], logg[k], feh[k], delta[k], derror[k],  parallax2[k], perror2[k], meanrad[k],pradii[k]]
#         writer.writerow(row)
#         print(k)
# =============================================================================
        

