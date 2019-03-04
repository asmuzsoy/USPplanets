#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 10:48:42 2018

@author: anasofiauzsoy
"""
from isochrones import StarModel
from isochrones.mist import MIST_Isochrone
import csv
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

plt.rcParams['axes.linewidth']=3
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['xtick.minor.width'] = 2
plt.rcParams['ytick.minor.width'] = 2
plt.rc('xtick.major', size=8, pad=8)
plt.rc('xtick.minor', size=6, pad=5)
plt.rc('ytick.major', size=8, pad=8)
plt.rc('ytick.minor', size=6, pad=5)

plt.style.use('dark_background')
# mpl.rcParams.update(mpl.rcParamsDefault)

Earthradius = 6371000 # in meters
Sunradius = 695700000 # in meters 


## Calculating radii

kepid = []
koi = []
adivr = []
per = []
teff = []
logg = []
feh = []
deltas = []
derror = []
parallax = []
perror = []
starrad = []
pradii = []

with open('USPHostProperties.csv') as File:
     reader = csv.reader(File, delimiter=',')
     rownum = 0
     for row in reader:
         if rownum > 0:
             kepid.append(int(row[0]))
             koi.append(int(row[1]))
             adivr.append(float(row[2]))
             per.append(float(row[3]))
             teff.append(float(row[4]))
             logg.append(float(row[5]))
             feh.append(float(row[6]))
             deltas.append(float(row[7]))
             derror.append(float(row[8]))
             parallax.append(float(row[9]))
             perror.append(float(row[10]))
             starrad.append(float(row[11]))
             pradii.append(float(row[12]))
         rownum += 1

plt.figure(1)
#n = plt.hist(pradii, 10, color = 'xkcd:sky blue')
x = plt.xlabel('Planet Radius ($R_{ \oplus}$)', fontsize = '13')
y = plt.ylabel('Number of planets', fontsize = '13')
#plt.title('USP Planet Radii')

# =============================================================================
# ## Uncertainty calculations
# array = []
# array2 = []
# radun = []
# weights = []
# num_samples = 1000
# mist = MIST_Isochrone()
# # =============================================================================
# # 
# # file = open('PlanetRadiiSamples.csv', 'w')
# # with file:
# #     writer = csv.writer(file)
# # =============================================================================
#     
# for i in range(len(starrad)):
#     radiisample = []
#     model = StarModel(mist, Teff=(teff[i],200), logg=(logg[i],.03), feh=(feh[i],.04), parallax = (parallax[i],perror[i]))
#     model.fit_multinest(n_live_points=1000, basename=None, verbose=True, refit=True, overwrite=True,
# test=False)
#     print(i)
#     array = np.random.choice(model.samples.radius_0_0, num_samples)
#     array2 = np.random.normal(deltas[i],derror[i], num_samples)
#     print(len(array),len(array2))
#     for j in range(len(array)):
#         if array2[j] > 0:
#             radun.append((array[j]*Sunradius*np.sqrt(array2[j]/1000000))/Earthradius)
#             radiisample.append((array[j]*Sunradius*np.sqrt(array2[j]/1000000))/Earthradius)
#             weights.append(1/num_samples)
#     print(len(radiisample))
# # =============================================================================
# #         row = [kepid[i],radiisample]
# #         writer.writerow(row)
# # =============================================================================
# =============================================================================

periodsamples = []
radiisamples = []
weights = []

with open('PlanetRadiiSamples.csv') as File:
     reader = csv.reader(File, delimiter=',')
     rownum = 0
     for row in reader:
         kepid2 = (int(row[0]))
         samplelist = list(str(row[1]).split(','))
         samplelist[0] = samplelist[0][1:] # remove brackets
         samplelist[len(samplelist)-1] = samplelist[0][:-1] #remove brackets
         sampleslist = [float(i) for i in samplelist]
         for i in range(len(samplelist)):
             periodsamples.append(float(per[kepid.index(kepid2)]))
         radiisamples = radiisamples + sampleslist
         rownum += 1
weights = [0.001 for i in radiisamples]
            
plt.figure(2)
n = plt.hist(radiisamples, 10, weights = weights, color = 'xkcd:sky blue')
plt.xlabel('Planet Radius ($R_{ \oplus}$)', fontsize = '13')
plt.ylabel('Number of planets', fontsize = '13')
#plt.title('USP Planet Radii (with uncertanties)')

## plotting isochrones vs CKS data
radii2 = []
koi2 = []
teff2 = []
terror = []
feh2 = []
feherror = []
parallax2 = []
perror2 = []
rerrorup = []
rerrordown =[]
with open('tab_star-machine.csv') as File:
     reader = csv.reader(File, delimiter=',')
     rownum = 0
     for row in reader:
         if rownum > 0:
             if row[0] != "":
                 if int(row[0][2:]) in koi:
                         koi2.append(int(row[0][2:]))
                         teff2.append(float(row[1]))
                         terror.append(float(row[2]))
                         feh2.append(float(row[4]))
                         feherror.append(float(row[5]))
                         radii2.append(float(row[13]))
                         rerrorup.append(float(row[14]))
                         rerrordown.append(float(row[15]))
                         parallax2.append(float(row[9]))
                         perror2.append(float(row[10]))
         rownum += 1
         
radii3 = []
teff3 = []
feh3 = []
terror2 = []
feherror2 = []
rerrorup2 = []
rerrordown2 = []
koi3 = []
for i in range(len(radii2)):
    for j in range(len(koi)):
        if koi2[i] == koi[j]:
            koi3.append(koi[j])
            radii3.append(starrad[j])
            feh3.append(feh[j])
            teff3.append(teff[j])
            terror2.append(terror[i]) #i, not j!
            feherror2.append(feherror[i])
            rerrorup2.append(rerrorup[i])
            rerrordown2.append(rerrordown[i])
           
x = np.linspace(0.5,1.7)
plt.figure(3)
plt.plot(radii3, radii2, 'o')
plt.plot(x,x)
plt.xlabel('isochrones radii')
plt.ylabel('new radii')
plt.title('isochrones radii vs. new radii')

x = np.linspace(0,6500)
plt.figure(4)
plt.plot(teff3,teff2, '.')
plt.errorbar(teff3, teff2,yerr = terror2, linestyle = 'None', ecolor = 'k')
plt.plot(x,x)
plt.xlabel('Winn et al Teff')
plt.ylabel('Petigure/Fulton Teff')
plt.title('Teff comparison')

x = np.linspace(-1,0.7)
plt.figure(5)
plt.plot(feh3,feh2,'.')
plt.errorbar(feh3, feh2,yerr = feherror2, linestyle = 'None', ecolor = 'k')
plt.plot(x,x)
plt.xlabel('Winn et al Fe/H')
plt.ylabel('Petigure/Fulton Fe/H')
plt.title('Fe/H comparison')

plt.figure(6)
n = plt.hist(radii2, 10)
plt.xlabel('Planet Radius (Rearth)')
plt.ylabel('frequency')
plt.title('CKS VII Radius distribution')
plt.show()


kepid4 = []
rad2 = []
raderrup = []
raderrdown = []
with open('USPCandidates.csv', encoding='latin-1') as File:
    reader = csv.reader(File, delimiter=',')
    rownum = 0
    for row in reader:
        if rownum > 0:
            if row[0] != "":
                kepid4.append(int(row[0]))
                rad2.append(float(row[23]))
                raderrup.append(float(row[24]))
                raderrdown.append(float(row[25]))
        rownum += 1 
        
        
rad3 = []
raderrup2 = []
raderrdown2 = []
for i in range(len(kepid)):
    for j in range(len(kepid4)):
        if kepid[i] == kepid4[j]:
            rad3.append(rad2[j])
            raderrup2.append(raderrup[j])
            raderrdown2.append(raderrdown[j])
            
errors3 = np.stack((raderrup2,raderrdown2))
            
plt.figure(7)
x = np.linspace(0.5,2.2)
plt.plot(pradii,rad3, '.')
#plt.errorbar(pradii,rad3,xerr = errors3,linestyle = 'None', ecolor = 'k')
plt.plot(x,x)
plt.xticks(fontsize = 14, fontweight = 'bold')
plt.yticks(fontsize = 14, fontweight = 'bold')
plt.xlabel('Old planetary radii ($R_{ \oplus}$) (Sanchis-Ojeda et al., 2014)', fontsize = 16, fontweight = 'bold')
plt.ylabel('New planetary radii ($R_{ \oplus}$)', fontsize = 16, fontweight = 'bold')
plt.show()
            
            
        
plt.show()