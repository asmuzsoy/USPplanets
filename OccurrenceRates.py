#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 13:06:46 2018

@author: anasofiauzsoy
"""
import csv
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl


## makes plot look nice for presentations/papers
plt.rcParams['axes.linewidth']=3
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['xtick.minor.width'] = 2
plt.rcParams['ytick.minor.width'] = 2
plt.rc('xtick.major', size=8, pad=8)
plt.rc('xtick.minor', size=6, pad=5)
plt.rc('ytick.major', size=8, pad=8)
plt.rc('ytick.minor', size=6, pad=5)

#plt.style.use('dark_background')
mpl.rcParams.update(mpl.rcParamsDefault)

kepid = []
koi = []
adivr = []
per = []
starrad = []
pradii = []
numstars = []

radbins = [0.5,(1/np.sqrt(2)),1,np.sqrt(2),2,2*np.sqrt(2),4]
radbinsflip = np.flip(radbins,0) ## makes orientation correct on y axis

perbins = [4,4*(6**0.25),4*np.sqrt(6),4*(6**0.75),24]
array = np.zeros((len(radbins)-1,len(perbins)-1))
scaledarray = np.zeros((len(radbins)-1,len(perbins)-1))
freqarray = np.zeros((len(radbins)-1,len(perbins)-1))


with open('USPHostProperties.csv') as File:
     reader = csv.reader(File, delimiter=',')
     rownum = 0
     for row in reader:
         if rownum > 0:
             kepid.append(int(row[0]))
             koi.append(int(row[1]))
             adivr.append(float(row[2]))
             per.append(float(row[3]))
             starrad.append(float(row[11]))
             pradii.append(float(row[12]))
             numstars.append(float(row[13]))
         rownum += 1
        
        
## calculates occurrence rate contribution for each planet sample
f = [(adivr[i]/numstars[i]) for i in range(0,len(kepid))]

kepid2 = []

periodsamples = []
radiisamples = []

## get 1000 planetary radius samples per planet, set up period array
with open('PhotometryPlanetRadiiSamples.csv') as File:
#with open('PlanetRadiiSamples.csv') as File:
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

# =============================================================================
# print(len(radiisamples))
# print(len(periodsamples))
# =============================================================================

## sum occurrence rate in bins
for i in range(len(periodsamples)):
    for j in range(len(perbins)-1):
        for k in range(len(radbins)-1):
            if (perbins[j] <= 24*periodsamples[i] <= perbins[j+1]) and (radbinsflip[k] >= radiisamples[i] >= radbinsflip[k+1]):
                array[k][j] += (f[per.index(periodsamples[i])])/1000
                freqarray[k][j] += 1
        
# =============================================================================
# count = 0        
# for i in range(len(pradii)):
#     if  0.5  <= pradii[i] <= 1/np.sqrt(2):
#         count += 1
# =============================================================================
        

#print(f)
print(array)
print(freqarray)


## scales array for logarithmic scale
for j in range(0,len(perbins)-1):
    for i in range(0,len(radbins)-1):
        if array[i][j] != 0:
            scaledarray[i][j] = 1 + (np.log10(array[i][j]/(5*(10**-9)))) # moved up by 1 to make way for zero.
        if array[i][j] == 0:
            scaledarray[i][j] = 1 + -1
        
print(scaledarray)

## attempt at shading bottom row, doesn't work
grayarray = np.zeros((len(radbins)-1,len(perbins)-1))
grayarray[:1] = [0.1,0.1,0.1,0.1]
grayarray = np.flip(grayarray,0)
#print(grayarray)

## convert period to hours
per2 = [24 * i for i in per]

## plot with larger bins, linear scale
fig, ax = plt.subplots()
ax = sns.heatmap(1000*array, annot = True, annot_kws={"size": 14})
ax.set_xticks([0,1,2,3,4])
ax.set_xticklabels([4,6.3,9.8,15.3,24], fontsize = 14, fontweight = 'bold')
ax.set_yticks([0,1,2,3,4,5,6])
ax.set_yticklabels([0.5,0.707,1,1.414,2,2.83,4], fontsize = 14, fontweight = 'bold')
cbar = ax.collections[0].colorbar
cbar.set_label('Planet Occurence - # planets per 1000 stars', fontsize = 14, fontweight = 'bold')
plt.xlabel("Orbital Period (days)", fontsize = 16, fontweight = 'bold')
plt.ylabel("Radius (Earth radii)", fontsize = 16, fontweight = 'bold')

plt.show()

## plot with larger bins, logarithmic scale
fig, ax = plt.subplots()
ax = sns.heatmap(scaledarray)
ax.set_xticks([0,1,2,3,4])
ax.set_xticklabels([4,6.3,9.8,15.3,24])
ax.set_yticks([0,1,2,3,4,5,6])
ax.set_yticklabels([0.5,0.707,1,1.414,2,2.83,4])
cbar = ax.collections[0].colorbar
cbar.set_ticks([-1,0,1,2,3,4,5,6,7])
cbar.set_ticklabels(['0', '5e-9', '5e-8', '5e-7', '5e-6', '5e-5', '5e-4', '0.005'])
cbar.set_label('Planet Occurence - fcell')
plt.xlabel("Orbital Period (days)")
plt.ylabel("Radius (Earth radii)")

## list of smaller radius bins
radbins2 = [0.5,0.5*(2**0.25),(1/np.sqrt(2)),(1/np.sqrt(2)) * (2**0.25),1,(2**0.25),np.sqrt(2),(2**0.75), 2,2*(2**0.25),2*np.sqrt(2),2*(2**0.75),4]
radbins2flip = np.flip(radbins2,0)
#print(radbins2flip)

array2 = np.zeros((len(radbins2)-1,len(perbins)-1))

for i in range(len(periodsamples)):
    for j in range(len(perbins)-1):
        for k in range(len(radbins2)-1):
            if (perbins[j] <= 24*periodsamples[i] <= perbins[j+1]) and (radbins2flip[k] >= radiisamples[i] >= radbins2flip[k+1]):
                array2[k][j] += (f[per.index(periodsamples[i])])/1000
                
scaledarray2 = np.zeros((len(radbins2)-1,len(perbins)-1))

for j in range(0,len(perbins)-1):
    for i in range(0,len(radbins2)-1):
        if array2[i][j] != 0:
            scaledarray2[i][j] = 1 + (np.log10(array2[i][j]/(5*(10**-9))))
        if array2[i][j] == 0:
            scaledarray2[i][j] = 1 + -1 # moved up by 1 to make way for zero.
            
#print(array2)
               
## plot with smaller bins, linear scale            
fig, ax = plt.subplots()
ax = sns.heatmap(array2,annot = True, annot_kws={"size": 15})
ax.set_xticks([0,1,2,3,4])
ax.set_xticklabels([4,6.3,9.8,15.3,24])
ax.set_yticks([0,1,2,3,4,5,6,7,8,9,10,11,12])
ax.set_yticklabels([0.5,0.59,0.70,0.84,1,1.19,1.41,1.68,2,2.38,2.83,3.36,4])
cbar = ax.collections[0].colorbar
cbar.set_label('Planet Occurence - fcell')
plt.xlabel("Orbital Period (days)")
plt.ylabel("Radius (Earth radii)")

## plot with smaller bins, logarithmic scale
fig, ax = plt.subplots()
ax = sns.heatmap(scaledarray2)
ax.set_xticks([0,1,2,3,4])
ax.set_xticklabels([4,6.3,9.8,15.3,24])
ax.set_yticks([0,1,2,3,4,5,6,7,8,9,10,11,12])
ax.set_yticklabels([0.5,0.59,0.70,0.84,1,1.19,1.41,1.68,2,2.38,2.83,3.36,4])
cbar = ax.collections[0].colorbar
cbar.set_ticks([-1, 0, 1, 2, 3, 4,  5,6])
cbar.set_ticklabels(['0', '5e-9', '5e-8', '5e-7', '5e-6', '5e-5', '5e-4'])
cbar.set_label('Planet Occurence - fcell', fontsize = 16)
plt.xlabel("Orbital Period (days)", fontsize = 16)
plt.ylabel("Radius (Earth radii)", fontsize = 16)

plt.show()

persum = np.zeros(len(perbins)-1)
## sum by period bins
for i in range(len(radbins) - 1):
    for j in range(len(perbins)-1):
        persum[j] += array[i][j]*1000
        
#print(persum)

radsum = np.zeros(len(radbins2)-1)
## sum by radius bins
for i in range(len(radbins2) - 1):
    for j in range(len(perbins)-1):
        radsum[i] += array2[i][j]*1000
        
#print(radsum)
        
## data from Sanchis-Ojeda et al. (2014) Figure 9
paperperx = [5.019858,7.8650255, 12.217099,19.226887]
paperpery = [0.11733242,0.369,1.1060447,4.2267613]
papererrorup = [0.195-0.11733242,0.486-0.369,1.34-1.1060447,4.97-4.2267613]
papererrordown = [0.11733242-0.053,0.369-0.248,1.1060447-0.84,4.2267613-3.59]
errors = np.stack((papererrorup,papererrordown))

## period vs. occurrence rate plot
plt.loglog(paperperx, paperpery, label = 'Sanchis-Ojeda et al. (2014)', marker = 'o')
plt.errorbar(paperperx, paperpery,yerr = errors, linestyle = 'None', ecolor = 'k')
plt.plot(paperperx,persum, label = 'me', marker = 'o')
plt.xlabel('Orbital period (hours)')
plt.ylabel('Number of planets per thousand stars')
plt.legend()
plt.show()

paperradx = [0.9144446,1.0894922,1.29322,1.5407747,1.8288893, 2.1708794, 2.5960968,3.0815494,3.671435]
paperrady = [1.5357143,1.3333334,1.3035715,0.92261904,0.42857143,0,0,0,0]
raderrorup = [2.0130641-1.5357143,1.6983373-1.3333334,1.58711-1.3035715,1.1401426-0.92261904,0.5700713-0.42857143,0.149165,0.149165,0.149165,0.149165]
raderrordown = [1.5357143-1.0510689,1.3333334-0.9679335,1.3035715-1.0095012,0.92261904-0.70665085,0.42857143-0.2672209,0,0,0,0]
errors2 = np.stack((raderrorup,raderrordown))

## radius vs. occurrence rate plot
plt.loglog(paperradx, paperrady, label = 'Sanchis-Ojeda et al. (2014)', marker = 'o')
plt.errorbar(paperradx, paperrady,yerr = errors2, linestyle = 'None', ecolor = 'k')
plt.plot(paperradx, np.flip(radsum[:-3],0), label = 'me', marker = 'o')
plt.xlabel('Planet radius (Rearth')
plt.ylabel('Number of planets per thousand stars')
plt.legend()
plt.show()

## calculating transit occurrence rate for uncertainty calculation
# =============================================================================
# f2 = np.zeros(len(kepid))
# for i in range(len(kepid)):
#     if numstars[i] != 0:
#         f2[i] = 1/numstars[i]
#     else:
#         f2[i]= 0
# =============================================================================
f2 = [(1/numstars[i]) for i in range(0,len(kepid))]

transitarray = np.zeros((len(radbins)-1,len(perbins)-1))

for i in range(len(periodsamples)):
    for j in range(len(perbins)-1):
        for k in range(len(radbins)-1):
            if (perbins[j] <= 24*periodsamples[i] <= perbins[j+1]) and (radbinsflip[k] >= radiisamples[i] >= radbinsflip[k+1]):
                transitarray[k][j] += (f2[per.index(periodsamples[i])])/1000


# =============================================================================
# ## model as binomial distribution, take 1 sigma (like Howard et al 2012)
# uncertainties = np.empty((len(radbins)-1,len(perbins)-1,2),dtype = object)
# for i in range(len(radbins)-1):
#     for j in range(len(perbins)-1):
#         if array[i][j] != 0:
#             n = (freqarray[i][j])/(transitarray[i][j])
#             binom = np.random.binomial(n,transitarray[i][j],1000)
#             t = ((np.percentile(binom,15.9)),np.percentile(binom,84.1))
#             uncertainties[i][j] = "(" + str((np.percentile(binom,15.9)/n)) +"," + str(np.percentile(binom,84.1)/n) + ")"
# 
# print(uncertainties)
# =============================================================================

## calculate 68.3% Wilson Confidence Interval
wilsonhigh = np.empty((len(radbins)-1,len(perbins)-1), dtype = object)
wilsonlow = np.empty((len(radbins)-1,len(perbins)-1), dtype = object)
wilson3d = np.empty((len(radbins)-1,len(perbins)-1,2))
for i in range(len(radbins)-1):
    for j in range(len(perbins)-1):
        if transitarray[i][j] == 0:
            wilsonhigh[i][j] = "N=0"
            wilsonlow[i][j] = "N=0"
        else:
            n = (freqarray[i][j])/(transitarray[i][j])
            if n == 0:
                wilsonhigh[i][j] = "n=0"
                wilsonlow[i][j] = "n=0"
            else:
                z = 0.41
                high = ((transitarray[i][j] + (z**2)/(2*n))/(1+(z**2)/n)) + (z/(1+((z**2)/n)))*np.sqrt(((transitarray[i][j]*(1-transitarray[i][j]))/n)-((z**2)/(4*(n**2))))
                low = ((transitarray[i][j] + (z**2)/(2*n))/(1+(z**2)/n)) - (z/(1+((z**2)/n)))*np.sqrt(((transitarray[i][j]*(1-transitarray [i][j]))/n)-((z**2)/(4*(n**2))))
                wilsonhigh[i][j] = (1000 * (array[i][j]/transitarray[i][j])*high) - 1000*array[i][j]
                wilsonlow[i][j] = (1000 * (array[i][j]/transitarray[i][j])*low) - 1000*array[i][j]
                wilson3d[i][j] = [low,high]

print("high")
print(wilsonhigh)
print("low")
print(wilsonlow)

labels = np.array([[0,0,0.0012,0.0084],[0,0,0.0444,0.0786],[0.0706,0.0564,0.4258,1.1371],[0.0517,0.3323,0.4635,1.6564],[0,0.0489,0.0673,1.1179],[0,0,0,0.0425]],dtype = object)
fig, ax = plt.subplots()
ax = sns.heatmap(1000*array, annot = labels, fmt = '', annot_kws={"size": 14})
plt.text(2.21,5.75,"+0.00010", size = 8)
plt.text(2.21,5.16," -0.00009", size = 8)
plt.text(3.21,5.75,"+0.00040", size = 8)
plt.text(3.21,5.16," -0.00038", size = 8)
plt.text(2.21,4.75,"+0.00059", size = 8)
plt.text(2.21,4.16," -0.00058", size = 8)
plt.text(3.21,4.75,"+0.00115", size = 8)
plt.text(3.21,4.16," -0.00113", size = 8)
plt.text(0.21,3.75,"+0.00094", size = 8)
plt.text(0.21,3.16," -0.00093", size = 8)
plt.text(1.21,3.75,"+0.00053", size = 8)
plt.text(1.21,3.16," -0.00052", size = 8)
plt.text(2.21,3.75,"+0.00186", size = 8)
plt.text(2.21,3.16," -0.00185", size = 8)
plt.text(3.21,3.75,"+0.00349", color = 'w', size = 8)
plt.text(3.21,3.16," -0.00348", color = 'w', size = 8)
plt.text(0.21,2.75,"+0.00066", size = 8)
plt.text(0.21,2.16," -0.00065", size = 8)
plt.text(1.21,2.75,"+0.00170", size = 8)
plt.text(1.21,2.16," -0.00169", size = 8)
plt.text(2.21,2.75,"+0.00209", size = 8)
plt.text(2.21,2.16," -0.00208", size = 8)
plt.text(3.21,2.75,"+0.00512", color = 'w', size = 8)
plt.text(3.21,2.16," -0.00511", color = 'w', size = 8)
plt.text(1.21,1.75,"+0.00081", size = 8)
plt.text(1.21,1.16," -0.00079", size = 8)
plt.text(2.21,1.75,"+0.00097", size = 8)
plt.text(2.21,1.16," -0.00096", size = 8)
plt.text(3.21,1.75,"+0.00687", color = 'w', size = 8)
plt.text(3.21,1.16," -0.00683", color = 'w', size = 8)
plt.text(3.21,0.75,"+0.00156", size = 8)
plt.text(3.21,0.16," -0.00151", size = 8)
ax.set_xticks([0,1,2,3,4])
ax.set_xticklabels([4,6.3,9.8,15.3,24], fontsize = 14)
ax.set_yticks([0,1,2,3,4,5,6])
ax.set_yticklabels([0.5," 0.707",1,1.414,2,2.83,4], fontsize = 14, horizontalalignment = 'right', verticalalignment = 'baseline')
cbar = ax.collections[0].colorbar
cbar.set_label('Planet Occurence - # planets / 1000 stars', fontsize = 14)
plt.xlabel("Orbital Period (hours)", fontsize = 16)
plt.ylabel("Planet Radius ($R_{\oplus}$)", fontsize = 16)

#ax2 = sns.heatmap(wilsonhigh, annot=True, alpha=0.0)

plt.show()



error4 = []
meanradband = []

for i in range(0,71):
    list1 = radiisamples[1000*i:(1000*(i+1))-1]
    meanradband.append(np.percentile(list1,50))
    error4.append(np.std(list1))


## plot median radii vs. period
plt.loglog(per2,meanradband, 'o', basex = (6**0.25), basey = np.sqrt(2))
plt.errorbar(per2,meanradband, yerr = error4, linestyle = 'None', ecolor = 'k')
plt.yticks(radbins, [0.5,0.707,1,1.414,2,2.83,4])
plt.xticks(perbins,[4,6.3,9.8,15.3,24])
plt.xlabel("Orbital Period (hours)")
plt.ylabel("Planet Radius ($R_{ \oplus}$)")
#plt.title("Photometric band radii")
plt.show()

periodsamples2 = [24*i for i in periodsamples]

# =============================================================================
# ## plot radii samples vs. period
# plt.loglog(periodsamples2,radiisamples, 'o', basex = (6**0.25), basey = np.sqrt(2), alpha=0.25)
# plt.yticks(radbins, [0.5,0.707,1,1.414,2,2.83,4], fontsize = 14, fontweight = 'bold')
# plt.xticks(perbins,[4,6.3,9.8,15.3,24], fontsize = 14, fontweight = 'bold')
# plt.xlabel("Orbital Period (hours)", fontsize = 16, fontweight = 'bold')
# plt.ylabel("Planet Radius ($R_{ \oplus}$)",fontsize = 16, fontweight = 'bold')
# plt.show()
# 
# =============================================================================
periodsamples4 = []
radiisamples4 = []


## old spectroscopic properties radii
with open('PlanetRadiiSamples.csv') as File:
     reader = csv.reader(File, delimiter=',')
     rownum = 0
     for row in reader:
         #kepid2 = (int(row[0]))
         samplelist = list(str(row[1]).split(','))
         samplelist[0] = samplelist[0][1:] # remove brackets
         samplelist[len(samplelist)-1] = samplelist[0][:-1] #remove brackets
         sampleslist = [float(i) for i in samplelist]
         for i in range(len(samplelist)):
             periodsamples4.append(float(per[kepid.index(kepid2)]))
         radiisamples4 = radiisamples4 + sampleslist
         rownum += 1
         
error5 = []
meanrad = []
for i in range(0,71):
    list1 = radiisamples4[1000*i:(1000*(i+1))-1]
    meanrad.append(np.percentile(list1,50))
    error5.append(np.std(list1))
    
    
plt.loglog(per2,meanrad, 'o', basex = (6**0.25), basey = np.sqrt(2))
plt.errorbar(per2,meanrad, yerr = error5, linestyle = 'None', ecolor = 'k')
plt.yticks(radbins, [0.5,0.707,1,1.414,2,2.83,4])
plt.xticks(perbins,[4,6.3,9.8,15.3,24])
plt.xlabel("Orbital Period (hours)")
plt.ylabel("Planet Radius ($R_{ \oplus}$)")
plt.title("Spectroscopic Parameter radii")
plt.show()

koiwCKS = [72, 191, 577, 1128, 1150, 1169, 1300, 1360, 1367, 1428, 1442, 1655, 1688, 1875, 2079, 2119, 2202, 2250, 2281, 2393, 2396, 2409, 2492, 2517, 2571, 2607, 2668, 2694, 2753, 2756, 2763, 2874, 2875, 2916, 3009, 3032, 3065, 3867, 4002, 4018, 4070, 4144, 4366, 4430, 4441]

radiiwCKS = [meanradband[koi.index(i)] for i in koiwCKS ]
radiibuds = [meanrad[koi.index(i)] for i in koiwCKS]
otherradiiband = [i for i in meanradband if i not in radiiwCKS]
otherradii = [i for i in meanrad if i not in radiibuds]



x = np.linspace(0.5,2.2)   
plt.plot(radiibuds,radiiwCKS,'bo', label = 'photometric bands + CKS spectroscopic data' )
plt.plot(otherradii,otherradiiband, 'ro', label = 'only photometric bands')
plt.errorbar(meanrad,meanradband, xerr = error5,yerr = error4 , linestyle = 'None', ecolor = 'k')
plt.plot(x,x, color = 'y')
plt.xlabel("Old spectroscopic properties radii")
plt.ylabel("New photometric band + CKS spectroscopic for some radii")
plt.title("Planetary radii comparison")
plt.legend()
plt.show()

radii7 = np.zeros(len(kepid)) ## Sanchis-Ojeda radii
radii7errup = np.zeros(len(kepid))
radii7errdown = np.zeros(len(kepid))
with open('USPCandidates.csv', encoding='latin-1') as File:
     reader = csv.reader(File, delimiter=',')
     rownum = 0
     for row in reader:
         if rownum > 0:
             if row[24] != "":
                 for i in range(len(kepid)):
                     if int(row[0]) == kepid[i]:
                         radii7[i] = (float(row[9]))
                         radii7errup[i] = (float(row[10]))
                         radii7errdown[i] = (float(row[11]))

         rownum += 1

magradii = [1.4185062355258975, 1.0128577290505476, 0.7395642185973244, 0.8718636643528979, 0.74277899440982, 0.8251688587227919, 0.9170365400567101, 0.8110354169337991, 0.8930698076902541, 0.9656160980644677, 0.6451729946227232, 0.6722056366700062, 0.6613939035363536, 0.8073935123830203, 0.8609784068563593, 0.9565993670405968, 0.9975652542760644, 0.6634612906668464, 0.8821534748085418, 0.6944669394891578, 0.785799640462931, 0.9403186141664807, 0.8228863441088534, 1.3019868811668818, 1.084489070507936, 1.6433495565444605, 0.8555009121304379, 1.2408125387260347, 0.736661982870762, 1.0571785991048095, 0.8734215615118575, 0.7272064197247533, 0.7748439544289258, 0.8459174731212303, 1.4157487718277189, 0.8439722448368516, 2.0391469252439958, 0.9609444461851793, 0.978254005193265, 0.8333964607833079, 0.8031435869064345, 0.6450487437816237, 0.8032472970998629, 0.9845700567987361, 0.757230205369519, 0.69284633855208, 0.8775549260118853, 1.055835849899823, 0.8372480183762283, 0.8257959562932159, 0.75739153390834, 0.7111707406802878, 0.9179480853386884, 1.0510167557101537, 1.090912498238495, 0.6974564720552092, 1.0571629176381634, 0.723617581610305, 0.7741279941280357, 0.8027130898506116, 0.7242784868126759, 0.8398411202169214, 0.6270665839279399, 1.3385064572601753, 0.7807403678647473, 0.7638964778501294, 0.9560069029531193, 0.7484923725216783, 0.9745764735768543, 0.7652803094892944, 0.8230487211415367]
sigma = [0.04987787443578625, 0.021119636599107934, 0.012047172567476776, 0.01645788407276533, 0.012641365590110592, 0.01671466913214443, 0.03809329439315205, 0.01993302814025634, 0.025758640790493566, 0.02744681391087703, 0.009463183143611704, 0.01590901087869035, 0.013913315043218601, 0.011749271062555008, 0.012286352779378278, 0.017909973365102635, 0.020589546027321545, 0.015908332630537057, 0.015211129336216971, 0.012956466488615679, 0.01414280405744474, 0.02227407026852073, 0.03115468166255689, 0.029342425084388037, 0.0793925837673903, 0.05174225294814394, 0.013575758243414601, 0.03553641934908235, 0.010159923525233318, 0.027229004196906194, 0.01678110319081692, 0.012610666032096832, 0.02159333657505373, 0.02660606230120436, 0.1525149160724844, 0.012975708030965523, 0.33850938847599255, 0.025219961625041365, 0.04403168081256531, 0.01583082955202924, 0.010800629901046159, 0.010820868008239231, 0.024301494993612997, 0.015946416724393227, 0.014475983686927235, 0.009409240215827074, 0.020276830671083088, 0.01724794548302806, 0.023309676331285938, 0.019731986490652498, 0.02166577023378421, 0.00896886687109212, 0.01569408678343353, 0.01771524966216805, 0.021745453065290538, 0.018447733142068427, 0.016620428318208324, 0.009608187232514885, 0.018451627200353402, 0.012085877303140766, 0.011214980749553235, 0.015148069834151973, 0.009605941995319982, 0.050580184992134423, 0.011022025888956651, 0.010262424510071462, 0.02087163293932602, 0.010247520887608698, 0.01925581991613471, 0.009320682304407375, 0.012779118999560577]

starradiiwCKS = [magradii[koi.index(i)] for i in koiwCKS ]
radiibuds2 = [radii7[koi.index(i)] for i in koiwCKS]

y = np.linspace(0,2.5)
radii7error = np.stack((radii7errdown, radii7errup))
plt.plot(radii7,magradii,'bo', label = 'only photometric band magnitudes')
plt.plot(radiibuds2,starradiiwCKS,'ro', label = 'photometric band magnitudes + spectroscopic data')
plt.plot(y,y, color = 'tab:orange')
plt.xlim(0,2)
plt.ylim(0,2.5)
plt.errorbar(radii7,magradii,xerr = radii7error,yerr = sigma,linestyle = 'None', ecolor = 'k')
plt.xlabel("Sanchis-Ojeda et al. (2014) USP planet host radii (Solar radii)")
plt.ylabel("Stellar host radii derived using isochrones (Solar radii)")
plt.legend()
plt.show()

mortonkepid = np.zeros(len(kepid))
mortonstarradii = np.zeros(len(kepid))
mortonerrup = np.zeros(len(kepid))
mortonerrdown = np.zeros(len(kepid))
mortonror = np.zeros(len(kepid))

with open('q1_q17_dr25_koifpp.csv') as File: ## Morton radii
     reader = csv.reader(File, delimiter=',')
     rownum = 0
     for row in reader:
         if rownum > 39:
             if int(row[1]) in kepid:
                 mortonkepid[kepid.index(int(row[1]))] = int(row[1])
                 mortonstarradii[kepid.index(int(row[1]))] = float(row[21])
                 mortonerrup[kepid.index(int(row[1]))] = float(row[22])
                 mortonror[kepid.index(int(row[1]))] = float(row[17])
                 mortonerrdown[kepid.index(int(row[1]))] = -1*float(row[23])
                 
         rownum += 1


magradii = [1.4185062355258975, 1.0128577290505476, 0.7395642185973244, 0.8718636643528979, 0.74277899440982, 0.8251688587227919, 0.9170365400567101, 0.8110354169337991, 0.8930698076902541, 0.9656160980644677, 0.6451729946227232, 0.6722056366700062, 0.6613939035363536, 0.8073935123830203, 0.8609784068563593, 0.9565993670405968, 0.9975652542760644, 0.6634612906668464, 0.8821534748085418, 0.6944669394891578, 0.785799640462931, 0.9403186141664807, 0.8228863441088534, 1.3019868811668818, 1.084489070507936, 1.6433495565444605, 0.8555009121304379, 1.2408125387260347, 0.736661982870762, 1.0571785991048095, 0.8734215615118575, 0.7272064197247533, 0.7748439544289258, 0.8459174731212303, 1.4157487718277189, 0.8439722448368516, 2.0391469252439958, 0.9609444461851793, 0.978254005193265, 0.8333964607833079, 0.8031435869064345, 0.6450487437816237, 0.8032472970998629, 0.9845700567987361, 0.757230205369519, 0.69284633855208, 0.8775549260118853, 1.055835849899823, 0.8372480183762283, 0.8257959562932159, 0.75739153390834, 0.7111707406802878, 0.9179480853386884, 1.0510167557101537, 1.090912498238495, 0.6974564720552092, 1.0571629176381634, 0.723617581610305, 0.7741279941280357, 0.8027130898506116, 0.7242784868126759, 0.8398411202169214, 0.6270665839279399, 1.3385064572601753, 0.7807403678647473, 0.7638964778501294, 0.9560069029531193, 0.7484923725216783, 0.9745764735768543, 0.7652803094892944, 0.8230487211415367]


plotradmorton = [i for i in mortonstarradii if i != 0]
plotradmagradii = [magradii[i] for i in range(len(magradii)) if mortonstarradii[i] != 0]
mortonerrup = [mortonerrup[i] for i in range(len(magradii)) if mortonstarradii[i] != 0]
mortonerrdown = [mortonerrdown[i] for i in range(len(magradii)) if mortonstarradii[i] != 0]

        
mortonerror = np.stack((mortonerrup,mortonerrdown))
    
plt.plot(plotradmorton,plotradmagradii,'o')
plt.plot(y,y)
plt.errorbar(plotradmorton,plotradmagradii,xerr = mortonerror, linestyle = 'None', ecolor = 'k')
plt.xlabel("Morton stellar radii")
plt.ylabel("New photometric band + spectroscopy for some radii")
plt.title("Blue has bands + CKS spectroscopic data")
plt.show()

masserrors = [1.4193578747640816, 0.20527205368400664, 0.29621778996024711, 0.30218425229114243, 0.13570151082339552, 0.629029688580795, 0.82383449237583328, 0.22664492240356018, 0.64327471983283269, 0.41288129290256875, 0.22182848088679732, 0.32926508847148661, 0.4763733419252521, 0.21589131442167206, 0.16531524616011511, 0.41865113720372465, 1.0081895368513911, 0.74209240205613836, 0.36767774737901704, 0.28753267566765317, 0.75461158652858762, 0.40278010276485171, 0.8202176325734114, 0.22824734781818029, 2.040358261088199, 1.1564851316604248, 0.16063251481103091, 0.71959401199533479, 0.29249264212761955, 0.23117614830868852, 0.24415811189222669, 0.14777203456257865, 0.79764583541268774, 0.27128709593202288, 1.1345545914617301, 0.13096700333633279, 14.114844414860242, 0.23686772377593232, 1.9460346462741056, 0.18664911537548956, 0.059233273242532074, 0.31847046711360033, 0.67836293990591023, 0.029624017184777036, 0.75958563425587444, 0.29016837117999567, 0.819961597708459, 0.3535171321858282, 0.57623525098640993, 0.45396617336525824, 0.76676131647009182, 0.38384387987651147, 0.28550947316101599, 0.15746350988461827, 0.14276560675934752, 0.95584819751276173, 0.29498980025645427, 0.31508077702052434, 0.96098514285148595, 0.22754688766478467, 0.23281269168591245, 0.39382048286593097, 0.12879495380002673, 5.3309404036945534, 0.1285900824064419, 0.2330733696347892, 0.69487455466954529, 0.28339145512102554, 0.35756455394744702, 0.12607068989593062, 0.33681037918730472]
meanmasses = [4.629942823350425, 1.63903652361664, 2.6327540067735891, 2.7780555118610781, 2.1256904380084611, 2.9589460628542734, 3.9450958502115503, 1.5915625830006455, 3.1783036674852223, 1.8246642750422, 3.736882988661518, 1.9129674297445496, 3.8218084317273453, 3.2655890892246804, 1.9506286847120593, 1.9772598497337692, 6.7842691659708843, 6.1747565361669281, 4.584804786464785, 2.2440555626842325, 7.9999472680241448, 3.2405557621104095, 4.2154090839363789, 1.6060728955453456, 6.4756380572083492, 7.8656480677506559, 2.61568033747277, 5.0114198065738549, 4.2584900034388173, 1.6652949210358632, 0.70667385732457677, 0.70110570155302476, 5.1395374652916432, 1.4272713291353056, 2.6826460347009626, 0.89847691154103537, 24.136577294090294, 1.6137618182649756, 10.470769632986693, 1.729307889011976, 0.68754548375780333, 1.1250149648166885, 4.2259047270128294, 0.31002530478237755, 3.0762802472519701, 4.5143663691740086, 7.3959548626289786, 5.3671196203989266, 3.0568515620400554, 2.3983530598658591, 7.1037375148056237, 7.4709833289502425, 3.8199271701727255, 2.3982831548411667, 1.1215082384472697, 7.7786834806473113, 4.8795174392296161, 4.0771966285063508, 3.2635831494958696, 1.3555262507833654, 2.062666862162192, 3.6112935844403276, 1.0948691895933078, 34.04727262902616, 0.95396596018396673, 3.4031403206982858, 7.4298279493229131, 4.4821433128440766, 3.9095343771929918, 1.3621205646105849, 4.5561092633046512]

## prints Tex syntax for table
# =============================================================================
# for i in range(len(kepid)):
#     print(kepid[i],'&', koi[i], '&', round(magradii[i],3), '$\pm$', round(sigma[i],3), '&',round(meanradband[i],3), '$\pm$', round(error4[i],3),'&',round(meanmasses[i],3), '$\pm$' ,round(masserrors[i],3),'\\\\')
# 
# =============================================================================
