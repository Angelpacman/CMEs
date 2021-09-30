from astropy.io import fits
from glob import os  # para windows
import glob  # para linux
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#############################################################

# Ordenar (por nombre) los archivos en una lista para posterior consulta
pathToFiles = "../DATA"
#pathToFiles = "DATA"
orderedList = sorted(glob.glob(os.path.join(pathToFiles, '*.fts')))
contador = 0
for filename in orderedList:
    contador += 1
print("\nExisten {} archivos fts".format(contador))

#############################################################

# Obtener lista de datos que son 2048x2048
listaDatos = []
for filename in orderedList:
    hdulist = fits.open(filename)
    data = hdulist[0].data
    hdulist.close
    if data.shape == (2048, 2048):
        listaDatos.append(filename)
    else:
        break
print("\nDe los cuales {} son  de tamaño 2048 x 2048".format(len(listaDatos)))

#############################################################

# make Dataframe with headers
hdulist1 = fits.open(listaDatos[0])
headerList = hdulist1[0].header
hdulist1.close()
maph = pd.DataFrame(columns=headerList)

# fill Dataframe
cel = []
for i in range(40):
    hdulist1 = fits.open(listaDatos[i])
    headerList = hdulist1[0].header
    hdulist1.close()
    for headerInplace in headerList:
        cel.append(headerList[headerInplace])
    maph.loc[i] = cel
    cel = []

#############################################################

# Filtro por tiempo de exsposicion
tama = len(listaDatos)
expTime = maph.EXPTIME
mExpTime = np.mean(expTime)
maph_tem = maph[maph.EXPTIME >= mExpTime - 1.5][maph.EXPTIME <= mExpTime + 1.5]

# Filtro con standard deviation
stdExpTime = np.std(maph.EXPTIME)
maph_tem = maph[(maph.EXPTIME >= (mExpTime - stdExpTime)) & (maph.EXPTIME <= (mExpTime + stdExpTime))]
maph_tem.info()

#############################################################

mitad=1 #esto es para elegir el tamaño de la imagen 2 -> 512; 1 -> 1024
xc=maph.loc[0].CRPIX1/mitad
yc=maph.loc[0].CRPIX2/mitad

#############################################################

for i in maph:
    if i == 'R_SUN':
        print("R_SUN existe")
    elif i == 'RSUN':
        print("RSUN existe")

#############################################################

for header in maph:
    if header == 'R_SUN':
        r0=maph.loc[0].R_SUN/(maph.loc[0].CDELT1*mitad)
    elif header == 'RSUN':
        r0=maph.loc[0].RSUN/(maph.loc[0].CDELT1*mitad)

r0_20=20.*r0
r0_19=19.*r0
r0_15=15.*r0
r0_10=10.*r0

print("r0_20 = {}".format(r0_20))
print("r0_19 = {}".format(r0_19))
print("r0_15 = {}".format(r0_15))
print("r0_10 = {}".format(r0_10))

#############################################################

import datetime
#df_tiempo=pd.to_datetime(maph['DATE-OBS'])
tiempos = [maph['DATE-OBS']]
header = ["DATE-OBS"]
df_tiempo = pd.concat(tiempos, axis = 1, keys = header)
df_tiempo.head()

#############################################################

tiempos=pd.to_datetime(df_tiempo['DATE-OBS'])
t0 = tiempos[0] # Esto fallará si se filtro el FITfile con el indice 0
delta = []
for celda in tiempos:
    delta.append((celda - t0).total_seconds())
print("delta = {}".format(delta))

#############################################################

df_tiempo['SECONDS'] = delta
df_tiempo.head()
