from astropy.io import fits
from glob import os #para windows
import glob #para linux
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#############################################################

# Ordenar (por nombre) los archivos en una lista para posterior consulta
pathToFiles = "../DATA"
orderedList = sorted(glob.glob(os.path.join(pathToFiles, '*.fts')))
contador=0
for filename in orderedList:
    contador+=1
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
maph=pd.DataFrame(columns=headerList)

# fill Dataframe
cel =[]
for i in range(40):
    hdulist1 = fits.open(listaDatos[i])
    headerList = hdulist1[0].header
    hdulist1.close()
    for headerInplace in headerList:
        cel.append(headerList[headerInplace])
    maph.loc[i] = cel
    cel=[]

#############################################################

# Filtro por tiempo de exsposicion
tama = len(listaDatos)
expTime = maph.EXPTIME
mExpTime = np.mean(expTime)
maph_tem = maph[maph.EXPTIME>=mExpTime-1.5][maph.EXPTIME<=mExpTime+1.5]

# Filtro con Desviacion estandard
stdExpTime = np.std(maph.EXPTIME)
stdExpTime
maph_tem = maph[(maph.EXPTIME>=(mExpTime-stdExpTime)) & (maph.EXPTIME<=(mExpTime+stdExpTime))]
maph_tem.info()
