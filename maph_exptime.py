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
