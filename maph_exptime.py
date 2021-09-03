from astropy.io import fits
from glob import os #para windows
import glob #para linux
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# Ordenar (por nombre) los archivos en una lista para posterior consulta
pathToFiles = "DATA"
orderedList = sorted(glob.glob(os.path.join(pathToFiles, '*.fts')))
contador=0
for filename in orderedList:
    contador+=1
print("\nExisten {} archivos fts".format(contador))
