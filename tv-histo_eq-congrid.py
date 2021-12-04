#############################################################

from astropy.io import fits
import numpy as np
import datetime
import pandas as pd

import glob  # para linux
from glob import os  # para windows

import matplotlib.pyplot as plt
import subprocess

from time import sleep
import sys

#############################################################

# Ordenar (por nombre) los archivos en una lista para posterior consulta
pathToFiles = "../DATA"
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

from astropy.visualization import (HistEqStretch, simple_norm, ImageNormalize, imshow_norm, MinMaxInterval, LogStretch)
from sklearn import preprocessing
import matplotlib.animation as animation

# Funcion para hacer video e imágenes con la lista de archivos como entrada
print("\nLeyendo lista de archivos...")
print("\nGenerando imagenes .png...")


def diff(files):
    n = len(files)
    frames = []
    if n > 0:
        for i in range(n - 1):  # Para desplegar todas las imagenes hacer range(n-1)
            # Abrir 2 FITS para restarlos
            hdulist1 = fits.open(files[i])
            datos1 = hdulist1[0].data
            hdulist1.close()

            hdulist2 = fits.open(files[i + 1])
            datos2 = hdulist2[0].data
            date = hdulist2[0].header['DATE-OBS']
            hdulist2.close()

            imag_diff = datos2 - datos1

            # Normalizar la diferencia
            normalized_data = preprocessing.normalize(imag_diff)

            # Configurar figura con mejora de hist_eq
            fig = plt.figure(figsize=(19.2, 14.4))
            ax = fig.add_subplot(1, 1, 1)
            im, norm = imshow_norm(imag_diff, ax, origin='lower',
                                   cmap='CMRmap_r',
                                   interval=MinMaxInterval(),
                                   stretch=HistEqStretch(normalized_data))

            # Obtener figrua y archivos .png
            fig.suptitle(date, fontsize=24)
            fig.colorbar(im)
            plt.close(fig)
            fig.savefig(pathToFiles + "/file%02d.png" % i)

            # Mostrar porcentaje de avance en consola
            j = (i + 1) / (n - 1)
            sys.stdout.write('\r')
            sys.stdout.write("[%-20s] %d%%" % ('=' * int(20 * j), 100 * j))
            sys.stdout.flush()
            sleep(0.25)

    # Obtener video junto a las imágenes .FIT y .png
    print("\n\nArchivos .png generados, Renderizando video...\n")
    os.chdir(pathToFiles)
    subprocess.call([
        'ffmpeg', '-framerate', '4', '-i', 'file%02d.png', '-r', '30', '-pix_fmt', 'yuv420p',
        'video_name.mp4'
    ])
    os.chdir("../")

# return fig
datos_diff = diff(listaDatos)
print("\nPrograma Completado")

#############################################################
#############################################################
