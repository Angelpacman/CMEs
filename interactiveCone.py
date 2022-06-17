from astropy.io import fits
import numpy as np
import datetime
import pandas as pd
import glob #para linux
from glob import os #para windows
import matplotlib.pyplot as plt
import subprocess
from astropy.visualization import (HistEqStretch, 
                                   simple_norm, 
                                   ImageNormalize, 
                                   imshow_norm, 
                                   MinMaxInterval, 
                                   LogStretch)
from sklearn import preprocessing
import matplotlib.animation as animation
from scipy import ndimage, misc, interpolate
from scipy.optimize import curve_fit

pathToFiles = '../DATA'
orderedList = sorted(glob.glob(os.path.join(pathToFiles, '*.fts')))
contador=0
for filename in orderedList:
    contador+=1
print("\nExisten {} archivos fts".format(contador))


# Obtener lista de datos que son 2048x2048
listaDatos = []
for filename in orderedList:
  hdulist = fits.open(filename)
  data = hdulist[0].data
  hdulist.close
  if data.shape == (2048, 2048):
    listaDatos.append(filename)
print("\nDe los cuales {} son  de tamaño 2048 x 2048".format(len(listaDatos)))

# Obtener headers
hdulist1 = fits.open(listaDatos[0])
headerList = hdulist1[0].header
hdulist1.close()
maph=pd.DataFrame(columns=headerList)

# Generar mapa de imagenes
cel = []
for i in range(len(listaDatos)):
    hdulist1 = fits.open(listaDatos[i])
    headerList = hdulist1[0].header
    hdulist1.close()
    for headerInplace in headerList:
        cel.append(headerList[headerInplace])
    maph.loc[i]=cel
    cel=[]

# Filtrar lista de datos
expTime=maph.EXPTIME
mExpTime=np.mean(expTime)
stdExpTime=np.std(maph.EXPTIME)
#maph = maph[(maph.EXPTIME >= (mExpTime - stdExpTime)) & (maph.EXPTIME <= (mExpTime + stdExpTime))]
maph = maph[(maph.EXPTIME >= (mExpTime - 1.5)) & (maph.EXPTIME <= (mExpTime + 1.5))]
tama = len(maph)
listaDatosFiltrados = []

for fit_index in maph.index:
    listaDatosFiltrados.append(listaDatos[fit_index])

mapa = []

for fit in listaDatosFiltrados:
    hdulist = fits.open(fit)
    datatemp = hdulist[0].data
    hdulist.close()
    mapa.append(datatemp)


#Suavizar mapa conf filtro
map_s = ndimage.uniform_filter(mapa, size=3, mode='reflect')

# Obtener constantes de los fits
mitad=1 #esto es para elegir el tamaño de la imagen 2 -> 512; 1 -> 1024
xc=maph.loc[0].CRPIX1/mitad
yc=maph.loc[0].CRPIX2/mitad
print("xc = {}".format(xc))
print("yc = {}".format(yc))

for i in maph:
    if i == 'R_SUN':
        print("R_SUN existe")
    elif i == 'RSUN':
        print("RSUN existe")

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

# Deal with times
df_tiempo=pd.to_datetime(maph['DATE-OBS'])
df_tiempo.head()

tiempos = [maph['DATE-OBS']]
header = ["DATE-OBS"]
df_tiempo = pd.concat(tiempos, axis = 1, keys = header)
df_tiempo.head()

tiempos=pd.to_datetime(df_tiempo['DATE-OBS'])
t0 = tiempos[0] # Esto fallará si se filtro el FITfile con el indice 0
delta = []
for celda in tiempos:
    delta.append((celda - t0).total_seconds())
print("Seconds on the array delta from zero to final: ")
delta

df_tiempo['SECONDS'] = delta
df_tiempo.head()
dx = maph.loc[0]['CDELT1']*1
print("dx = .{}".format(dx))

###########################################################################################
# Set Cone
n_puntos=100
ang_inc=3.0
print('centro = {}, {}'.format(xc,yc))

# radio inicial y final
if maph.loc[0].DETECTOR == 'COR2': #para stereo
    rrin=3
    rrfin=15
else:
    rrin=6
    rrfin=15
rrin

# funcion para dibujar x, y del cono (angulo principal, espacio direcciones, numero de direcciones)
def setAngle(angulo, ang_inc, gap):
    gap_angle = np.arange(gap*(-1),gap+1)
    rr=rrt=np.zeros((1+gap*2,2,n_puntos),dtype=int)
    PA_m = angulo + 90
    rsol = r0
    index = np.arange(100)
    rads = np.linspace(rrin,rrfin,100)
    radios = np.zeros(n_puntos)
    for j in gap_angle:
        for i in index:
            rd = rads[i]
            radios[i] = rd
            radio=rsol * rd
            teta = (PA_m + ang_inc * j)*np.pi/180
            x = radio * np.cos(teta) + xc
            y = radio * np.sin(teta) + yc
            if x < 0 or y < 0:
                print(i, x, y)
            rr[j+gap][0][i]=x
            rr[j+gap][1][i]=y
    return rr

rr = setAngle(80,3,5)
print(rr.shape)

# Librerias para iterar colores xd
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import matplotlib.cm as mplcm
import matplotlib.colors as colors
# Dar formato a las fechas en formato datatime
from matplotlib import dates as mpl_dates


# Constantes
rSolarcsec=32*60
disterr=56./1920./2. ##mitad de la minima escala en Rs para LASCO C3
if (maph['DETECTOR'][0] == 'COR2'):
    disterr = maph['CDELT1'][0]/rSolarcsec/2
radios = np.linspace(rrin,rrfin,100)
derr= disterr*radios**2

# Tomar valores para el cono dibujado
def fltarr(a,b,c):
    """
    FLTARR(a, b, c)    -> np.zeros((c, b, a))
    """        
    return np.zeros((c,b,a))

output_py=fltarr(n_puntos,rr.shape[0],tama) ## con el cambio en la dir radial
#output_idl=fltarr(n_puntos,rr.shape[0],tama) ## con el cambio en la dir radial

for i in range(tama):
    for j in range(rr.shape[0]):
        for k in range(100):
            output_py[i][j][k]  = map_s[     i  ][ rr[j,1,k] ][ rr[j,0,k] ]
            #output_idl[i][j][k] = s_idl[ str(i) ][ rr[j,1,k] ][ rr[j,0,k] ]

# Esta variable se usa para tomar manualmente indices para definir tiempo de background
indtbkg = np.arange(25,44)

# Salida del cono
def salidaCono():
    global output
    output=fltarr(n_puntos,rr.shape[0],tama)
    for i in range(tama):
        for j in range(rr.shape[0]):
            for k in range(100):
                output[i][j][k]  = map_s[     i  ][ rr[j,1,k] ][ rr[j,0,k] ]

# obtener el error en un arreglo
def backgroundError():
    global den_bg_er, den_bgs_er
    den_bg_er  = np.zeros((n_puntos,2))
    den_bgs_er = np.zeros((n_puntos,2))
    for i in range (n_puntos):
        den_bg_er[i][0] = np.mean(output.T[i])
        den_bg_er[i][1] = np.std( output.T[i])
        den_bgs_er[i][0] = np.mean(output[min(indtbkg):max(indtbkg)].T[i])
        den_bgs_er[i][1] = np.std( output[min(indtbkg):max(indtbkg)].T[i])

########################################################################################
# obtener error de manera manual usando output_py
den_bg_er_py = np.zeros((100,2))
for i in range(n_puntos):
    den_bg_er_py[i][0] = np.mean(output_py.T[i])
    den_bg_er_py[i][1] = np.std(output_py.T[i])

# 
den_bgs_er=np.zeros((n_puntos,2))  #la den y el sigma de la den del background
for i in range(n_puntos):
    #den_bgs_er[i][0] = np.mean(output_py.T[i][:][min(indtbkg):max(indtbkg)])
    #den_bgs_er[i][1] = np.std( output_py.T[i][:][min(indtbkg):max(indtbkg)])
    den_bgs_er[i][0] = np.mean(output_py[ min(indtbkg):max(indtbkg) ].T[ i ])
    den_bgs_er[i][1] = np.std( output_py[ min(indtbkg):max(indtbkg) ].T[ i ])
#den_bgs_er
#######################################################################################
