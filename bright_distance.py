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
from scipy import ndimage, misc

pathToFiles = 'DATA'
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


# Manejo de los tiempos
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

df_tiempo['SECONDS'] = delta
df_tiempo.head()

# Definiendo el cono
dx = maph.loc[0]['CDELT1']*1
n_puntos=100
ang_inc=3.0
print('centro = {}, {}'.format(xc,yc))

# Distancia de las direcciones del con en radios Solares
if maph.loc[0].DETECTOR == 'COR2': #para stereo
    rrin=3
    rrfin=15
else:
    rrin=6
    rrfin=15
rrin

# Funcion para generar rr (dibujo de las direcciones del cono)
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

# bosquejo
rr = setAngle(80,3,5)
print(rr.shape)
for i in np.arange(len(rr)):
    plt.plot(rr[i,0,:],rr[i,1,:])
#plt.xlim([0,2040])
#plt.ylim([750,1600])
plt.scatter(xc, yc, s=2*maph.RSUN.iloc[30], facecolors='none', edgecolors='b')
plt.grid()
plt.show()

#librerias para iterar colores xd
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import matplotlib.cm as mplcm
import matplotlib.colors as colors
#dar formato a las fechas en formato datatime
from matplotlib import dates as mpl_dates

# radio central PA analizado en el tiempo
fig = plt.figure(figsize = (18,6))
ax = fig.add_subplot(1, 1, 1)
dic_radios={}

NUM_COLORS = 100
cm = plt.get_cmap('gnuplot2_r')
cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)

for radio in range(0,n_puntos,5):
    valorRadioCentral = np.array([])
    for indice in range(len(listaDatosFiltrados)):
        map = map_s[indice]
        valorRadioCentral = np.append(valorRadioCentral, [ map[ rr[len(rr)//2, 1, radio], rr[len(rr)//2, 0, radio] ] ])
    brilloProm = np.mean(valorRadioCentral)
    valorRadioCentral = valorRadioCentral - brilloProm
    dic_radios["radio_No.%s" %radio]=valorRadioCentral

    #iteración de colores
    ax.set_prop_cycle(color=[scalarMap.to_rgba(radio)])
    ax.plot_date(tiempos, dic_radios["radio_No.%s" %radio], linestyle='solid',label='Punto %s' %radio)
    #plt.xticks(rotation=45,  ha='right')

date_format = mpl_dates.DateFormatter('%H:%M \n%d-%b-%Y')
plt.gca().xaxis.set_major_formatter(date_format)
ax.grid(color = "#242326")
ax.legend(ncol=2)
ax.set_facecolor('black')
plt.ylabel("Brillo - $\mu_b$")
plt.title("$Brillo - \mu_{brillo}$ vs Tiempo (en la dirección central PA del cono)")
plt.show()
