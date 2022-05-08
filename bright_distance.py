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
