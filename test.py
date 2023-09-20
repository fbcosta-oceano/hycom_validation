import datetime
import remo
import time
import numpy as np
import scipy.io
import os
import sys
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
import matplotlib.dates as mdates
import cartopy.mpl.ticker as cticker
import cartopy.feature as cfeature
import warnings

warnings.filterwarnings("ignore")

cfg_file = sys.argv[1]
opt = open(cfg_file, "r")
content_list = opt.readlines()
expt = list(content_list[0].rstrip("\n").split(","))
rodada = list(content_list[1].rstrip("\n").split(","))
legenda = list(content_list[2].rstrip("\n").split(","))
data_inicial = str(content_list[3].rstrip("\n"))
data_final = str(content_list[4].rstrip("\n"))
inc_tempo = float(content_list[5].rstrip("\n"))

INPUT_DIR = "/home/filipe.costa/resultados"
fig_output_dir = "/home/filipe.costa/resultados/figs"

data_inicial = datetime.datetime(*time.strptime(data_inicial, "%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final, "%d/%m/%Y")[0:3])







under = "#000033"
over = "#330000"
bad = "#FFFFFF"
cmap = matplotlib.cm.get_cmap("jet")
