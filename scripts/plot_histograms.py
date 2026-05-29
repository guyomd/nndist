import numpy as np
import pandas as pd
import pygmt


df = pd.read_csv('output_declust.txt', sep=";")
nn = df[' nn_distance'].values.astype(np.float64)  # Nearest-neighbor distances
if ' norm_prox' in df.columns:
    nd = df[' norm_prox'].values.astype(np.float64)  # Normalized proximities

# Make histograms:
xmin = -6.0
xmax = 6.0
ymin = 0.0
ymax = 30.0

fig = pygmt.Figure()
fig.plot(x=[0.0, 0.0], y=[ymin, ymax], pen='1p,black,dotted')
fig.histogram(data=np.log10(nn[1:]), 
              histtype=1, 
              frame=['WSne'],
              series=0.3, 
              fill='lightblue', 
              pen='1p', 
              region=[xmin, xmax, ymin, ymax],
              transparency=30, 
              label='Nearest-neighbor distance')
if ' norm_prox' in df.columns:
    fig.histogram(data=np.log10(nd[1:]), 
                  histtype=1, 
                  frame=['WSne', 'x+llog@-10@- distance', 'y+lFrequency (%)'],
                  series=0.3, 
                  fill='lightorange', 
                  pen='1p', 
              region=[xmin, xmax, ymin, ymax], 
              transparency=40, 
              label='Normalized proximity')
fig.legend()

#fig.show()
fig.savefig("histograms.png")