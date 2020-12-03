import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import pandas as pd

x = []
y = []
z = []
res = ["GLY", "PRO", "ASP", "GLU", "LYS", "ARG", "HIS", "SER", "THR", "ASN", "GLN", "ALA", "MET", "TYR", "TRP", "VAL", "ILE", "LEU", "PHE", "CYS"]
frame = []
for k in range(20):
    for i in range(20):
        for j in range(20):
            x.append(i+1)
            y.append(j+1)
            z.append(k)
            frame.append(res[k])
heat = np.random.randint(100, size=(20, 20, 20)).flatten()
df = pd.DataFrame({"x": x, "y": y, "z": z, "frame": frame, "heat": heat})

fig = px.scatter_3d(df, x="x", y="y", z="z", animation_frame="frame", color="heat", range_x=[0,21], range_y=[0,21], range_z=[0, 21])
fig.update_layout(scene= dict(
    xaxis = dict(
        tickmode = 'array',
        ticktext = res,
        tickvals = list(range(1, 21))
    ),
yaxis = dict(
        tickmode = 'array',
        ticktext = res,
        tickvals = list(range(1, 21))
    ),
zaxis = dict(
        tickmode = 'array',
        ticktext = res,
        tickvals = list(range(0, 20))
    )

    )
)
fig["layout"].pop("updatemenus")

fig.show()