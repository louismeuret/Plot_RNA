
import barnaba as bb
import pandas as pd
import matplotlib.pyplot as plt

traj_xtc = "/home/cya/ratchet_26_rna/md_noPBC.xtc"
top_pdb = "/home/cya/ratchet_26_rna/init_conf.gro" 
native = "/home/cya/ratchet_26_rna/native.pdb"

stackings, pairings, res = bb.annotate(traj_xtc,top_pdb)
df = pd.DataFrame(stackings)
print(df)
"""
angles,res = bb.backbone_angles(traj_xtc,native_pdb)
df = pd.DataFrame(angles)
"""
ermsd = bb.ermsd(native,traj_xtc,topology=top_pdb)
df = pd.DataFrame(ermsd)
#plt.plot(df)
#plt.show()


rvecs, res = bb.dump_rvec(traj_xtc, top_pdb)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import seaborn as sns

import numpy as np
import plotly.graph_objects as go

# Calculate nonzero elements
nonzero = np.where(np.sum(rvecs**2, axis=3) > 0.01)
rr = rvecs[nonzero]

# Calculate rho and zeta
z = rr[:, 2]
rho = np.sqrt(rr[:, 0]**2 + rr[:, 1]**2)

# Create scatter plot
fig = go.Figure()

fig.add_trace(go.Scatter(
    x=rho,
    y=z,
    mode='markers',
    marker=dict(
        size=5,
        color='blue',
        opacity=0.5
    )
))

# Add ellipses
f1 = 0.5
f2 = 0.3

fig.add_shape(
    type="circle",
    xref="x",
    yref="y",
    x0=-f1,
    y0=-f2,
    x1=f1,
    y1=f2,
    line=dict(
        color="red",
        width=2,
        dash="dashdot",
    ),
)

fig.add_shape(
    type="circle",
    xref="x",
    yref="y",
    x0=-2*f1,
    y0=-2*f2,
    x1=2*f1,
    y1=2*f2,
    line=dict(
        color="orange",
        width=2,
        dash="dashdot",
    ),
)

fig.add_shape(
    type="circle",
    xref="x",
    yref="y",
    x0=-3*f1,
    y0=-3*f2,
    x1=3*f1,
    y1=3*f2,
    line=dict(
        color="yellow",
        width=2,
        dash="dashdot",
    ),
)

# Update layout
fig.update_layout(
    title="Scatter Plot with Ellipses",
    xaxis_title=r'$\rho$ (nm)',
    yaxis_title='z (nm)',
    yaxis=dict(
        tickvals=[-1.2, -0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9, 1.2]
    ),
    xaxis=dict(
        tickvals=[0, 0.5, 1.0, 1.5]
    )
)

fig.show()



