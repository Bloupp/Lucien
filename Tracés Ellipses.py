## Tracé ellipses

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse

e=Ellipse((0,0),5,2,0,linewidth=1, fill=False, zorder=1)

fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})
ax.add_artist(e)

ax.set_xlim(-10, 10)
ax.set_ylim(-10, 10)

plt.show()