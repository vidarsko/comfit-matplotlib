import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import comfit_matplotlib as cfm
import comfit as cf

import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import scipy as sp

# 2D system
bs = cf.BaseSystem(2,xRes=31,yRes=31)
field = bs.x**2+bs.y**2
fig = bs.plot_field(field)
plt.show()