import matplotlib.pyplot as plt
import numpy as np

def errorfill(x, y, yerr, color=None, alpha_fill=0.3, ax=None, label = None):
    ax = ax if ax is not None else plt.gca()
    if color is None:
        color = 'k'
#        color = ax._get_lines.color_cycle.next()
    if np.isscalar(yerr) or len(yerr) == len(y):
        ymin = y - yerr
        ymax = y + yerr
    elif len(yerr) == 2:
        ymin, ymax = yerr
    ax.plot(x, y, color=color, label = label)
    ax.fill_between(x, ymax, ymin, color=color, alpha=alpha_fill)