import matplotlib.pyplot as plt
import numpy as np
from itertools import chain
from scipy.interpolate import PchipInterpolator
import matplotlib.ticker as mtick
from pathlib import Path

x_brentq = [0.00001,0.0001,0.001,0.01,0.1,1]
y_average_brentq = [0.000000008618732327292,0.0000000637209367242553,0.000000881455431390134,0.0000160408051608856,0.000189779954948825,0.0020986584340042]
x_mod_secant = [0.00001,0.0001,0.001,0.01,0.1,1]
y_average_mod_secant = [0.00000000282806588522127,0.0000000430687213592053,0.000000208226108101891,0.00000128597607835323,0.0000329581554488887,0.000758061976014223]


x_smooth = np.linspace(min(x_mod_secant), max(x_mod_secant), 300)
pchip_mod_secant = PchipInterpolator(x_mod_secant, y_average_mod_secant)
y_smooth_mod_secant = pchip_mod_secant(x_smooth)

pchip_brentq = PchipInterpolator(x_brentq, y_average_brentq)
y_smooth_brentq = pchip_brentq(x_smooth)

fig, ax = plt.subplots()

ax.plot(x_smooth, y_smooth_mod_secant, 'r', label='This study')
ax.plot(x_smooth, y_smooth_brentq, 'b', label='Brent Method')

min_y = min(chain(y_average_mod_secant, y_average_brentq))
max_y = max(chain(y_average_mod_secant, y_average_brentq))
ax.set_xlim(min(x_mod_secant),max(x_mod_secant))
ax.set_ylim(min_y, max_y*1.1)

ax.grid(linewidth=0.5)

ax.set_title("Average Error and Tolerance (L-shape section)")
ax.set_xlabel("Tolerance")
ax.set_ylabel("Average error")

ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1))

ax.set_xscale('log')

ax.legend(loc='upper right', ncols=2, bbox_to_anchor=(0.75, -0.16), framealpha=0.4)
plt.tight_layout()

script_dir = Path(__file__).resolve().parent
save_path = script_dir / "error-vs-tolerance-L.svg"
plt.savefig(save_path, bbox_inches="tight", pad_inches=0)


plt.show()
