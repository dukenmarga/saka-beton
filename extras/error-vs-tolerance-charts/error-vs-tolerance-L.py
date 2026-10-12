import matplotlib.pyplot as plt
import numpy as np
from itertools import chain
from scipy.interpolate import PchipInterpolator
import matplotlib.ticker as mtick
from pathlib import Path

x_brentq = [0.00001,0.0001,0.001,0.01,0.1,1]
y_average_brentq = [0.00000000909219598176471,0.0000000622084178842779,0.000000655026771379396,0.00000966971890658108,0.000149557860021294,0.00148590577291983]
x_mod_secant = [0.00001,0.0001,0.001,0.01,0.1,1]
y_average_mod_secant = [0.00000000194728782860987,0.0000000386015528925806,0.000000299295289787715,0.00000193678951370792,0.0000787409104755844,0.000416659678208043]


x_smooth = np.logspace(np.log10(min(x_mod_secant)), np.log10(max(x_mod_secant)), 300)
pchip_mod_secant = PchipInterpolator(x_mod_secant, y_average_mod_secant)
y_smooth_mod_secant = pchip_mod_secant(x_smooth)

pchip_brentq = PchipInterpolator(x_brentq, y_average_brentq)
y_smooth_brentq = pchip_brentq(x_smooth)

fig, ax = plt.subplots()

ax.plot(x_smooth, y_smooth_mod_secant, 'r', label='This study')
ax.plot(x_smooth, y_smooth_brentq, 'b', label='Brent Method')
ax.plot(x_mod_secant, y_average_mod_secant, 'ro', ms=4)
ax.plot(x_brentq, y_average_brentq, 'bo', ms=4)

min_y = min(chain(y_average_mod_secant, y_average_brentq))
max_y = max(chain(y_average_mod_secant, y_average_brentq))
ax.set_xlim(min(x_mod_secant),max(x_mod_secant))
ax.set_ylim(min_y, max_y*1.1)

ax.grid(linewidth=0.5)

ax.set_title("Average Error and Tolerance (L-shape section)")
ax.set_xlabel("Tolerance (kN)")
ax.set_ylabel("Average error")

ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1))

ax.set_xscale('log')

ax.legend(loc='upper right', ncols=2, bbox_to_anchor=(0.75, -0.12), framealpha=0.4)
plt.tight_layout()

script_dir = Path(__file__).resolve().parent
save_path = script_dir / "error-vs-tolerance-L.svg"
plt.savefig(save_path, bbox_inches="tight", pad_inches=0)


plt.show()
