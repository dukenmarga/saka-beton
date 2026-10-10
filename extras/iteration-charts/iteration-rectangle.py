import matplotlib.pyplot as plt
from itertools import chain
from pathlib import Path

x_mod_secant_e1 = [-0.44,-0.38,-0.32,-0.26,-0.2,-0.14,-0.08,-0.02,0.03,0.09,0.15,0.21,0.27,0.33,0.39,0.45,0.51,0.57,0.63,0.69,0.75,0.81,0.87,0.93,1]
y_average_mod_secant_e1 = [2.11111111111111,7.44444444444445,7.61111111111111,8.44444444444445,8.38888888888889,9.61111111111111,9.11111111111111,9,8.83333333333333,8.22222222222222,7.94444444444445,8.11111111111111,7.77777777777778,8,7.11111111111111,8.11111111111111,8.27777777777778,8.77777777777778,9,8.83333333333333,8.33333333333333,7.94444444444445,7.22222222222222,6.61111111111111,1]
y_min_mod_secant_e1 = [2,6,6,7,8,8,8,8,8,8,7,8,7,7,6,8,7,7,7,6,6,5,5,6,1]
y_max_mod_secant_e1 = [3,10,10,13,10,11,10,10,12,9,9,9,8,9,8,9,9,9,10,10,10,9,8,7,1]

x_mod_secant_e2 = [-0.44,-0.38,-0.32,-0.26,-0.2,-0.14,-0.08,-0.02,0.03,0.09,0.15,0.21,0.27,0.33,0.39,0.45,0.51,0.57,0.63,0.69,0.75,0.81,0.87,0.93,1]
y_average_mod_secant_e2 = [12.7777777777778,7.83333333333333,7.88888888888889,9.22222222222222,9.22222222222222,10.3333333333333,9.72222222222222,9.22222222222222,9.16666666666667,8.61111111111111,8.44444444444445,8.44444444444445,8.27777777777778,8,7.5,8.27777777777778,8.83333333333333,9.16666666666667,9.22222222222222,8.88888888888889,8.83333333333333,8.33333333333333,7.55555555555556,6.94444444444445,1]
y_min_mod_secant_e2 = [4,6,6,8,8,8,9,8,8,8,8,8,8,7,7,8,7,7,7,7,6,6,6,6,1]
y_max_mod_secant_e2 = [100,11,10,14,10,12,10,10,12,10,9,9,9,9,9,9,9,10,10,10,10,9,8,8,1]

x_brentq_e2 = [-0.44,-0.38,-0.32,-0.26,-0.2,-0.14,-0.08,-0.02,0.03,0.09,0.15,0.21,0.27,0.33,0.39,0.45,0.51,0.57,0.63,0.69,0.75,0.81,0.87,0.93,1]
y_average_brentq_e2 = [11.4722222222222,9.77777777777778,9,7.83333333333333,7,8.05555555555556,8.05555555555556,7.66666666666667,7.77777777777778,8.05555555555556,8.11111111111111,8.38888888888889,9.22222222222222,9.55555555555556,9.88888888888889,9.94444444444445,10.1666666666667,10.2222222222222,10.4444444444444,10.2222222222222,10.2777777777778,10.6111111111111,10.7222222222222,10.0555555555556,1]
y_min_brentq_e2 = [3,6,7,7,6,7,6,7,7,7,8,7,9,9,9,9,9,10,10,9,8,8,7,8,1]
y_max_brentq_e2 = [13,12,10,9,8,10,9,9,9,9,9,9,10,10,10,10,11,11,11,11,11,11,12,12,1]

x_bisect_e2 = [-0.44,-0.38,-0.32,-0.26,-0.2,-0.14,-0.08,-0.02,0.03,0.09,0.13,0.15,0.21,0.27,0.33,0.39,0.45,0.51,0.57,0.63,0.69,0.75,0.81,0.87,0.93,1]
y_average_bisect_e2 = [14.4571428571429,19.3333333333333,19.4444444444444,19.3888888888889,19.4444444444444,19.4444444444444,19.4444444444444,19.4444444444444,19.3333333333333,19.2222222222222,19,19,19.1111111111111,19.4444444444444,19.4444444444444,19.4444444444444,18.5555555555556,19.2222222222222,19.3333333333333,19.1111111111111,19.3888888888889,19.4444444444444,19,19.3333333333333,18.7777777777778,1]
y_min_bisect_e2 = [13,19,19,18,19,19,19,19,19,18,19,15,17,19,19,19,15,18,17,18,18,19,15,18,18,1]
y_max_bisect_e2 = [18,20,20,20,20,20,20,20,20,20,19,20,20,20,20,20,20,20,20,20,20,20,20,20,20,1]


fig, ax = plt.subplots()

ax.plot(x_mod_secant_e1, y_average_mod_secant_e1, color='navy', marker="*", linestyle='dashed', label='This study (avg) tol=1e-1')
ax.plot(x_mod_secant_e2, y_average_mod_secant_e2, color='red', marker=None, linestyle='dashed', label='This study (avg) tol=1e-2')
ax.plot(x_brentq_e2, y_average_brentq_e2, color='blue', marker="o", linestyle='dashed', label='Brent Method (avg) tol=1e-2')
ax.plot(x_bisect_e2, y_average_bisect_e2, color='green', marker=".", linestyle='dashed', label='Bisection Method (avg) tol=1e-2')

# ax.fill_between(x_mod_secant_e1, y_min_mod_secant_e1, y_max_mod_secant_e1,color="#9498B8", edgecolor="#6B76F2", alpha=0.5, label='This study (min/max) tol=1e-1')
# ax.fill_between(x_mod_secant_e2, y_min_mod_secant_e2, y_max_mod_secant_e2,color="#FACBCB", edgecolor="#FDADAD", alpha=0.5, label='This study (min/max) tol=1e-2')
# ax.fill_between(x_brentq_e2, y_min_brentq_e2, y_max_brentq_e2,color="#CBDBFA", edgecolor="#9AB9F8", alpha=0.5, label='Brent Method (min/max) tol=1e-2')
# ax.fill_between(x_bisect_e2, y_min_bisect_e2, y_max_bisect_e2,color="#CFFACB", edgecolor="#9AF8A5", alpha=0.5, label='Bisection Method (min/max) tol=1e-2')

min_y = min(chain(y_min_mod_secant_e1, y_min_mod_secant_e2, y_min_brentq_e2, y_min_bisect_e2))
max_y = max(chain(y_max_mod_secant_e1, y_max_mod_secant_e2, y_max_brentq_e2, y_max_bisect_e2))
ax.set_xlim(min(x_mod_secant_e2),max(x_mod_secant_e2))
ax.set_ylim(min_y, 23)

ax.grid(linewidth=0.5)

ax.set_title("Rectangle Section - Iteration Count and Axial Load Ratio")
ax.set_xlabel("Pn-target / Po")
ax.xaxis.set_label_coords(0.3, -0.1)
ax.set_ylabel("Iteration to convergence")

ax.legend(loc='upper right', ncols=1, bbox_to_anchor=(1, -0.08), framealpha=0.4)
plt.tight_layout()

script_dir = Path(__file__).resolve().parent
save_path = script_dir / "iteration-rectangle.svg"
plt.savefig(save_path, bbox_inches="tight", pad_inches=0)

plt.show()
