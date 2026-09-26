import matplotlib.pyplot as plt
from itertools import chain
from pathlib import Path

x_mod_secant_e1 = [-0.29,-0.24,-0.19,-0.13,-0.08,-0.02,0.02,0.07,0.13,0.18,0.24,0.29,0.35,0.4,0.45,0.51,0.56,0.62,0.67,0.72,0.78,0.83,0.89,0.94,1]
y_average_mod_secant_e1 = [7.55555555555556,8.44444444444445,8.44444444444445,8.66666666666667,8.44444444444445,8.66666666666667,9.44444444444445,8.88888888888889,8.88888888888889,8.33333333333333,8.11111111111111,8.11111111111111,8.11111111111111,7.5,8.22222222222222,8.41666666666667,8.63888888888889,8.44444444444445,8.55555555555556,8.52777777777778,8.55555555555556,8.5,7.55555555555556,7.55555555555556,11]
y_min_mod_secant_e1 = [2,7,8,8,8,8,9,8,8,8,8,8,8,7,7,7,7,7,6,6,6,6,6,7,11]
y_max_mod_secant_e1 = [100,14,10,14,10,10,10,11,12,11,10,10,10,9,10,11,11,11,11,11,11,10,8,11,11]

x_mod_secant_e2 = [-0.29,-0.24,-0.19,-0.13,-0.08,-0.02,0.02,0.07,0.13,0.18,0.24,0.29,0.35,0.4,0.45,0.51,0.56,0.62,0.67,0.72,0.78,0.83,0.89,0.94,0.99]
y_average_mod_secant_e2 = [15.3333333333333,8.44444444444445,8.44444444444445,8.66666666666667,8.44444444444445,8.66666666666667,9.44444444444445,8.88888888888889,8.88888888888889,8.33333333333333,8.11111111111111,8.11111111111111,8.11111111111111,7.5,8.22222222222222,8.41666666666667,8.63888888888889,8.44444444444445,8.55555555555556,8.52777777777778,8.55555555555556,8.5,7.55555555555556,7.55555555555556,12.3142857142857]
y_min_mod_secant_e2 = [4,7,8,8,8,8,9,8,8,8,8,8,8,7,7,7,7,7,6,6,6,6,6,7,11]
y_max_mod_secant_e2 = [100,14,10,14,10,10,10,11,12,11,10,10,10,9,10,11,11,11,11,11,11,10,8,11,15]

x_brentq_e2 = [-0.29,-0.24,-0.19,-0.13,-0.08,-0.02,0.02,0.07,0.13,0.18,0.24,0.29,0.35,0.4,0.45,0.51,0.56,0.62,0.67,0.72,0.78,0.83,0.89,0.94,0.99,1]
y_average_brentq_e2 = [13.5,9.44444444444445,8.44444444444445,8,7.77777777777778,8,8.55555555555556,8.44444444444445,7.88888888888889,8.77777777777778,8.88888888888889,9.11111111111111,9.44444444444445,9.88888888888889,10.1111111111111,9.88888888888889,10.1111111111111,10,9.88888888888889,10.3333333333333,10.8888888888889,11.3333333333333,11.6666666666667,10.1111111111111,21.4571428571429,21]
y_min_brentq_e2 = [3,9,8,7,7,8,8,7,7,8,8,9,9,9,9,9,9,9,8,9,8,8,9,10,21,21]
y_max_brentq_e2 = [15,10,9,9,8,8,9,10,9,9,10,10,10,10,11,11,11,11,11,11,12,12,12,11,25,21]

x_bisect_e2 = [-0.29,-0.24,-0.19,-0.13,-0.08,-0.02,0.02,0.07,0.13,0.18,0.24,0.29,0.35,0.4,0.41,0.45,0.51,0.56,0.62,0.67,0.72,0.78,0.83,0.89,0.94,0.99]
y_average_bisect_e2 = [20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20]
y_min_bisect_e2 = [20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20]
y_max_bisect_e2 = [20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20]


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
ax.set_ylim(min_y, 25)

ax.grid(linewidth=0.5)

ax.set_title("Square with Circular Opening Section - Iteration Count and Axial Load Ratio")
ax.set_xlabel("Pn-target / Po")
ax.xaxis.set_label_coords(0.3, -0.1)
ax.set_ylabel("Iteration to convergence")

ax.legend(loc='upper right', ncols=1, bbox_to_anchor=(1, -0.08), framealpha=0.4)
plt.tight_layout()

script_dir = Path(__file__).resolve().parent
save_path = script_dir / "iteration-square-with-circular-opening.svg"
plt.savefig(save_path, bbox_inches="tight", pad_inches=0)

plt.show()
