import matplotlib.pyplot as plt
from itertools import chain
from pathlib import Path

x_mod_secant_e1 = [-0.307,-0.253,-0.198,-0.144,-0.089,-0.035,0.019,0.073,0.128,0.182,0.237,0.291,0.346,0.4,0.455,0.509,0.564,0.618,0.673,0.727,0.782,0.836,0.891,0.945,1]
y_average_mod_secant_e1 = [2.22222222222222,8.44444444444445,7,7.44444444444445,7.22222222222222,8.11111111111111,9.55555555555556,8.55555555555556,8.44444444444445,8,8,8,7.77777777777778,6.66666666666667,7.22222222222222,7.66666666666667,8.33333333333333,8.66666666666667,8.88888888888889,8.88888888888889,8.22222222222222,8.11111111111111,6.77777777777778,6.77777777777778,1]
y_min_mod_secant_e1 = [2,7,6,6,6,7,8,8,8,8,8,8,7,6,7,7,7,6,6,6,6,5,5,6,1]
y_max_mod_secant_e1 = [4,14,9,9,9,9,10,9,9,8,8,8,8,8,8,9,9,9,10,10,9,9,7,7,1]

x_mod_secant_e2 = [-0.307,-0.253,-0.198,-0.144,-0.089,-0.035,0.019,0.073,0.128,0.182,0.237,0.291,0.346,0.4,0.455,0.509,0.564,0.618,0.673,0.727,0.782,0.836,0.891,0.945,1]
y_average_mod_secant_e2 = [9.55555555555556,8.94444444444445,7.55555555555556,7.88888888888889,7.66666666666667,8.55555555555556,9.77777777777778,8.77777777777778,8.66666666666667,8,8.22222222222222,8,8,7.11111111111111,8,8.44444444444445,8.33333333333333,8.77777777777778,9.22222222222222,8.88888888888889,8.88888888888889,8.44444444444445,7.33333333333333,7.11111111111111,1]
y_min_mod_secant_e2 = [4,8,7,7,7,8,8,8,8,8,8,8,8,7,8,8,7,7,7,6,6,6,6,7,1]
y_max_mod_secant_e2 = [11,15,10,9,9,9,11,9,9,8,9,8,8,8,8,9,9,9,10,10,10,9,8,8,1]

x_brentq_e2 = [-0.3,-0.25,-0.19,-0.14,-0.08,-0.03,0.01,0.07,0.12,0.18,0.23,0.29,0.34,0.4,0.45,0.5,0.56,0.61,0.67,0.72,0.78,0.83,0.89,0.94,1]
y_average_brentq_e2 = [10.7777777777778,8.83333333333333,8.11111111111111,7.66666666666667,7.11111111111111,7.44444444444445,7.88888888888889,8,7.77777777777778,8,8.22222222222222,8.55555555555556,9,9,9.44444444444445,9.44444444444445,9.66666666666667,10.1111111111111,10,10,10.4444444444444,10.6666666666667,10.6666666666667,9,1]
y_min_brentq_e2 = [3,8,7,7,6,7,7,7,7,8,8,8,9,9,9,9,9,9,8,8,8,8,8,9,1]
y_max_brentq_e2 = [12,11,9,8,8,8,9,9,8,8,9,9,9,9,10,10,10,11,11,11,11,11,11,9,1]

x_bisect_e2 = [-0.3,-0.25,-0.19,-0.14,-0.08,-0.03,0.01,0.07,0.12,0.18,0.23,0.29,0.34,0.4,0.45,0.5,0.56,0.61,0.67,0.72,0.78,0.83,0.89,0.94,1]
y_average_bisect_e2 = [14.2941176470588,19.1666666666667,19.3333333333333,19.3888888888889,19,19.1111111111111,19.2777777777778,19.1111111111111,19.3888888888889,19.2777777777778,19.2222222222222,19.1666666666667,19.0555555555556,19.2105263157895,19.2222222222222,19.1666666666667,19.0555555555556,19.3333333333333,19.1666666666667,19.0555555555556,18.7777777777778,19,18.7777777777778,17.2222222222222,1]
y_min_bisect_e2 = [13,17,18,19,17,17,18,14,18,17,16,17,17,16,16,17,17,18,17,17,16,17,17,14,1]
y_max_bisect_e2 = [19,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,19,1]


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

ax.set_title("Square Section - Iteration Count and Axial Load Ratio")
ax.set_xlabel("Pn-target / Po")
ax.xaxis.set_label_coords(0.3, -0.1)
ax.set_ylabel("Iteration to convergence")

ax.legend(loc='upper right', ncols=1, bbox_to_anchor=(1, -0.08), framealpha=0.4)
plt.tight_layout()

script_dir = Path(__file__).resolve().parent
save_path = script_dir / "iteration-square.svg"
plt.savefig(save_path, bbox_inches="tight", pad_inches=0)

plt.show()
