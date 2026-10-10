import matplotlib.pyplot as plt
from itertools import chain
from pathlib import Path

x_mod_secant_e1 = [-0.48,-0.42,-0.36,-0.29,-0.23,-0.17,-0.11,-0.05,0.01,0.07,0.13,0.19,0.25,0.31,0.38,0.44,0.5,0.56,0.62,0.69,0.75,0.81,0.87,0.93,1]
y_average_mod_secant_e1 = [2.11111111111111,7.25,7.66666666666667,8.33333333333333,8,8.11111111111111,8.38888888888889,8.27777777777778,8.16666666666667,8,7.66666666666667,7.80555555555556,7.80555555555556,7.66666666666667,7.94444444444445,7.69444444444445,8.61111111111111,8.91666666666667,8.75,8.22222222222222,7.66666666666667,6.52777777777778,5.88888888888889,7.05555555555556,1]
y_min_mod_secant_e1 = [2,6,6,7,7,7,7,7,7,7,6,6,7,7,6,6,7,6,6,5,5,4,4,6,1]
y_max_mod_secant_e1 = [4,13,12,10,10,10,10,10,9,9,9,9,9,9,10,9,11,10,10,9,8,8,7,8,1]

x_mod_secant_e2 = [-0.48,-0.42,-0.36,-0.29,-0.23,-0.17,-0.11,-0.05,0.01,0.07,0.13,0.19,0.25,0.31,0.38,0.44,0.5,0.56,0.62,0.69,0.75,0.81,0.87,0.93,1]
y_average_mod_secant_e2 = [8.5,7.55555555555556,8.19444444444445,8.61111111111111,8.77777777777778,8.75,8.88888888888889,8.72222222222222,8.63888888888889,8.22222222222222,7.97222222222222,8.16666666666667,8,8.13888888888889,8.30555555555556,8.22222222222222,8.88888888888889,9.27777777777778,9.19444444444445,8.69444444444445,8,7.16666666666667,6.61111111111111,7.63888888888889,1]
y_min_mod_secant_e2 = [3,7,7,7,8,7,8,7,7,7,6,6,7,7,7,7,7,6,6,5,6,5,5,6,1]
y_max_mod_secant_e2 = [10,14,13,11,10,11,10,10,10,10,9,9,9,10,10,9,12,11,10,10,9,8,8,9,1]

x_brentq_e2 = [-0.48,-0.42,-0.36,-0.29,-0.23,-0.17,-0.11,-0.05,0.01,0.07,0.13,0.19,0.25,0.31,0.38,0.44,0.5,0.56,0.62,0.69,0.75,0.81,0.87,0.93,1]
y_average_brentq_e2 = [10.3333333333333,8.5,8.41666666666667,7,6.63888888888889,7.61111111111111,8.16666666666667,7.94444444444445,7.86111111111111,7.88888888888889,8.02777777777778,8.63888888888889,9.22222222222222,9.22222222222222,9.5,9.66666666666667,9.77777777777778,9.86111111111111,9.97222222222222,9.80555555555556,9.72222222222222,9.58333333333333,9.47222222222222,9.77777777777778,1]
y_min_brentq_e2 = [3,5,7,6,5,7,7,7,7,7,7,7,8,7,7,8,8,7,8,8,8,8,8,8,1]
y_max_brentq_e2 = [12,11,11,8,8,9,11,9,9,9,9,11,11,11,11,11,11,11,11,10,11,12,11,11,1]

x_bisect_e2 = [-0.48,-0.42,-0.36,-0.29,-0.23,-0.17,-0.11,-0.05,0.01,0.07,0.13,0.19,0.25,0.31,0.38,0.44,0.5,0.56,0.62,0.69,0.75,0.81,0.87,0.93,1]
y_average_bisect_e2 = [13.5,17.7222222222222,17.9444444444444,17.4444444444444,17.5,18,18,17.8888888888889,17.7222222222222,17.7777777777778,17.6666666666667,17.8888888888889,17.9444444444444,17.6111111111111,17.5555555555556,18,17.8888888888889,17.6388888888889,17.9444444444444,17.3333333333333,16.9166666666667,17.5,17.0277777777778,15.3055555555556,1]
y_min_bisect_e2 = [12,16,17,11,13,18,18,16,15,14,14,14,17,15,14,18,16,16,17,13,10,16,14,11,1]
y_max_bisect_e2 = [18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,17,1]


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
ax.set_ylim(min_y, max_y+1)

ax.grid(linewidth=0.5)

ax.set_title("T-shape Section - Iteration Count and Axial Load Ratio")
ax.set_xlabel("Pn-target / Po")
ax.xaxis.set_label_coords(0.3, -0.1)
ax.set_ylabel("Iteration to convergence")

ax.legend(loc='upper right', ncols=1, bbox_to_anchor=(1, -0.08), framealpha=0.4)
plt.tight_layout()

script_dir = Path(__file__).resolve().parent
save_path = script_dir / "iteration-T.svg"
plt.savefig(save_path, bbox_inches="tight", pad_inches=0)

plt.show()
