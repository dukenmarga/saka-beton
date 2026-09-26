import matplotlib.pyplot as plt
from itertools import chain
from pathlib import Path

x_mod_secant_e1 = [-0.52,-0.46,-0.39,-0.33,-0.27,-0.2,-0.14,-0.07,-0.01,0.04,0.11,0.17,0.23,0.3,0.36,0.42,0.49,0.55,0.61,0.68,0.74,0.8,0.87,0.93,1]
y_average_mod_secant_e1 = [2.05555555555556,7.22222222222222,7.13888888888889,7.77777777777778,8.22222222222222,8.33333333333333,8.25,8.13888888888889,8.25,8.16666666666667,8.52777777777778,8.63888888888889,8.97222222222222,9.02777777777778,9.25,9.41666666666667,9.61111111111111,9.72222222222222,9.69444444444445,9.16666666666667,6.44444444444445,5.91666666666667,7.13888888888889,8.61111111111111,1]
y_min_mod_secant_e1 = [2,6,5,6,7,7,7,7,6,7,7,7,7,7,7,7,7,8,7,8,4,5,5,8,1]
y_max_mod_secant_e1 = [4,11,10,11,10,10,10,10,9,10,10,10,10,11,11,11,12,12,11,11,8,7,8,10,1]

x_mod_secant_e2 = [-0.52,-0.46,-0.39,-0.33,-0.27,-0.2,-0.14,-0.07,-0.01,0.04,0.11,0.17,0.23,0.3,0.36,0.42,0.49,0.55,0.61,0.68,0.74,0.8,0.87,0.93,1]
y_average_mod_secant_e2 = [9.27777777777778,7.72222222222222,7.69444444444445,8.33333333333333,8.61111111111111,9.08333333333333,8.69444444444445,8.63888888888889,8.75,8.63888888888889,8.72222222222222,9.05555555555556,9.38888888888889,9.69444444444445,9.69444444444445,9.77777777777778,10.0555555555556,10.1388888888889,10.1111111111111,9.83333333333333,6.88888888888889,6.5,7.61111111111111,9.05555555555556,1]
y_min_mod_secant_e2 = [4,7,6,7,7,8,7,8,7,7,7,8,8,7,7,7,8,8,7,8,5,5,5,8,1]
y_max_mod_secant_e2 = [11,12,11,12,10,11,10,10,10,10,10,10,11,11,11,12,12,12,12,11,8,8,9,10,1]

x_brentq_e2 = [-0.52,-0.46,-0.39,-0.33,-0.27,-0.2,-0.14,-0.07,-0.01,0.04,0.11,0.17,0.23,0.3,0.36,0.42,0.49,0.55,0.61,0.68,0.74,0.8,0.87,0.93,0.99,1]
y_average_brentq_e2 = [12.3888888888889,9,8.69444444444445,7.36111111111111,7.36111111111111,8.25,8.38888888888889,8.13888888888889,8.19444444444445,8.80555555555556,8.77777777777778,9.25,9.91666666666667,9.77777777777778,10.0277777777778,10.1944444444444,10.1944444444444,10.2777777777778,10.0833333333333,9.83333333333333,10,10,10.5277777777778,9.05555555555556,19.4117647058824,20]
y_min_brentq_e2 = [3,7,8,6,6,8,7,7,7,8,8,8,9,8,8,8,8,9,8,9,9,9,10,8,14,20]
y_max_brentq_e2 = [15,11,11,8,9,9,9,10,9,9,10,11,11,11,11,12,12,11,11,11,15,11,11,10,23,20]

x_bisect_e2 = [-0.52,-0.46,-0.39,-0.33,-0.27,-0.2,-0.14,-0.07,-0.01,0.04,0.11,0.17,0.23,0.3,0.36,0.42,0.49,0.55,0.61,0.68,0.74,0.8,0.87,0.93,0.99,1]
y_average_bisect_e2 = [18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18]
y_min_bisect_e2 = [18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18]
y_max_bisect_e2 = [18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18]


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
