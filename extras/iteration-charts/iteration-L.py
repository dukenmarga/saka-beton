import matplotlib.pyplot as plt
from itertools import chain
from pathlib import Path


x_mod_secant_e1 = [-0.47,-0.41,-0.35,-0.28,-0.22,-0.16,-0.1,-0.04,0.01,0.07,0.14,0.2,0.26,0.32,0.38,0.44,0.5,0.57,0.63,0.69,0.75,0.81,0.87,0.93,1]
y_average_mod_secant_e1 = [2.13888888888889,7.41666666666667,7.44444444444445,7.61111111111111,8.52777777777778,8.63888888888889,8.58333333333333,8.30555555555556,8.33333333333333,8.08333333333333,7.80555555555556,7.77777777777778,7.83333333333333,7.58333333333333,7.61111111111111,7.77777777777778,7.88888888888889,8.30555555555556,8.97222222222222,8.72222222222222,8.41666666666667,8.02777777777778,7.08333333333333,6.25,1]
y_min_mod_secant_e1 = [2,5,5,6,7,6,7,7,7,7,6,7,6,6,6,6,6,6,6,6,5,5,5,5,1]
y_max_mod_secant_e1 = [4,17,11,9,11,10,12,11,10,10,12,11,12,10,10,10,10,10,11,10,10,12,8,9,1]

x_mod_secant_e2 = [-0.47,-0.41,-0.35,-0.28,-0.22,-0.16,-0.1,-0.04,0.01,0.07,0.14,0.2,0.26,0.32,0.38,0.44,0.5,0.57,0.63,0.69,0.75,0.81,0.87,0.93,1]
y_average_mod_secant_e2 = [8,7.80555555555556,7.94444444444445,8.19444444444445,9,9.13888888888889,8.91666666666667,8.80555555555556,8.86111111111111,8.69444444444445,8.30555555555556,8.16666666666667,8.25,8.08333333333333,8.05555555555556,8.13888888888889,8.19444444444445,8.77777777777778,9.44444444444445,9.19444444444445,8.97222222222222,8.44444444444445,7.52777777777778,6.97222222222222,1]
y_min_mod_secant_e2 = [4,6,6,7,7,7,7,7,8,7,7,7,6,6,6,7,6,7,7,6,6,6,5,5,1]
y_max_mod_secant_e2 = [9,17,11,10,11,11,12,11,11,11,12,11,12,10,10,10,10,10,12,11,11,12,9,10,1]

x_brentq_e2 = [-0.47,-0.41,-0.35,-0.28,-0.22,-0.16,-0.1,-0.04,0.01,0.07,0.14,0.2,0.26,0.32,0.38,0.44,0.5,0.57,0.63,0.69,0.75,0.81,0.87,0.93,1]
y_average_brentq_e2 = [11.3611111111111,9.52777777777778,8.77777777777778,7.94444444444445,7.44444444444445,8.30555555555556,8.63888888888889,8.58333333333333,8.36111111111111,8.36111111111111,8.88888888888889,9.27777777777778,9.69444444444445,9.97222222222222,10.2222222222222,10.3611111111111,10.4722222222222,10.6388888888889,10.8333333333333,10.8333333333333,11.1388888888889,11.2777777777778,11.6111111111111,12.0277777777778,28]
y_min_brentq_e2 = [3,8,7,6,6,7,7,7,7,7,8,8,8,8,9,8,8,9,9,9,9,9,9,9,26]
y_max_brentq_e2 = [14,12,11,10,10,10,10,10,10,10,10,11,11,11,11,12,11,12,12,12,13,13,13,13,30]

x_bisect_e2 = [-0.47,-0.41,-0.35,-0.28,-0.22,-0.16,-0.1,-0.04,0.01,0.07,0.14,0.2,0.26,0.32,0.38,0.44,0.5,0.57,0.63,0.69,0.75,0.81,0.87,0.93,0.99]
y_average_bisect_e2 = [18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18]
y_min_bisect_e2 = [18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18]
y_max_bisect_e2 = [18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18]


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

ax.set_title("L-shape Section - Iteration Count vs Axial Load Ratio")
ax.set_xlabel("Pn-target / Po")
ax.xaxis.set_label_coords(0.3, -0.1)
ax.set_ylabel("Iteration to convergence")

ax.legend(loc='upper right', ncols=1, bbox_to_anchor=(1, -0.08), framealpha=0.4)
plt.tight_layout()

script_dir = Path(__file__).resolve().parent
save_path = script_dir / "iteration-L.svg"
plt.savefig(save_path, bbox_inches="tight", pad_inches=0)

plt.show()
