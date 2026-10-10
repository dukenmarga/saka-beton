import matplotlib.pyplot as plt
from itertools import chain
from pathlib import Path


x_mod_secant_e1 = [-0.43,-0.37,-0.31,-0.25,-0.19,-0.13,-0.07,-0.01,0.04,0.1,0.16,0.22,0.28,0.34,0.4,0.46,0.52,0.58,0.64,0.7,0.76,0.82,0.88,0.94,1]
y_average_mod_secant_e1 = [2.16666666666667,7.30555555555556,7.63888888888889,7.80555555555556,8.52777777777778,8.58333333333333,8.47222222222222,8.25,8.38888888888889,8.11111111111111,7.91666666666667,7.69444444444445,7.77777777777778,7.61111111111111,7.58333333333333,7.91666666666667,8.02777777777778,8.25,8.83333333333333,8.66666666666667,8.52777777777778,7.91666666666667,7.13888888888889,6.27777777777778,1]
y_min_mod_secant_e1 = [2,5,6,5,7,6,7,7,7,7,6,6,6,6,6,6,6,6,6,6,5,5,5,5,1]
y_max_mod_secant_e1 = [4,18,11,11,10,11,11,10,12,10,11,11,12,10,9,10,10,10,11,11,13,10,9,9,1]

x_mod_secant_e2 = [-0.43,-0.37,-0.31,-0.25,-0.19,-0.13,-0.07,-0.01,0.04,0.1,0.16,0.22,0.28,0.34,0.4,0.46,0.52,0.58,0.64,0.7,0.76,0.82,0.88,0.94,1]
y_average_mod_secant_e2 = [7.86111111111111,7.75,8.05555555555556,8.33333333333333,8.86111111111111,9,8.86111111111111,8.88888888888889,8.72222222222222,8.58333333333333,8.27777777777778,8.08333333333333,8.13888888888889,8.05555555555556,7.91666666666667,8.22222222222222,8.30555555555556,8.75,9.19444444444445,9.19444444444445,9.02777777777778,8.36111111111111,7.55555555555556,6.94444444444445,1]
y_min_mod_secant_e2 = [3,6,6,6,7,7,7,7,7,7,6,6,6,6,6,7,6,7,7,6,6,6,5,5,1]
y_max_mod_secant_e2 = [9,19,12,11,11,12,11,11,12,11,11,11,12,10,10,10,10,10,11,12,14,10,9,10,1]

x_brentq_e2 = [-0.43,-0.37,-0.31,-0.25,-0.19,-0.13,-0.07,-0.01,0.04,0.1,0.16,0.22,0.28,0.34,0.4,0.46,0.52,0.58,0.64,0.7,0.76,0.82,0.88,0.94,1]
y_average_brentq_e2 = [9.27777777777778,8.61111111111111,7.91666666666667,7.30555555555556,6.80555555555556,7.52777777777778,8.05555555555556,7.88888888888889,7.58333333333333,7.88888888888889,8.11111111111111,8.5,9.08333333333333,9.16666666666667,9.44444444444445,9.55555555555556,9.72222222222222,9.88888888888889,10.1111111111111,10.1944444444444,10.1388888888889,10.5555555555556,10.6944444444444,11.3611111111111,1]
y_min_brentq_e2 = [3,8,6,6,5,6,6,7,6,7,6,7,8,7,8,8,7,7,8,8,8,8,8,9,1]
y_max_brentq_e2 = [11,10,9,9,9,10,10,9,10,10,10,10,11,10,11,11,11,11,11,13,13,12,12,13,1]

x_bisect_e2 = [-0.43,-0.37,-0.31,-0.25,-0.19,-0.13,-0.07,-0.01,0.04,0.1,0.16,0.22,0.28,0.34,0.4,0.46,0.52,0.58,0.64,0.7,0.76,0.82,0.88,0.94,1]
y_average_bisect_e2 = [12.7777777777778,17.4166666666667,17.3888888888889,17.6666666666667,17.3333333333333,17.6666666666667,17.75,17.4166666666667,17.7777777777778,17.7222222222222,17.8055555555556,17.7222222222222,17.8333333333333,17.5,17.6388888888889,17.4722222222222,17.5833333333333,17.75,17.5,17.3888888888889,17.3055555555556,17.2777777777778,17.3611111111111,16.4444444444444,1]
y_min_bisect_e2 = [10,14,8,16,11,13,15,12,16,16,15,14,16,14,13,14,13,16,12,14,13,14,16,11,1]
y_max_bisect_e2 = [17,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,1]


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
