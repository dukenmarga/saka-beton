import matplotlib.pyplot as plt
from itertools import chain
from pathlib import Path

x_mod_secant_e1 = [-0.486,-0.424,-0.362,-0.3,-0.238,-0.176,-0.114,-0.052,0.009,0.07,0.132,0.194,0.256,0.318,0.38,0.442,0.504,0.566,0.628,0.69,0.752,0.814,0.876,0.938,1]
y_average_mod_secant_e1 = [2.11111111111111,8.33333333333333,7.66666666666667,8.16666666666667,8.66666666666667,9.88888888888889,9.22222222222222,9.11111111111111,8.55555555555556,8.16666666666667,8.16666666666667,8.33333333333333,7.66666666666667,7.77777777777778,7.16666666666667,8.11111111111111,8.27777777777778,11.3888888888889,9.22222222222222,8.94444444444445,8.33333333333333,7.94444444444445,6.88888888888889,6.27777777777778,1]
y_min_mod_secant_e1 = [2,6,6,7,8,8,8,8,7,7,8,7,7,7,7,8,7,7,7,6,6,5,5,5,1]
y_max_mod_secant_e1 = [3,17,10,10,10,11,10,10,9,9,9,9,8,9,8,9,9,31,10,10,10,9,8,7,1]

x_mod_secant_e2 = [-0.486,-0.424,-0.362,-0.3,-0.238,-0.176,-0.114,-0.052,0.009,0.07,0.132,0.194,0.256,0.318,0.38,0.442,0.504,0.566,0.628,0.69,0.752,0.814,0.876,0.938,1]
y_average_mod_secant_e2 = [10.1111111111111,8.72222222222222,8.22222222222222,8.66666666666667,9.16666666666667,10.2222222222222,9.72222222222222,9.5,9.05555555555556,8.61111111111111,8.55555555555556,8.61111111111111,8.38888888888889,8.11111111111111,7.44444444444445,8.16666666666667,8.88888888888889,11.7777777777778,9.38888888888889,9.16666666666667,8.83333333333333,8.33333333333333,7.55555555555556,6.94444444444445,1]
y_min_mod_secant_e2 = [4,6,7,8,8,9,8,9,8,8,8,8,8,8,7,8,8,8,7,7,6,6,6,6,1]
y_max_mod_secant_e2 = [11,18,11,10,10,11,11,10,10,10,9,9,9,9,9,9,9,32,11,10,10,9,8,8,1]

x_brentq_e2 = [-0.48,-0.42,-0.36,-0.3,-0.23,-0.17,-0.11,-0.05,0,0.07,0.13,0.19,0.25,0.31,0.38,0.44,0.5,0.56,0.62,0.69,0.75,0.81,0.87,0.93,0.99,1]
y_average_brentq_e2 = [13.2222222222222,10.3888888888889,9.77777777777778,8.27777777777778,7.66666666666667,8.05555555555556,8.55555555555556,8.55555555555556,8.22222222222222,8.66666666666667,8.94444444444445,9.33333333333333,9.77777777777778,10.1111111111111,10.6666666666667,10.8888888888889,11.0555555555556,10.9444444444444,10.7777777777778,10.8888888888889,11.2777777777778,11.3888888888889,11.4444444444444,10.6111111111111,28.7647058823529,26]
y_min_brentq_e2 = [3,9,8,7,7,7,8,8,7,8,8,9,9,10,10,10,10,10,10,10,9,9,9,10,25,26]
y_max_brentq_e2 = [15,12,11,10,8,9,9,10,10,10,9,10,11,11,11,11,12,11,12,12,12,12,12,12,31,26]

x_bisect_e2 = [-0.48,-0.42,-0.36,-0.3,-0.23,-0.17,-0.11,-0.05,0,0.05,0.07,0.13,0.19,0.25,0.31,0.38,0.44,0.5,0.56,0.62,0.69,0.75,0.81,0.87,0.93,0.99]
y_average_bisect_e2 = [19.4571428571429,19.4444444444444,19.4444444444444,19.4444444444444,19.4444444444444,19.4444444444444,19.4444444444444,19.4444444444444,19.4444444444444,19,19.4444444444444,19.4444444444444,19.4444444444444,19.4444444444444,19.4444444444444,19.4444444444444,19.4444444444444,19.4444444444444,19.4444444444444,19.4444444444444,19.4444444444444,19.4444444444444,19.4444444444444,19.4444444444444,19.4444444444444,19.6153846153846]
y_min_bisect_e2 = [19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19]
y_max_bisect_e2 = [20,20,20,20,20,20,20,20,20,19,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20]


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

ax.set_title("Rectangle Section - Iteration Count and Axial Load Ratio")
ax.set_xlabel("Pn-target / Po")
ax.xaxis.set_label_coords(0.3, -0.1)
ax.set_ylabel("Iteration to convergence")

ax.legend(loc='upper right', ncols=1, bbox_to_anchor=(1, -0.08), framealpha=0.4)
plt.tight_layout()

script_dir = Path(__file__).resolve().parent
save_path = script_dir / "iteration-rectangle.svg"
plt.savefig(save_path, bbox_inches="tight", pad_inches=0)

plt.savefig("iteration-rectangle.svg", bbox_inches="tight", pad_inches=0)

plt.show()
