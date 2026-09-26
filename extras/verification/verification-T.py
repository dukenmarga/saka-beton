# Add parent directory to path so that biaxial_analysis can be imported.
# This is necessary because biaxial_analysis is not installed as a package.
import os
import sys
import timeit
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir) # traverse one level up
parent_dir = os.path.dirname(parent_dir) # traverse one level up
sys.path.append(parent_dir)

from biaxial_analysis import ConcreteInteractionDiagram, run # noqa: E402
import pandas as pd # noqa: E402

def example5():
    return ConcreteInteractionDiagram(
        concrete_boundary = pd.DataFrame({
            'x': [76.2, 2*76.2, 2*76.2, 3*76.2, 3*76.2, 0,     0,     76.2],
            'y': [0,    0,      114.3,  114.3,  190.5 , 190.5, 114.3, 114.3]
        }),
        rebar = pd.DataFrame({
            'diameter': 18*[9.525],
            'x': [93.6625, 134.9375,
                  93.6625, 134.9375,
                  93.6625, 134.9375,
                  17.4625, 17.4625+38.735, 93.6625, 134.9375, 17.4625+38.735*4, 17.4625+38.735*5,
                  17.4625, 17.4625+38.735, 93.6625, 134.9375, 17.4625+38.735*4, 17.4625+38.735*5,
                  ],
            'y': [17.4625, 17.4625,
                  17.4625+38.89, 17.4625+38.89,
                  17.4625+38.89*2, 17.4625+38.89*2,
                  17.4625+38.89*3, 17.4625+38.89*3, 17.4625+38.89*3, 17.4625+38.89*3, 17.4625+38.89*3, 17.4625+38.89*3, 
                  17.4625+38.89*4, 17.4625+38.89*4, 17.4625+38.89*4, 17.4625+38.89*4, 17.4625+38.89*4, 17.4625+38.89*4, 
                  ],
        }),
        fc_ = 29.3,
        epsilon_0 = 0.002,
        epsilon_u = 0.003,
        fy = 532.28,
        Es = 200000,
        transverse="ties",
        loads = {
            "Pu": [10],
            "Mux": [51.45],
            "Muy": [10],
        },
        method = "secant",
        tol = 1e-2
    )

if __name__ == "__main__":
    example = example5()
    example.check_clear_concrete_cover()
    run(example)

    # Benchmark time
    # example.preanalysis()
    # execution_times = timeit.repeat(
    #     stmt=lambda: example.solve_biaxial_analysis(), 
    #     globals=globals(), 
    #     number=1,
    #     repeat=5
    # )
    # best_time = min(execution_times)
    # average_time = sum(execution_times) / len(execution_times)
    # print(f"All trial execution times: {execution_times}")
    # print(f"Best trial execution time: {best_time:.5f} seconds")
    # print(f"Average trial execution time: {average_time:.5f} seconds")
