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
from section import Section, Reinforcement # noqa: E402
import pandas as pd # noqa: E402

def example5():
    section = Section()
    reinforcement = Reinforcement(25.4)
    a = section.create_circle(400, 350, 350)
    r = reinforcement.build_rectangle_rebar(700, 700, 5, 5, 5, 5, 64)
    return ConcreteInteractionDiagram(
        concrete_boundary = pd.DataFrame({
            'x': [0, 0, 700, 700],
            'y': [0, 700, 700, 0]
        }),
        rebar = r,
        openings = [a],
        fc_ = 25,
        epsilon_0 = 0.002,
        epsilon_u = 0.003,
        fy = 400,
        Es = 200000,
        transverse="ties",
        loads = {
            "Pu": [10],
            "Mux": [10],
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
