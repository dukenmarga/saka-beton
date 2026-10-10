# Add parent directory to path so that biaxial_analysis can be imported.
# This is necessary because biaxial_analysis is not installed as a package.
import os
import sys
import timeit
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir) # traverse one level up
parent_dir = os.path.dirname(parent_dir) # traverse one level up
sys.path.append(parent_dir)

from biaxial_analysis import ConcreteInteractionDiagram, run  # noqa: E402
from section import Section, Reinforcement # noqa: E402

def example5():
    section = Section()
    reinforcement = Reinforcement(9.525)
    return ConcreteInteractionDiagram(
        concrete_boundary = section.create_l_shape("L-shape-1", 152.4, 190.5, 76.2, 76.2),
        rebar = reinforcement.build_lshape_rebar(
            "L-shape-1", 152.4, 190.5, 76.2, 76.2, 4, 5, 2, 2, 4, 3, 17.4625
        ),
        fc_ = 24.13166,
        epsilon_0 = 0.002,
        epsilon_u = 0.0038,
        fy = 357.148568,
        Es = 200000,
        transverse="ties",
        loads = {
            "Pu": [10],
            "Mux": [10],
            "Muy": [10],
        },
        max_stress_modifier=1,
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
