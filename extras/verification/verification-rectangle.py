# Add parent directory to path so that biaxial_analysis can be imported.
# This is necessary because biaxial_analysis is not installed as a package.
import os
import sys
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir) # traverse one level up
parent_dir = os.path.dirname(parent_dir) # traverse one level up
sys.path.append(parent_dir)

from biaxial_analysis import ConcreteInteractionDiagram, run # noqa: E402
import pandas as pd # noqa: E402

def example5():
    return ConcreteInteractionDiagram(
        concrete_boundary = pd.DataFrame({
            'x': [0, 0, 457.2, 457.2],
            'y': [0, 304.8, 304.8, 0]
        }),
        rebar = pd.DataFrame({
            'diameter': 12*[25.4],
            'x': [60.452, 127, 330.2, 396.748,
                  60.452,             396.748,
                  60.452,             396.748,
                  60.452, 127, 330.2, 396.748
                  ],
            'y': [60.96, 60.96, 60.96, 60.96,
                  121.92, 121.92,
                  182.88, 182.88,
                  243.84, 243.84, 243.84, 243.84, 
                  ],
        }),
        fc_ = 27.58,
        epsilon_0 = 0.002,
        epsilon_u = 0.0038,
        fy = 413.69,
        Es = 200000,
        transverse="ties",
        loads = {
            "Pu": [10],
            "Mux": [10],
            "Muy": [25.1],
        },
        method = "secant",
        tol = 1e-2
    )


if __name__ == "__main__":
    example = example5()
    example.check_clear_concrete_cover()
    run(example)