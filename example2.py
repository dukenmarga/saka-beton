import pandas as pd
from biaxial_analysis import ConcreteInteractionDiagram, run

def example2():
    return ConcreteInteractionDiagram(
        concrete_boundary = pd.DataFrame({
            'x': [0, -216.5, -216.5, 0, 216.5, 216.5],
            'y': [250, 125, -125, -250, -125, 125]
        }),
        rebar = pd.DataFrame({
            'diameter': [16, 16, 16, 16, 16, 16, 16, 16],
            'x': [180, 127.28, 0, -127.28, -180, -127.28, 0 , 127.8],
            'y': [0, 127.28, 180, 127.28, 0, -127.28, -180, -127.28]
        }),
        fc_ = 27.58,
        epsilon_0 = 0.002,
        epsilon_u = 0.0038,
        gamma = 0.15,
        fy = 413.69,
        Es = 200000,
        eps_sh = 0.015,
        eps_su = 0.05,
        transverse="spiral",
        loads = {
            "Pu": [1500, 1100, -300],
            "Mux": [120, 210, 200],
            "Muy": [170, 150, 400],
        }
    )

if __name__ == "__main__":
    run(example2())