from biaxial_analysis import ConcreteInteractionDiagram, run
from section import Section, Reinforcement

def example5():
    section = Section()
    reinforcement = Reinforcement(25)
    return ConcreteInteractionDiagram(
        concrete_boundary = section.create_rectangle(500, 800),
        rebar = reinforcement.build_rectangle_rebar(
            500, 800, 6, 5, 6, 5, 62.5
        ),
        fc_ = 35,
        epsilon_0 = 0.002,
        epsilon_u = 0.003,
        gamma = 0.25,
        fy = 390,
        Es = 200000,
        transverse="ties",
        loads = {
            "Pu": [3700],
            "Mux": [900],
            "Muy": [0],
        }
    )

if __name__ == "__main__":
    example = example5()
    example.check_clear_concrete_cover()
    run(example)