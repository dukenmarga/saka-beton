from biaxial_analysis import ConcreteInteractionDiagram, run
from section import Section, Reinforcement

def example5():
    rebar_dia = 32
    diameter = 508
    n = 8
    clear_cover = 67.056
    section = Section()
    reinforcement = Reinforcement(rebar_dia)
    return ConcreteInteractionDiagram(
        concrete_boundary = section.create_circle(diameter),
        rebar = reinforcement.build_circle_rebar(diameter, n, clear_cover),
        fc_ = 35,
        epsilon_0 = 0.002,
        epsilon_u = 0.003,
        gamma = 0.25,
        fy = 413.68,
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