from dataclasses import dataclass

@dataclass
class Capacity:
    # Whether the capacity values has been reduced (e.g multiplied by phi)
    reduced: bool
    # Capacity values
    capacity: list