import copy

import numpy as np

from orbitalsim.calculations import mechanical_energy_calculator, momentum_calculator, update_extrema
from orbitalsim.simulation import NBodySimulation


def two_simulation_position_comparison(bodies, steps, integrator1, integrator2, dt1, dt2):
    sim1 = NBodySimulation(integrator1, copy.deepcopy(bodies))
    sim2 = NBodySimulation(integrator2, copy.deepcopy(bodies))

    _, positions_1, *_ = sim1.simulator(steps,dt1)
    _, positions_2, *_ = sim2.simulator(int(round(steps * dt1/dt2)),dt2)

    return np.linalg.norm(positions_1 - positions_2, axis=1)

def horizon_data_position_comparison():
    # blah
    print("WIP")

