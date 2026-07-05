from orbitalsim.universal_constants import UniversalConstants
from orbitalsim.simulation import calculate_accelerations

class ForwardEuler:

    @staticmethod
    def step(bodies, positions, velocities):

        accelerations = calculate_accelerations(bodies, positions)

        positions += UniversalConstants.dt * velocities

        velocities += UniversalConstants.dt * accelerations

        return positions, velocities