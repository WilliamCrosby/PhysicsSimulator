from orbitalsim.universal_constants import UniversalConstants
from orbitalsim.simulation import calculate_accelerations

class ForwardEuler:

    @staticmethod
    def step(bodies, positions, velocities, dt):

        accelerations = calculate_accelerations(bodies, positions)

        positions += dt * velocities

        velocities += dt * accelerations

        return positions, velocities