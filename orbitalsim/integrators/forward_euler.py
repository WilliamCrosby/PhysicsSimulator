from orbitalsim.universal_constants import UniversalConstants
from orbitalsim.simulation import calculate_accelerations
import orbitalsim.simulation as sim


class ForwardEuler:

    @staticmethod
    def step(positions, velocities, GMs, dt):

        accelerations = calculate_accelerations(positions, GMs)

        positions += dt * velocities

        velocities += dt * accelerations

        return positions, velocities