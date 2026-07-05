from orbitalsim.universal_constants import UniversalConstants
from orbitalsim.simulation import  calculate_accelerations

class VelocityVerletIntegrator:

    @staticmethod
    def step(bodies, positions, velocities):

        accelerations = calculate_accelerations(bodies, positions)

        positions += velocities * UniversalConstants.dt + 0.5 * accelerations * UniversalConstants.dt**2

        new_accelerations = calculate_accelerations(bodies, positions)

        velocities += 0.5 * (accelerations + new_accelerations) * UniversalConstants.dt

        return positions, velocities