from orbitalsim.universal_constants import UniversalConstants
from orbitalsim.simulation import  calculate_accelerations

# poor for relativity considerations

class VelocityVerletIntegrator:

    @staticmethod
    def step(positions, velocities, GMs, dt):

        accelerations = calculate_accelerations(positions, GMs, velocities)

        positions += velocities * dt + 0.5 * accelerations * dt**2

        new_accelerations = calculate_accelerations(positions, GMs, velocities)

        velocities += 0.5 * (accelerations + new_accelerations) * dt

        return positions, velocities