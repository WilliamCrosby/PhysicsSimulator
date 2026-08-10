import numpy as np

from orbitalsim.calculations import mechanical_energy_calculator, momentum_calculator, update_extrema
from orbitalsim.celestial_body import CelestialBody
from orbitalsim.universal_constants import UniversalConstants

def calculate_accelerations(positions, masses):
    accelerations = np.zeros_like(positions)
    n = len(masses)
    G = UniversalConstants.G

    for i in range(n):
        for j in range(n):
            if i == j:
                continue

            r = positions[j] - positions[i]
            dist_sq = np.dot(r, r)

            accelerations[i] += (G * masses[j] * r / (dist_sq * np.sqrt(dist_sq)))

    return accelerations

class NBodySimulation:

    def __init__(self, integrator, bodies: list[CelestialBody]):
        self.integrator = integrator
        self.bodies = bodies

    def simulator(self, steps, dt):

        positions = np.array([b.position.copy() for b in self.bodies], dtype = float)
        velocities = np.array([b.velocity.copy() for b in self.bodies], dtype = float)
        masses = [body.mass for body in self.bodies]

        trajectories = np.empty((steps + 1, len(self.bodies), positions.shape[1]), dtype = float) # for visuals later
        trajectories[0] = positions

        energy = mechanical_energy_calculator(self.bodies, positions, velocities)
        maximum_energy = energy
        minimum_energy = energy

        lin_momentum, ang_momentum = momentum_calculator(self.bodies, positions, velocities)
        maximum_linear_momentum = np.linalg.norm(lin_momentum)
        minimum_linear_momentum = np.linalg.norm(lin_momentum)
        maximum_angular_momentum = np.linalg.norm(ang_momentum)
        minimum_angular_momentum = np.linalg.norm(ang_momentum)


        # primary loop
        for step in range(steps):

            self.integrator.step(positions, velocities, masses, dt)

            trajectories[step + 1] = positions.copy()

            maximum_energy, minimum_energy, maximum_linear_momentum, minimum_linear_momentum, maximum_angular_momentum, minimum_angular_momentum = update_extrema(self.bodies, positions, velocities, maximum_energy, minimum_energy, maximum_linear_momentum, minimum_linear_momentum, maximum_angular_momentum, minimum_angular_momentum)


        for i, body in enumerate(self.bodies):
                body.position = positions[i]
                body.velocity = velocities[i]

        return trajectories, positions, maximum_energy, minimum_energy, maximum_linear_momentum, minimum_linear_momentum, maximum_angular_momentum, minimum_angular_momentum
