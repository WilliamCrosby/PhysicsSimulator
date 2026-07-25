import numpy as np

from orbitalsim.calculations import mechanical_energy_calculator, momentum_calculator, update_extrema
from orbitalsim.celestial_body import CelestialBody
from orbitalsim.universal_constants import UniversalConstants

def calculate_accelerations(bodies, positions):
    accelerations = np.zeros_like(positions)

    for i in range(len(bodies)):
        for j in range(len(bodies)):
            if i == j:
                continue

            r = positions[j] - positions[i]

            dist_sq = np.dot(r, r)
            dist = np.sqrt(dist_sq)

            accelerations[i] += (UniversalConstants.G * bodies[j].mass * r / dist ** 3)

    return accelerations

class NBodySimulation:

    def __init__(self, integrator, bodies: list[CelestialBody]):
        self.integrator = integrator
        self.bodies = bodies

    def simulator(self, steps, dt):

        positions = np.array([b.position.copy() for b in self.bodies], dtype = float)
        velocities = np.array([b.velocity.copy() for b in self.bodies], dtype = float)

        trajectories = [positions.copy()] # for visuals later

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

            self.integrator.step(self.bodies, positions, velocities, dt)

            trajectories.append(positions.copy())

            maximum_energy, minimum_energy, maximum_linear_momentum, minimum_linear_momentum, maximum_angular_momentum, minimum_angular_momentum = update_extrema(self.bodies, positions, velocities, maximum_energy, minimum_energy, maximum_linear_momentum, minimum_linear_momentum, maximum_angular_momentum, minimum_angular_momentum)


        for i, body in enumerate(self.bodies):
                body.position = positions[i]
                body.velocity = velocities[i]

        return np.array(trajectories), positions, maximum_energy, minimum_energy, maximum_linear_momentum, minimum_linear_momentum, maximum_angular_momentum, minimum_angular_momentum
