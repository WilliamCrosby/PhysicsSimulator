import numpy as np

from orbitalsim.celestial_body import CelestialBody
from orbitalsim.universal_constants import UniversalConstants


def mechanical_energy_calculator(bodies: list[CelestialBody], positions, velocities):
    total_kinetic_energy = 0.0
    total_gravitational_potential_energy = 0.0

    for i in range(len(bodies)):
        total_kinetic_energy += 0.5 * bodies[i].mass * np.linalg.norm(velocities[i]) ** 2

        for j in range(i + 1, len(bodies)):
            total_gravitational_potential_energy += (
                -UniversalConstants.G
                * bodies[j].mass
                * bodies[i].mass
                / np.linalg.norm(positions[i] - positions[j])
            )

    return total_kinetic_energy + total_gravitational_potential_energy


def momentum_calculator(bodies: list[CelestialBody], positions, velocities):
    total_linear_momentum = np.zeros(3)
    total_angular_momentum = np.zeros(3)

    for i in range(len(bodies)):
        momentum = bodies[i].mass * velocities[i]

        total_linear_momentum += momentum
        total_angular_momentum += np.cross(positions[i], momentum)

    return total_linear_momentum, total_angular_momentum


def update_extrema(
    bodies,
    positions,
    velocities,
    maximum_energy,
    minimum_energy,
    maximum_linear_momentum,
    minimum_linear_momentum,
    maximum_angular_momentum,
    minimum_angular_momentum,
):
    energy = mechanical_energy_calculator(bodies, positions, velocities)
    linear_momentum, angular_momentum = momentum_calculator(bodies, positions, velocities)
    linear_momentum_norm = np.linalg.norm(linear_momentum)
    angular_momentum_norm = np.linalg.norm(angular_momentum)

    maximum_energy = max(maximum_energy, energy)
    minimum_energy = min(minimum_energy, energy)
    maximum_linear_momentum = max(maximum_linear_momentum, linear_momentum_norm)
    minimum_linear_momentum = min(minimum_linear_momentum, linear_momentum_norm)
    maximum_angular_momentum = max(maximum_angular_momentum, angular_momentum_norm)
    minimum_angular_momentum = min(minimum_angular_momentum, angular_momentum_norm)

    return (
        maximum_energy,
        minimum_energy,
        maximum_linear_momentum,
        minimum_linear_momentum,
        maximum_angular_momentum,
        minimum_angular_momentum,
    )
