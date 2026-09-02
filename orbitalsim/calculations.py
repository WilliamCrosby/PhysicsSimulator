import numpy as np

from orbitalsim.celestial_body import CelestialBody
from orbitalsim.universal_constants import UniversalConstants

''' this one may be better for large simulations with a large number of bodies, when we do that we'll test it
def mechanical_energy_calculator(masses, positions, velocities):
    total_kinetic_energy = 0.5 * np.sum(masses * np.sum(velocities * velocities, axis = 1))
    total_gravitational_potential_energy = 0.0

    for i in range(len(masses)):

        for j in range(i + 1, len(masses)):
            total_gravitational_potential_energy += (
                -UniversalConstants.G
                * masses[j]
                * masses[i]
                / np.linalg.norm(positions[i] - positions[j])
            )

    return total_kinetic_energy + total_gravitational_potential_energy
'''

def mechanical_energy_calculator(masses, positions, velocities):
    G = UniversalConstants.G

    kinetic_energy = 0.5 * np.sum(
        masses * np.sum(velocities * velocities, axis=1)
    )

    diff = positions[:, None, :] - positions[None, :, :]
    distances = np.linalg.norm(diff, axis=2)

    i, j = np.triu_indices(len(masses), k=1)

    potential_energy = -G * np.sum(
        masses[i] * masses[j] / distances[i, j]
    )

    return kinetic_energy + potential_energy

def momentum_calculator(masses, positions, velocities):
    momenta = masses[:, None] * velocities

    total_linear_momentum = np.sum(momenta, axis=0)
    total_angular_momentum = np.sum(np.cross(positions, momenta), axis=0)

    return total_linear_momentum, total_angular_momentum


def update_extrema(
    masses,
    positions,
    velocities,
    maximum_energy,
    minimum_energy,
    maximum_linear_momentum,
    minimum_linear_momentum,
    maximum_angular_momentum,
    minimum_angular_momentum,
):
    energy = mechanical_energy_calculator(masses, positions, velocities)
    linear_momentum, angular_momentum = momentum_calculator(masses, positions, velocities)
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
