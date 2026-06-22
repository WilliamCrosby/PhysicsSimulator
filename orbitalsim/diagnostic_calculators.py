import numpy as np
from orbitalsim.celestial_body import CelestialBody
from orbitalsim.universal_constants import UniversalConstants


def mechanical_energy_calculator(bodies:list[CelestialBody], positions, velocities):

    gravitational_potential_energies = []
    kinetic_energies = []

    total_kinetic_energy = 0.0
    total_gravitational_potential_energy = 0.0

    for i in range(len(bodies)):

        kinetic_energies.append(0.5 * bodies[i].mass * np.linalg.norm(velocities[i])**2)

        gpe = 0.0

        for j in range(i + 1, len(bodies)):

            if j == i:
                continue

            gpe += (-1 * UniversalConstants.G * bodies[j].mass * bodies[i].mass / np.linalg.norm(positions[i] - positions[j]))

        gravitational_potential_energies.append(gpe)

    for i in range(len(gravitational_potential_energies)):
        total_gravitational_potential_energy += gravitational_potential_energies[i]
    for j in range (len(kinetic_energies)):
        total_kinetic_energy += kinetic_energies[j]

    return total_kinetic_energy + total_gravitational_potential_energy


def momentum_calculator(bodies:list[CelestialBody], positions, velocities):

    total_linear_momentum = np.zeros(3)
    total_angular_momentum = np.zeros(3)

    for i in range(len(bodies)):
        p = bodies[i].mass * velocities[i]

        total_linear_momentum += p
        total_angular_momentum += np.cross(positions[i], p)

    return total_linear_momentum, total_angular_momentum


def update_extrema(bodies, positions, velocities, maximum_energy, minimum_energy, maximum_linear_momentum, minimum_linear_momentum, maximum_angular_momentum, minimum_angular_momentum):
    if mechanical_energy_calculator(bodies, positions, velocities) > maximum_energy:
        maximum_energy = mechanical_energy_calculator(bodies, positions, velocities)
    if mechanical_energy_calculator(bodies, positions, velocities) < minimum_energy:
        minimum_energy = mechanical_energy_calculator(bodies, positions, velocities)
    if np.linalg.norm(momentum_calculator(bodies, positions, velocities)[0]) > maximum_linear_momentum:
        maximum_linear_momentum = np.linalg.norm(momentum_calculator(bodies, positions, velocities)[0])
    if np.linalg.norm(momentum_calculator(bodies, positions, velocities)[0]) < minimum_linear_momentum:
        minimum_linear_momentum = np.linalg.norm(momentum_calculator(bodies, positions, velocities)[0])
    if np.linalg.norm(momentum_calculator(bodies, positions, velocities)[1]) > maximum_angular_momentum:
        maximum_angular_momentum = np.linalg.norm(momentum_calculator(bodies, positions, velocities)[1])
    if np.linalg.norm(momentum_calculator(bodies, positions, velocities)[1]) < minimum_angular_momentum:
        minimum_angular_momentum = np.linalg.norm(momentum_calculator(bodies, positions, velocities)[1])

    return maximum_energy, minimum_energy, maximum_linear_momentum, minimum_linear_momentum, maximum_angular_momentum, minimum_angular_momentum
