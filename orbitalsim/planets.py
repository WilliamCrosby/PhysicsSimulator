import numpy as np

from orbitalsim.celestial_body import CelestialBody
from orbitalsim.universal_constants import UniversalConstants

class Planets:

    @staticmethod
    def circular_velocities(r: float):
        return np.sqrt(UniversalConstants.G * 1.989e30 / r)

    @staticmethod
    def true_circular_velocities(r: float, m1: float, m2: float):
        return np.sqrt(UniversalConstants.G * (m1 + m2) / r)

    mercury = CelestialBody(
        3.301e23,
        np.array([0.39 * UniversalConstants.AU, 0.0, 0.0]),
        np.array([0.0, circular_velocities(0.39 * UniversalConstants.AU), 0.0]),
        name="Mercury",
        horizons_id='199', barycenter_id='1'
    )

    venus = CelestialBody(
        4.867e24,
        np.array([0.72 * UniversalConstants.AU, 0.0, 0.0]),
        np.array([0.0, circular_velocities(0.72 * UniversalConstants.AU), 0.0]),
        name="Venus",
        horizons_id='299',barycenter_id='2'
    )

    earth = CelestialBody(
        5.972e24,
        np.array([UniversalConstants.AU, 0.0, 0.0]),
        np.array([0.0, true_circular_velocities(UniversalConstants.AU, 5.972e24, 1.32712440018e20 / UniversalConstants.G), 0.0]),
        name="Earth",
        horizons_id = '399', barycenter_id = '3'
    )

    mars = CelestialBody(
        6.417e23,
        np.array([1.52 * UniversalConstants.AU, 0.0, 0.0]),
        np.array([0.0, circular_velocities(1.52 * UniversalConstants.AU), 0.0]),
        name="Mars",
        horizons_id = '499', barycenter_id = '4'
    )

    jupiter = CelestialBody(
        1.898e27,
        np.array([5.20 * UniversalConstants.AU, 0.0, 0.0]),
        np.array([0.0, circular_velocities(5.20 * UniversalConstants.AU), 0.0]),
        name="Jupiter",
        horizons_id = '599', barycenter_id = '5'
    )

    saturn = CelestialBody(
        5.683e26,
        np.array([9.58 * UniversalConstants.AU, 0.0, 0.0]),
        np.array([0.0, circular_velocities(9.58 * UniversalConstants.AU), 0.0]),
        name="Saturn",
        horizons_id = '699', barycenter_id = '6'
    )

    uranus = CelestialBody(
        8.681e25,
        np.array([19.2 * UniversalConstants.AU, 0.0, 0.0]),
        np.array([0.0, circular_velocities(19.2 * UniversalConstants.AU), 0.0]),
        name="Uranus",
        horizons_id = '799', barycenter_id = '7'
    )

    neptune = CelestialBody(
        1.024e26,
        np.array([30.05 * UniversalConstants.AU, 0.0, 0.0]),
        np.array([0.0, circular_velocities(30.05 * UniversalConstants.AU), 0.0]),
        name="Neptune",
        horizons_id = '899', barycenter_id = '8'
    )

    allPlanets = [mercury, venus, earth, mars, jupiter, saturn, uranus, neptune]
