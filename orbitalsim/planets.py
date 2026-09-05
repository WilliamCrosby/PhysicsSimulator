import numpy as np

from orbitalsim.celestial_body import CelestialBody
from orbitalsim.universal_constants import UniversalConstants

SUN_GM = 1.32712440018e20


def gm_from_km3_per_s2(gm: float):
    return gm * 1e9


def mass_from_gm(gm: float):
    return gm / UniversalConstants.G


def circular_velocity(r: float):
    return np.sqrt(SUN_GM / r)


class Planets:

    @staticmethod
    def circular_velocities(r: float):
        return circular_velocity(r)

    @staticmethod
    def true_circular_velocities(r: float, m1: float, m2: float):
        return np.sqrt(UniversalConstants.G * (m1 + m2) / r)

    mercury = CelestialBody(
        gm_from_km3_per_s2(22031.868551),
        mass_from_gm(gm_from_km3_per_s2(22031.868551)),
        np.array([0.39 * UniversalConstants.AU_TO_M, 0.0, 0.0]),
        np.array([0.0, circular_velocity(0.39 * UniversalConstants.AU_TO_M), 0.0]),
        name="Mercury",
        barycenter_id='1'
    )

    venus = CelestialBody(
        gm_from_km3_per_s2(324858.592),
        mass_from_gm(gm_from_km3_per_s2(324858.592)),
        np.array([0.72 * UniversalConstants.AU_TO_M, 0.0, 0.0]),
        np.array([0.0, circular_velocity(0.72 * UniversalConstants.AU_TO_M), 0.0]),
        name="Venus",
        barycenter_id='2'
    )

    earth = CelestialBody(
        gm_from_km3_per_s2(398600.435507),
        mass_from_gm(gm_from_km3_per_s2(398600.435507)),
        np.array([UniversalConstants.AU_TO_M, 0.0, 0.0]),
        np.array([0.0, np.sqrt((SUN_GM + gm_from_km3_per_s2(398600.435507)) / UniversalConstants.AU_TO_M), 0.0]),
        name="Earth",
        horizons_id = '399', barycenter_id = '3'
    )

    moon = CelestialBody(
        gm_from_km3_per_s2(4902.800118),
        mass_from_gm(gm_from_km3_per_s2(4902.800118)),
        np.array([1.009 * UniversalConstants.AU_TO_M, 0.0, 0.0]),
        np.array([0.0, circular_velocity(1.009 * UniversalConstants.AU_TO_M), 0.0]),
        name="Moon",
        horizons_id='301', barycenter_id='3'
    )

    mars = CelestialBody(
        gm_from_km3_per_s2(42828.375816),
        mass_from_gm(gm_from_km3_per_s2(42828.375816)),
        np.array([1.52 * UniversalConstants.AU_TO_M, 0.0, 0.0]),
        np.array([0.0, circular_velocity(1.52 * UniversalConstants.AU_TO_M), 0.0]),
        name="Mars",
        barycenter_id = '4'
    )

    jupiter = CelestialBody(
        gm_from_km3_per_s2(126712764.1),
        mass_from_gm(gm_from_km3_per_s2(126712764.1)),
        np.array([5.20 * UniversalConstants.AU_TO_M, 0.0, 0.0]),
        np.array([0.0, circular_velocity(5.20 * UniversalConstants.AU_TO_M), 0.0]),
        name="Jupiter",
        barycenter_id = '5'
    )

    saturn = CelestialBody(
        gm_from_km3_per_s2(37940584.8418),
        mass_from_gm(gm_from_km3_per_s2(37940584.8418)),
        np.array([9.58 * UniversalConstants.AU_TO_M, 0.0, 0.0]),
        np.array([0.0, circular_velocity(9.58 * UniversalConstants.AU_TO_M), 0.0]),
        name="Saturn",
        barycenter_id = '6'
    )

    uranus = CelestialBody(
        gm_from_km3_per_s2(5794556.4),
        mass_from_gm(gm_from_km3_per_s2(5794556.4)),
        np.array([19.2 * UniversalConstants.AU_TO_M, 0.0, 0.0]),
        np.array([0.0, circular_velocity(19.2 * UniversalConstants.AU_TO_M), 0.0]),
        name="Uranus",
        barycenter_id = '7'
    )

    neptune = CelestialBody(
        gm_from_km3_per_s2(6836527.10058),
        mass_from_gm(gm_from_km3_per_s2(6836527.10058)),
        np.array([30.05 * UniversalConstants.AU_TO_M, 0.0, 0.0]),
        np.array([0.0, circular_velocity(30.05 * UniversalConstants.AU_TO_M), 0.0]),
        name="Neptune",
        barycenter_id = '8'
    )

    allPlanets = [mercury, venus, earth, mars, jupiter, saturn, uranus, neptune]
    solar_system_bodies = [mercury, venus, earth, moon, mars, jupiter, saturn, uranus, neptune]
