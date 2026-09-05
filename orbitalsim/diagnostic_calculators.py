import copy
from pathlib import Path
import numpy as np
import pandas as pd

from orbitalsim.calculations import mechanical_energy_calculator, momentum_calculator, update_extrema
from orbitalsim.simulation import NBodySimulation
from orbitalsim.universal_constants import UniversalConstants

PCK_KERNEL_PATH = Path(__file__).resolve().parents[1] / "pck00010.tpc"
GM_KERNEL_PATH = Path(__file__).resolve().parents[1] / "gm_de431.tpc"

def _to_time(time):
    from astropy.time import Time

    if isinstance(time, Time):
        return time
    if isinstance(time, (int, float)):
        return Time(time, format='jd', scale='tdb')
    return Time(time)

def _horizons_target_id(body):
    target_id = body.horizons_id if body.horizons_id is not None else body.barycenter_id
    if target_id is None:
        raise ValueError(f"{body.name} must define either horizons_id or barycenter_id")
    return target_id

def two_simulation_position_comparison(bodies, steps, integrator1, integrator2, dt1, dt2):
    sim1 = NBodySimulation(integrator1, copy.deepcopy(bodies))
    sim2 = NBodySimulation(integrator2, copy.deepcopy(bodies))

    _, _, positions_1, max_energy1, min_energy1, *_ = sim1.simulator(steps,dt1)
    _, _, positions_2, max_energy2, min_energy2, *_ = sim2.simulator(int(round(steps * dt1/dt2)),dt2)

    return np.round(np.linalg.norm(positions_1 - positions_2, axis=1), 2).tolist(), f"integrator 1 energy fluctuation: {max_energy1 - min_energy1}   "f"integrator 2 energy fluctuation: {max_energy2 - min_energy2}"

def horizon_data_position_comparison(bodies, start_time, steps, integrator, dt):
    from astroquery.jplhorizons import Horizons
    from astropy.time import TimeDelta
    import spiceypy as sp

    if sp.ktotal("ALL") == 0:
        sp.furnsh(str(PCK_KERNEL_PATH))
        sp.furnsh(str(GM_KERNEL_PATH))

    bodies_copy = copy.deepcopy(bodies)
    start_time_converted = _to_time(start_time)
    deltaTime = TimeDelta(steps * dt, format='sec')
    end_time = start_time_converted + deltaTime

    for body in bodies_copy:
        spice_id = _horizons_target_id(body)
        _, gm = sp.bodvcd(int(spice_id), "GM", 1)
        gm_si = gm[0] * 1e9
        body.GM = gm_si
        body.mass = gm_si / UniversalConstants.G

        vec = Horizons(id=spice_id, location='@0', epochs=start_time_converted.tdb.jd).vectors()
        body.position = np.array([vec['x'][0], vec['y'][0], vec['z'][0]]) * UniversalConstants.AU_TO_M
        body.velocity = np.array([vec['vx'][0], vec['vy'][0], vec['vz'][0]]) * UniversalConstants.AU_TO_M / 86400

    sim = NBodySimulation(integrator, bodies_copy)
    _, _, positions, *_ = sim.simulator(steps, dt)

    objs = []
    for body in bodies_copy:
        target_id = _horizons_target_id(body)
        objs.append(Horizons(id=target_id, location='@0', epochs=end_time.tdb.jd))

    positional_differences = np.zeros(len(objs))
    for i in range(len(objs)):
        vec = objs[i].vectors()
        horizons_position = np.array([vec['x'][0], vec['y'][0], vec['z'][0]]) * UniversalConstants.AU_TO_M
        positional_differences[i] = np.linalg.norm(positions[i] - horizons_position)

    return positional_differences

def analyze_two_body(traj, saved_velocities, bodies, dt, orbital_period):

    if len(bodies) != 2:
        raise ValueError("analyze_two_body expects exactly two bodies")

    total_time = (len(traj) - 1) * dt
    number_of_periods = int(total_time / orbital_period)

    initial_relative_position = traj[0, 0] - traj[0, 1]
    initial_relative_velocity = saved_velocities[0, 0] - saved_velocities[0, 1]
    initial_radius = np.linalg.norm(initial_relative_position)
    gravitational_parameter = UniversalConstants.G * (bodies[0].mass + bodies[1].mass)
    analytical_speed = np.sqrt(gravitational_parameter / initial_radius)

    orbital_axis = np.cross(initial_relative_position, initial_relative_velocity)
    orbital_axis_norm = np.linalg.norm(orbital_axis)
    if orbital_axis_norm == 0:
        raise ValueError("initial relative position and velocity cannot be parallel")
    orbital_axis = orbital_axis / orbital_axis_norm

    radial_unit = initial_relative_position / initial_radius
    tangential_unit = np.cross(orbital_axis, radial_unit)
    initial_analytical_velocity = analytical_speed * tangential_unit

    rows = []

    for i in range(number_of_periods + 1):
        k = int(round(i * orbital_period / dt))
        if k >= len(traj):
            break

        time = k * dt
        theta = 2 * np.pi * time / orbital_period
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)

        analytical_relative_position = (
            initial_relative_position * cos_theta
            + np.cross(orbital_axis, initial_relative_position) * sin_theta
            + orbital_axis * np.dot(orbital_axis, initial_relative_position) * (1 - cos_theta)
        )
        analytical_relative_velocity = (
            initial_analytical_velocity * cos_theta
            + np.cross(orbital_axis, initial_analytical_velocity) * sin_theta
            + orbital_axis * np.dot(orbital_axis, initial_analytical_velocity) * (1 - cos_theta)
        )

        radius = np.linalg.norm(traj[k, 0] - traj[k, 1])
        relative_velocity = (saved_velocities[k, 0] - saved_velocities[k, 1])
        relative_position = traj[k, 0] - traj[k, 1]
        radius_error = radius - initial_radius
        position_error = np.linalg.norm(relative_position - analytical_relative_position)
        velocity_vector_error = np.linalg.norm(relative_velocity - analytical_relative_velocity)
        speed_error = np.linalg.norm(relative_velocity) - analytical_speed

        rows.append({
            "period": i,
            "time_years": time / (365.25 * 24 * 3600),
            "radius_meters": radius,
            "radius_error_meters": radius_error,
            "velocity_error_meters": speed_error,
            "speed_error_m_per_s": speed_error,
            "position_error_meters": position_error,
            "velocity_vector_error_m_per_s": velocity_vector_error,
        })

    return pd.DataFrame(rows)
