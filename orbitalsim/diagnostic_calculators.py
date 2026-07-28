import copy
from pathlib import Path
from astroquery.jplhorizons import Horizons
import numpy as np
from astropy.time import Time, TimeDelta
import spiceypy as sp

from orbitalsim.calculations import mechanical_energy_calculator, momentum_calculator, update_extrema
from orbitalsim.simulation import NBodySimulation
from orbitalsim.universal_constants import UniversalConstants

PCK_KERNEL_PATH = Path(__file__).resolve().parents[1] / "pck00010.tpc"
GM_KERNEL_PATH = Path(__file__).resolve().parents[1] / "gm_de431.tpc"

def _to_time(time):
    if isinstance(time, Time):
        return time
    if isinstance(time, (int, float)):
        return Time(time, format='jd', scale='tdb')
    return Time(time)


def two_simulation_position_comparison(bodies, steps, integrator1, integrator2, dt1, dt2):
    sim1 = NBodySimulation(integrator1, copy.deepcopy(bodies))
    sim2 = NBodySimulation(integrator2, copy.deepcopy(bodies))

    _, positions_1, *_ = sim1.simulator(steps,dt1)
    _, positions_2, *_ = sim2.simulator(int(round(steps * dt1/dt2)),dt2)

    return np.linalg.norm(positions_1 - positions_2, axis=1)

def horizon_data_position_comparison(bodies, start_time, steps, integrator, dt):

    if sp.ktotal("ALL") == 0:
        sp.furnsh(str(PCK_KERNEL_PATH))
        sp.furnsh(str(GM_KERNEL_PATH))


    bodies_copy = copy.deepcopy(bodies)

    for body in bodies_copy:
        if body.barycenter_id is not None:
            _, gm = sp.bodvcd(int(body.barycenter_id), "GM", 1)
            gm_si = gm[0] * 1e9
            body.mass = gm_si / UniversalConstants.G

    sim = NBodySimulation(integrator, copy.deepcopy(bodies_copy))

    start_time_converted = _to_time(start_time)
    deltaTime = TimeDelta(steps * dt, format='sec')
    end_time = start_time_converted + deltaTime

    _, positions, *_ = sim.simulator(steps,dt)

    objs = []
    for body in bodies_copy:
        if body.barycenter_id is not None:
            objs.append(Horizons(id=body.barycenter_id, location='@0', epochs=end_time.tdb.jd))
        else:
            objs.append(Horizons(id=body.horizons_id, location='@0', epochs=end_time.tdb.jd))

    positional_differences = np.zeros(len(objs))
    for i in range(len(objs)):
        vec = objs[i].vectors()
        horizons_position = np.array([vec['x'][0], vec['y'][0], vec['z'][0]]) * UniversalConstants.AU_TO_M
        positional_differences[i] = np.linalg.norm(positions[i] - horizons_position)

    return positional_differences
