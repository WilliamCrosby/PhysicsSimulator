from astroquery.jplhorizons import Horizons
import numpy as np
from astropy.time import Time

import orbitalsim.planets
from orbitalsim.celestial_body import CelestialBody
from orbitalsim.universal_constants import UniversalConstants


def _to_horizons_epoch(time):
    if isinstance(time, Time):
        return time.tdb.jd
    if isinstance(time, (int, float)):
        return time
    return Time(time).tdb.jd


class PositionsSetter:

    @staticmethod
    def set_positions_at_time(bodies: list[CelestialBody], start_time):
        epoch = _to_horizons_epoch(start_time)

        for body in bodies:
            if body.horizons_id is not None:
                obj = Horizons(id=body.horizons_id, location='@0', epochs=epoch)
            else:
                obj = Horizons(id=body.barycenter_id, location='@0', epochs=epoch)

            vec = obj.vectors()

            body.position = np.array([vec['x'][0], vec['y'][0], vec['z'][0]]) * UniversalConstants.AU_TO_M
            body.velocity = np.array([vec['vx'][0], vec['vy'][0], vec['vz'][0]]) * UniversalConstants.AU_TO_M / 86400
