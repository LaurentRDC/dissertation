"""
Generate Debye-Waller factor calculations
"""

import numpy as np
from pathlib import Path
from crystals import Crystal
from mkoneph import MODES, ncells, prepare_modes, debye_waller_factors
from skued import detector_scattvectors
from scipy.interpolate import griddata

INPUT = Path(__file__).parent.parent / "data" / "graphite"
OUTPUT = INPUT / "debye-waller"
OUTPUT.mkdir(exist_ok=True)


def dwaller(qx, qy, reflections, temperatures):
    """
    Calculation of the DW sum
    """
    modes = prepare_modes(reflections)

    Ms = debye_waller_factors(modes, temperatures=temperatures)
    dw = sum(Ms).reshape((-1, 1))

    interpolated = griddata(
        points=modes["LA"].q_points[:, 0:2],
        values=dw,
        xi=(qx, qy),
        method="nearest",
        fill_value=0.0,
    )
    return np.squeeze(interpolated)


if __name__ == "__main__":
    hot = {m: 300 for m in MODES}
    hot["LA"] = 1500

    in_plane_refls = filter(
        lambda tup: tup[2] == 0, Crystal.from_database("C").bounded_reflections(12)
    )
    reflections = list(in_plane_refls)

    qx, qy, _ = detector_scattvectors(
        keV=90, camera_length=0.25, shape=(1024, 1024), pixel_size=28e-6, center=None
    )

    room_temp_dw = dwaller(qx, qy, reflections, temperatures=None)
    hot_dw = dwaller(qx, qy, reflections, temperatures=hot)

    np.save(OUTPUT / "qx.npy", qx)
    np.save(OUTPUT / "qy.npy", qy)
    np.save(OUTPUT / "hot_dw.npy", hot_dw)
    np.save(OUTPUT / "room_temp_dw.npy", room_temp_dw)
