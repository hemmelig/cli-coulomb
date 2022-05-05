# -*- coding: utf-8 -*-
"""
Utility functions to help with running Okada deformation modelling in MATLAB
from Python.

"""

import pandas as pd

def shmax_txt2xy(text_file, proj, refx, refy):
    """
    Convert .txt file containing SHmax vectors in grid coords to lon/lat .xy for
    plotting in GMT.

    Parameters
    ----------
    text_file : pathlib.Path object
        Path to .txt file containing SHmax vectors.
    proj : pyproj.Proj object
        Map projection in which to convert coordinates.
    refx : float
        Reference position in x dimension.
    refy : float
        Reference position in y dimension.

    """

    vectors = pd.read_csv(text_file, header=None, delim_whitespace=True)
    xy_file = text_file.with_suffix(".xy")
    with open(xy_file, "w") as f:
        for i in range(len(vectors)):
            if i % 3 == 1 or i % 3 == 2:
                pass
            else:
                # Convert vector to lat/lon
                p11, p12 = (float(vectors.iloc[i][0])*1000 + refx,
                            float(vectors.iloc[i][1])*1000 + refy)
                p21, p22 = (float(vectors.iloc[i+1][0])*1000 + refx,
                            float(vectors.iloc[i+1][1])*1000 + refy)
                lon1, lat1 = proj(p11, p12, inverse=True)
                lon2, lat2 = proj(p21, p22, inverse=True)
                f.write(f">\n{lon1} {lat1}\n{lon2} {lat2}\n")
