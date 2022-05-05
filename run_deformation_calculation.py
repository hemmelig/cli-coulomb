# coding: utf-8
"""
Wrapper script for running the MATLAB Okada source codes provided with Coulomb.

Author: Conor Bacon
Last edited: 05/05/2022
"""

import pathlib

from matlab import engine
import pyproj

import utilities as utils


# Define a projection
proj = pyproj.Proj(
    proj="utm", zone="20W", north=True, ellps="WGS84", units="m", no_defs=True
)

# Set a reference point
refx, refy = proj(-17.6, 64.5)

p = pathlib.Path.cwd()
(p / "strain_files").mkdir(parents=True, exist_ok=True)
(p / "SHmax_files").mkdir(parents=True, exist_ok=True)

session = engine.start_matlab()
session.addpath(session.genpath("okada_source"), nargout=0)

for depth in range(0, 16):
    session.calculate_deformation(
        str(p / "input_files/askja.inp"),
        str(p / f"strain_files/deformation_{depth}km.strain"),
        depth,
        "strain",
        nargout=0
    )
    session.calculate_SHmax(
        str(p / "input_files/askja.inp"),
        f"SHmax_files/SHmax_{depth}km.txt",
        depth,
        nargout=0
    )

    # Convert SHmax .txt file to .xy
    utils.shmax_txt2xy(p / f"SHmax_files/SHmax_{depth}km.txt", proj, refx, refy)
