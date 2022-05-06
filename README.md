# Python wrapper to run deformation modelling code based on Okada, 1992
A small wrapper to run the deformation modelling code batch style, rather than using GUI provided by Coulomb 3.3.

The core implementation of the Okada 1992 equations have been copied directly from Coulomb, which is described in detail here:

https://pubs.usgs.gov/of/2011/1060/

I then took some of the `.m` scripts for running Coulomb through the GUI and refactored them.

I could not find a proper license file with the distribution, or find any reference to a software license in the manual.

For now, this wrapper only performs stress/strain modelling. Resolution of Coulomb stresses on receiver faults may be added in the future.

Installation
------------
1. Clone this repo - `git clone https://github.com/hemmelig/cli-coulomb` - and navigate to it
2. Install MATLAB - unfortunately, this is proprietary, licensed software! :(
3. Install the packages listed in the environment.yml file, either manually, or using (for example) conda:

```
conda env create -f environment.yml
```
4. Activate your environment, and install the MATLAB Python engine - on *nix systems, this is done as (where `MATLAB_ROOT` should be replaced with the full path to MATLAB):
```
cd "MATLAB_ROOT/extern/engines/python"
python setup.py install
```

Workflow overview
-----------------
1. Create a Coulomb-compatible input file
2. Edit `run_deformation_calculation.py` to cover the depth range of interest

Internally, this script loops over depths, performing the following actions at each depth:
1. Run matlab `calculate_deformation.m` to get the strain field for a given model at a given depth
2. Run matlab `calculate_SHmax.m` to get the orientation of SHmax
3. Convert SHmax `.txt` file to gmt-plottable `.xy` file
4. Convert deformation `.strain` file to gmt-plottable `.grd` file
5. Save everything to be plotted.

An example of a figure showing the outputs of this workflow can be found here - https://github.com/hemmelig/2021JB022655 (Figure 8).

License
-------
This collection is shared here "as is", with no license. Do with it what you will.

Contact
-------
Any additional comments/questions can be directed to:
* **Conor Bacon** - conor.bacon@esc.cam.ac.uk
