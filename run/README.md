# Stacking analysis to search for DM from EGs


## Running stacking analysis 

### Setup
1. Define which targets to consider (Consider viable sources that pass overlap removal in `../overlap/`)
2. Download Fermi data for each target and save in a dir named with the target name in `../data/`
	- Include 15deg search radius
	- 12 years of Fermi data, MET time range (239557417,691545605)
	- Energy range MeV (1000,1000000)
	- Ensure one spacecraft (SC) file is downloaded, but only one is necessary. This is stored in `../data/` outside of a target directory
3. In each directory with source data, define the file `ft1.fits` that lists all data files
4. If no `ltcube` file exists, generate one 
	- Using the SC file in `../data/` and one example target with known RA/dec, run:
		```
		cd run
		python make_ltcube.py <target> <RA> <dec>
		```
		
	- Save `ltcube` output in `../data/` with the SC file

### Fermipy joint likelihood fitting and SED
1. Move to `run/`

	```cd run```
2. Ensure Fermipy/fermitools are installed/usable (eg. activate a conda environment)
3. Update bash steering file `run_preprocess_sed.sh`
	- Update arrays of targets and their associated RA / dec
4. Update hardcoded file paths in `preprocess_likelihoodFitting.py`
	- Point `homedir` to this git repo space (where `run/` and `data/` are stored)
	- Point paths to diffuse models (for variables `galdiff` and `isodiff`)
5. Update hardcoded paths in `get_sed.py`
	- Point `preprocessing_path` to the `run/` dir
5. Run bash file which calls `preprocess_likelihoodFitting.py` and `get_sed.py`

	```source run_preprocess_sed.sh```
	- All outputs are saved to directories in `run/` named after the target. If the directory does not exist, it will be created
	- These directories are organized as follows:
		- main dir space stores configuration file (<target>.yaml), fermipy log file, and txt file of output fit parameters
		- `output/` contains all output files from fermipy including fit results (*.fits, *.npy), source finding, the resulting SED results, etc
		- `output/plots/` contains debugging plots from likelihood fitting and SED creation

Outputs: target directory containing fit parameters, SED fit, many plots, many fit files

### Consider contribution of DM model and convert SED output from energy/flux space to DM mass/sigmav space. Then stack likelihoods and plot 
1. Create .txt file listing all targets to stack (eg. `stackingTargets.txt`)
2. Update bash steering file `run_interpret.sh` or `run_interpret_test.sh`
	- Update arrays containing targets to stack and their associated black hole masses and distances
3. Update hardcoded file paths in `schnittman_likelihood_DM.py`
	- Point `homedir` to `run/` space
	- Check that SED output .npy files will be found by `sedfits`
4. Update hardcoded file path for `homedir` in `plot_TSmaps_stack.py`
5. Define whether you want to save or view plots with `savePlots` bool in `plot_TSmaps_stack.py`
6. Run bash file which calls `schnittman_likelihood_DM.py` for each individual source and then calls `plot_TSmaps_stack.py` to visualize them all and stack them together

	```source run_interpret.sh```
Outputs: (`schnittman_likelihood_DM.py`) npy files with 2D arrays of dloglike values in DM mass/sigmav space
Outputs: (`plot_TSmaps_stack.py`) plots of TS map (in DM space) for individual targets and for full stack (saved either in individual target dir OR in `run/stack/` for full stack)


### Troubleshooting scripts

```tsmap.py```
	- Creates TS maps in the style of the DM TS maps but for SED outputs (in flux/energy space)
	- Make sure to update hardcoded paths and savePlot variable before running
	- Run for one single input like
		python tsmap.py <target>	

```spectrum/spectrum.py```
	- Visualizes DM flux spectrum (dN/dE) associated with Schnittman model for unbound DM near black holes (can compare to Schnittman paper considering DM energy cutoff E_com>3m_X)
	- Uses parameters from test galaxy NGC4649 - can update this by updating global variable values in header
	- Saves plots (dN/dE spectrum as a function of E and as a function of E/m) to the spectrum/ directory. Can comment out the save and replace with plt.show() if no saving is desired
	- Run with no inputs like
		python spectrum/spectrum.py

### Notes
- Original scripts provided by A. McDaniel at Clemson stored in `clemsonScripts`