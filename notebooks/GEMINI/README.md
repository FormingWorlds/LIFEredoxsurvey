## GEMINI 
**PROTEUS SUBMODULE FOR OBSERVATIONAL SIMULATION OF EXOPLANETS**

GEMINI uses the simulator tool LIFEsim of the LIFE collaboration (Large Interferometer For Exoplanets) to simulate observation of exoplanets in the mid-infrared region.

>[!IMPORTANT]
>Running this module requires the user to have the following Python packages installed:
>
>**lifesim**: to install, follow https://lifesim.readthedocs.io
>
>**spectres**: https://pypi.org/project/spectres/
>
>**configparser**: https://pypi.org/project/configparser/

### Contributors

LC - Lorenzo Cesario (l.cesario@rug.nl)

### Repository Structure

`README.md` - Overview file

`lifesim_setup.py` - LIFE instrument set-up in LIFEsim

`functions.py` - functions for relevant spectrum units conversion and (if requested) the planet's blackbody

`init_simulation.cfg` - configuration file of inputs for simulation parameters

`gemini.py` - main GEMINI Python script

`outputdir` - output directory for plots and text files
