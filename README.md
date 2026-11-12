This is the code repository used for the MSc thesis 'Detecting atmospheric trends due to mantle redox states of rocky exoplanets orbiting nearby Sun-like stars using LIFE' written by Sanne Borgmann.

--Output locations--
All important output data gets written to the output folder. A lot of software packages are used, and sometimes it's necessary to copy in outputs from the software's own output directory. Because of this, some software installation locations need to be specified in config/masterconfig.yaml.
config/masterconfig.yaml also lets you define major input settings that are used in the notebooks.

--Workflow--
1. EVo_run.ipynb - This uses the EVo package written by Pip Liggins to simulate atmospheres. The masterconfig file defines what fo2 and C/H values are used to generate the atmospheres. This in turn determines how many different atmospheres are generates (default value 12). chem.yaml, env.yaml and output.yaml are also used by EVo, but you shouldn't have to change anything in here.
2. agni_run.ipynb - This part of the pipeline has deprecated in a weird way. AGNI was used in the past as the forward model (and also applies FastChem in the process). The forward model was later changed to ARCiS. However, AGNI is still run to apply FastChem. This means a bunch of compute power is wasted on running AGNI while you could just only use FastChem instead. I didn't get around to changing this. So now only the atmospheric contents from AGNI are used.
3. ARCiS_forward_model.ipynb - This takes the AGNI atmospheric compositions and generates emission spectra. Note that observing distance is an important input here (default = 5, 10, 15, 20 pc). I do all ARCiS runs on an external computer cluster (referred to as norma2 in the code) through SSH. This is because (especially the retrievals) can take very long. You will probably have to slightly rewrite this code to have it run on whatever machine you want.

--Extra notebooks--
In the Extra folder there are two notebooks, p-pop_run.ipynb and PT_profiles.ipynb. These notebooks are not used for the thesis. The PT profiles provide extra insight but were not discussed in the thesis. p-pop is called in survey.ipynb and it is not necessary to run p-pop_run.ipynb first. The notebooks outputs some extra interesting images.

--Acknowledgements--
1. Please note that the GEMINI folder has been copied in from https://github.com/FormingWorlds/GEMINI.
I have not written any of the code in this folder. Copying it was necessary because it is not an installable package. The GEMINI wrapper is written by Lorenzo Cesario.
