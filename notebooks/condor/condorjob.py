
import subprocess

inputfolder = input
outputfolder = output

atmosphere_key = input.split("/")[-1]
splits = atmosphere_key.split("_")
distance = int(splits[-2])
config_ARCiS_loc = "/Users/users/borgmann/Documents/masterproject/LIFEredoxsurvey/config/config_ARCiS.in"
save_location = output
atmosphere_file_loc = input

command = f"ARCiS {config_ARCiS_loc} -o {save_location} -s obs1:file={atmosphere_file_loc} -s distance={distance}"
prompt = (
f"module unload ARCiS;"
f"module load ARCiS;"
f"module unload cfitsio;"
f"module load cfitsio;"
f"{command}"
)
subprocess.run(prompt, shell=True, check=True) 