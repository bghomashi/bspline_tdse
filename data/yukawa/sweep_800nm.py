import json
import os
import os.path
from os import path

working_directory = os.getcwd()

photons = False
ns = [600, 800]
dts = [0.1, 0.05, 0.02]
Rs = [500, 750]
Ls = [30, 40]

if photons:
  repo_dir = "/home/becker/begh0305/Documents/bspline_code/bspline_tdse"
else:
  repo_dir = "/data/becker/begh0305/Research/bspline_tdse"
params = []

os.system("mkdir -p " + working_directory)
for r in Rs:
  for l in Ls:
    for dt in dts:
      for n in ns:
        dir  = working_directory + "/n_" + str(n) + "_dt_" + str(dt) + "_bs_" + str(r) + "_l_" + str(l) + "/"
        params.append({
          "n" : n, 
          "dt" : dt, 
          "R" : r, 
          "L" : l, 
          "dir" : dir})

if photons:
  params.reverse()

for ip in range(len(params)):
  param = params[ip]
  n = param["n"]
  dt = param["dt"]
  r = param["R"]
  l = param["L"]
  dir = param["dir"]

  # start from working_directory 
  print("mkdir " + dir)
  os.chdir(working_directory)                                    # back to main directory
  os.system("rm " + dir + "input.json")
  os.system("rm -r " + dir)
  os.system("mkdir " + dir)
  os.system("cp base_input.json " + dir + "input.json")
  os.chdir(dir)                        # go to working directory

  # open the input file
  with open('input.json', 'r+') as f:
    data = json.load(f)
    
    ## input file
    data["time_step"] = dt
    data["basis"]["num_nodes"] = n
    data["basis"]["x_max"] = r
    data["basis"]["lmax"] = l

    f.seek(0)        # <--- should reset file position to the beginning.
    json.dump(data, f, indent=4)
    f.truncate()     # remove remaining part

    with open('Queue_h.bash', 'w') as f:
      f.write("#!/usr/bin/bash\n")
      f.write(f'#SBATCH --job-name Ar_conv_ip_{ip:1.3f}\n')
      f.write(f'#SBATCH --output Ar_conv_ip_{ip:1.3f}.o%j\n')

      if photons:
        f.write("#SBATCH --partition=compute --qos=normal\n")
        f.write("#SBATCH --exclude=photon13,photon16,photon11\n")
        f.write("#SBATCH --nodes 1\n")
        f.write("#SBATCH --ntasks 8\n")
        f.write("#SBATCH --nice=1.0\n")
        f.write("#SBATCH --mem=16G\n")
      else:
        f.write("#SBATCH --partition=jila\n")
        f.write("#SBATCH --nodes 1\n")
        f.write("#SBATCH --ntasks 16\n")
        f.write("#SBATCH --nice=1.0\n")
        f.write("#SBATCH --mem=16G\n")
      
      f.write("\n")
      f.write("source ~/.bashrc\n\n")
      f.write("REPO_DIR=\"" + repo_dir + "\"\n")
      f.write("hostname\n")
      f.write("pwd\n")
      f.write("\n")
      f.write("mpirun -np $SLURM_NTASKS $REPO_DIR/bspline_tdse.out tise >> results.log\n")
      f.write("mpirun -np $SLURM_NTASKS $REPO_DIR/bspline_tdse.out >> results.log\n")
    
    # execute job
    os.system("rm *.log")
    os.system("sbatch Queue_h.bash")

    os.chdir(working_directory)

