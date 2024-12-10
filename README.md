# cosmo-ug-dev
----
## **Installation Instructions**

### 1. **Set up the environment**
To use the scripts create the required environment, following these steps. (You don't have to create an environment, just run the code below)

```bash
module load python
source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
```
-----

### 2. **Ask for 4 hours of gpu at NERSC**
Once you activate the environment, you can request 4 hours of NERSC GPU time. With this, the emulator takes approximately 20 minutes to run for all tracers.

```bash
salloc --nodes 1 --qos interactive --time 4:00:00 --constraint gpu --gpus 4 --account=desi
```

### 3. **Run the desired script**

```bash
srun -n 4 python script_name.py
```

The chains will be saved in a folder named 'Chains'. You can find a notebook for plotting the chains in the 'Plots' folder. 
