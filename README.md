# openmpFV
A OMP parallel second order accurate FV code to solve the Euler equations on Cartesian and Voronoi meshes. 

Includes grid generation methods from [vmp](https://github.com/lucas56098/voronoi_mesh_project) and a second order slope limited MUSCL Hankock scheme to solve the Euler equations as in ([Springel 2010](https://arxiv.org/abs/0901.4107)). We employ an HLL Riemann solver. Additionally there exist various plotting and analysis routines in a python visualization toolkit.

The original version of this code was part of my bachelor thesis and can be found [here](https://github.com/lucas56098/hydro_bsc_project). This version has various improvements and will recieve further updates in the future.

---
> [!Warning]  
> Construction Zone: Expect large changes in single commits. Thins may break, change or evolve rapidly. No documentation so far...

#### Improvements:
- openMP parallelization of hydro calculations
- CFL timestepping and revised main loop
- general code cleanup + removal of advection, SWE ...
- improved KH and RT initial conditions for direct comparison to ([Springel 2010](https://arxiv.org/abs/0901.4107))
- improved file naming: ```c_n30_FV2_BC1_0_1s_test_step148.csv```
- initial conditions into seperate file
- run & compile options in bash file ```run.sh```
- option to restart code from snapshot (voronoi only so far)
- customizable box length for cartesian and voronoi
- Added HLLC riemann solver instead of HLL
- Getting started guide

#### Todo:
- Long runs: KH1024

#### Further ideas
- Profiling for Ahmdals law calculations (parallel fration calc in optional runs `-t`)
- Profiling for performance improvements
- direct speed comparison with AREPO
- FV moving mesh
- eventually DG improvements

---
### Simulation Examples
KH: cartesian FV 2nd-order, HLLC
<p align="center">
  <img src="/figures/KH_HLLC_50.png" alt="1" width="39%">
  <img src="/figures/KH_HLLC_1024_full_crop.gif" alt="1" width="39%">
</p>


RT: N = 48x144, cartesian FV 2nd-order, HLLC
<p align="center">
  <img src="/figures/rt_ompFV48.png" alt="1" width="80%">
</p>

2D-Riemann Problem as in ([Kurganov and Tadmor, 2002](https://www.semanticscholar.org/paper/Solution-of-two%E2%80%90dimensional-Riemann-problems-for-Kurganov-Tadmor/a44da75f9a36ab879fb9073f2571801eb7bc74a3))
<p align="center">
  <img src="/figures/quadshock2_1024_lres_crop.gif" alt="1" width="60%">
</p>


---
### Getting started
Before starting make sure you have the following installed:

- C++ with a working compiler that supports OpenMP
- CMake
- Git (alternatively one can manually download the files)
- Python packages for ```vis_tk.py```:
```bash
pip install pandas numpy matplotlib tqdm scipy sodshock aeropy geopandas 
```


Start by going into the folder where you want to clone the repository into and do:

```bash
git clone https://github.com/lucas56098/openmpFV.git
```

After that run the install script

```bash
cd openmpFV
chmod +x install.sh
./install.sh
```

This will install the Eigen header Library into ```src/Eigen``` and make sure all python packages are installed.

You can specify your simulation options in the ```run.sh``` and build the code with

```bash
./run.sh -b
```

To run the simulation simply do

```bash
./run.sh
```

If everything works correctly your output will look something like this

```
Starting program...
grid generated
file storage format: src/files/testfolder/v_n10_FV2_BC-1_1_0s_testname_step0.csv

snap nr. : delta_t : t_sim, Time: [ELAPSED < ETA]
---------------------------------------------------
0 : 0.000984895 : 0, Time: [00:00<35791394:07]
250 : 0.00109206 : 0.252848, Time: [00:00<00:00]
500 : 0.0010486 : 0.531439, Time: [00:00<00:00]
750 : 0.000893345 : 0.770378, Time: [00:00<00:00]
1000 : 0.000808048 : 0.98032, Time: [00:00<00:00]
1023 : 0.000813236 : 0.999675, Time: [00:00<00:00]
---------------------------------------------------
Total time: 0.197646
max RSS memory size: 11.9062 MB
done
```

As a next step you can visualize a snapshot using ```vis_tk.py```. Simply make a new python file in src and run. To find your snapshots look at the file storage format above.

```python
import vis_tk as v

# load snapshot
s, p, q = v.process_file("files/your_folder_name/your_filename.csv")

# plot density
v.plot_2D(s, p, q[:, 1], vmin = 0, vmax = 3)
```
For now just look through the ```vis_tk.py``` package for further plotting options and through the ```run.sh``` file for changed run options. Eventually a documentation might follow.