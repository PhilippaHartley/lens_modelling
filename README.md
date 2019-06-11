
# Lens modelling with MULTINEST

Lens modelling is performed using the PyMultiNest package (http://www.aanda.org/articles/aa/abs/2014/04/aa22971-13/aa22971-13.html). Priors must be set in prior(), using a cube with a length equal to the number of parameters. The ordering of the cube will match the order of parameters passed into pelim() as spar and gpar. These parameters are given names (snames and gnames) in pelim(); change these names accordingly to make sure that the output files contain the correct parameter names. 



### Usage

```
python peluv_multinest.py

```

### Using with MPI 

Install the mpi4py python library via package manage, pip or conda and run

```
mpiexec -n N python peluv_multinest.py
```

where N is the number of cores you would like to use.


### Model parameter uncertainties

After the modelling has finished, marginalised PDF plots and associated uncertainties are obtained by running 

```
python multinest_marginals.py
```

```multinest_marginals2``` has been written for the purpose of making plot fonts easier to read for publication-quality figures. The settings in this script could be tweaked according to the number of parameters used, to produce nice-looking plots.
