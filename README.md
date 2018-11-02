# reprinz

The goal of this project is to generate lots of models of networks and neurons using `procrustes`, a powerful optimizer, and letting it vary all the parameters of the model.

# Installation

Reprinz uses the [xolotl](https://github.com/sg-s/xolotl) simulation environment.
This can be downloaded directly as a github repository or as a [MATLAB package](https://drive.google.com/uc?id=19xddT00ObfsKHaa2T1YobZeA2C2qfuIh&export=download).
Optimization is provided through the [procrustes](https://github.com/sg-s/procrustes)
package.

If you are downloading the repositories manually, you will need
* [xolotl](https://github.com/sg-s/xolotl)
* [procrustes](https://github.com/sg-s/procrustes)
* [srinivas.gs_mtools](https://github.com/sg-s/srinivas.gs_mtools)
* [cpplab](https://github.com/sg-s/cpplab)

# Use

Reprinz provides a series of scripts to generate lots of models of networks and models.
Models are specified by the kinetics, network structure, and maximal conductances.

For example, to find single-compartment bursting models using the kinetics from
[Liu _et al._ 1998](http://www.jneurosci.org/content/18/7/2309), run the script
`find_bursting_neurons_liu.m` in the `bursting-single-compartment-neurons/liu` folder.
This script will begin the optimization, automatically recruiting multiple cores
(if parallel computing is set up). The optimizer will vary the maximal conductances
to find models that produce bursting neurons like the ones shown in the paper.

Each time the optimization finishes (generating a new model), the maximal conductances are 
saved in a `.mat` file beginning with `reprinz` and ending with the computer name.
This allows the user to keep track of simulations run on different machines.

When a simulation is restarted, `reprinz` will pick up where it left off, adding
to the save file rather than overwriting it. By default, `reprinz` will stop after
10,000 models have been generated.
