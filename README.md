## Content of the repository

This repository provides a training algorithm for Lagrangian multipliers generation for the Held-Karp TSP relaxation, and an associated branch-and-bound TSP solver.

```bash
.
├── conda_env.yml  # configuration file for the conda environment
└── solver/ 
	├── src/ # solver source code
	├── models/  # GNN models
	├── testGraphs/  # test instances
	├── loadNN.py  # GNN loading and calling functions
	├── run_tsp.sh # to run the solver
	├── run_test.sh
	├── makefile  # to compile the solver
└── training/ 
	├── src/ # solver source code
	├── trained_models/  # trained GNN models
	├── training_graphs/  # training instances
	├── trainHKgnn.py  # GNN loading and calling functions
	├── run_training.sh # to run the training
```

## Installation instructions

### 1. Importing the repository

```shell
git clone https://github.com/corail-research/learning-hk-bound.git
```

### 2. Setting up the conda virtual environment

```shell
conda env create -f environment.yml 
```

### 3. Compiling the solver

A makefile is available in the solver and trainer. First, add your python path. Then, you can compile the project as follows:

```shell
cd ./training
make
cd ../solver
make
```


## Basic use

### 1. Training a model

```shell
cd ./training
```

Edit the configuration in training/trainHKgnn, then you can start the training as follows:
```shell
cd ./training
./run_training.sh
```

### 2. Solving instances

```shell
# For TSPTW
./test.sh
```


## Technologies and tools used

* The TSP solver is implemented in C++ and is based on the [solver](https://hal.science/hal-01344070/document) of Pascal Benchimol.
* The code handling the training and the GNN inferences is implemented in Python3.
* The graph neural network architecture has been implemented in Pytorch together with DGL.
