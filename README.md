# Introduction

This is a repository for the coursework in Computational Neuroscience in Cog-Sup 2025. It implements a HH Model, simulating a neuron model that recreate graded persistent activities.

# Usage

- Run `main.py` to execute the main script. The script will load the necessary data and run the computations. 
- For custom configurations, go to `parametres` folder to see initialization files. 
- For custom input current, go to the `main.m` to change the types variable or go to `util/create_protocol.m` to create your own protocol. There are already three setups in the `main.m` file. Uncomment the one you want to use. Default is all depolarizing current.
