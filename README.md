# Agent Based Simulation of Biological Systems 
Internship Project - Sapienza (Computer Science) 2021/2022 
---

## Abstract
The purpose of the project is to create an Agent-Based simulator for biological systems that functions by taking SBML files as input and returning executable LAMMPS code as output. This code will define a simulation based on data extrapolated from the model, and when executed, it will display the evolution of the system over time. Through this approach, it is possible to interpret even highly complex biochemical models, rich in diverse species and reactions among elements belonging to these species.

The Agent-Based modeling provides an intuitive view of the system, where each species is considered as an autonomous agent independent of other agents until they come into direct contact. The contacts between agents are the central point of the simulation, as they conceptually represent reactions between species and enable the system itself to change. The independence of individual agents, especially regarding movement within the simulation environment, allows the introduction of a level of controlled and reproducible randomness, adding realism to the simulations and making them relevant for observing the progression of a system defined by certain rules and initial conditions.

This mechanism allows for the verification or refutation of hypotheses about the final state of the system, formulated solely based on the definition of an initial state, without considering all possible intermediate states.

---

## How to run
To test the code run...
```code
python main.py
```
... which is contained inside the *prj* directory, and just follow the instruction.
