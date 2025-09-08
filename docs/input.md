# **Input Files**

## **Gas Chemistry Input File**

The YAML input file for atmospheric chemistry consists of five sections:

* `environmental_conditions`: pressure and temperature of individual cells.
* `initial_state`: initial species concentrations within the cells.
* `reactions`: list of reactions with rate parameters and reaction type.
* `constant_species`: species that are part of the reaction mechanism but are assumed constant in time.
    * E.g., invariant species.
* `species`: list of species names.

For example, given the toy reaction $A \rightarrow B$, simulation with N cells, one reaction, and three species, we have

```yaml
NCAR-version: v1.0
environmental_conditions:
  pressure:
    evolving: false
    initial_value: [P1, ..., PN]
    units: Pa
  temperature:
    evolving: false
    initial_value: [T1, ..., TN]
    units: K
initial_state:
  A:
    initial_value: [A1, .., AN]
    units: mol m-3
  M:
    initial_value: [M1, .., MN]
    units: mol m-3
reactions:
- coefficients:
  products:
    B: 1.0
  reactants:
    A: 1.0
  type: ARRHENIUS
constant_species:
- description: tracer-CONSTANT
  name: M
species:
- description: A
  name: A
- description: B
  name: B
```

A description of the reaction types currently implemented in TChem-atm is presented in [Methodology section](methodology.md). In addition, a set of examples of input files is presented under `/src/examples/runs/atmopheric_chemistry`.

## **Gas-Aerosol Chemistry Input File**

To execute a gas-aerosol case in TChem-atm, one needs an input file with the aerosol mechanisms, `aerofile`, and a file where scenario particle information is given, `inputfile_particles`. The aerosol mechanism utilizes a YAML file and follows to the [CAMP format](https://github.com/open-atmos/camp). The scenario particle information follows this format:

```yaml
particles:
  initial_state:
    species_name_1: 
      initial_value: [1.0e-8]
    species_name_2: 
      initial_value: [1.4e-2]
  num_concentration:
   initial_value: [1.3e6]  

```
In this structure, `species_name_1` denotes the initial concentration of species 1 within a single particle. It's important to note that the initial composition for species 1 is consistent across all particles.

Additionally, to incorporate gas-phase reactions into the simulation, a gas reaction mechanism can be seamlessly integrated into the gas-aerosol case.
<!-- Future work -->
<!-- ## Gas-Particle chemistry input file. -->

<!-- ## Initial condition input file for particles. -->
