For the TChem-atm Development branch, please visit: https://github.com/PCLAeroParams/TChem-atm.

[![Github Actions Status](https://github.com/sandialabs/TChem-atm/actions/workflows/auto_test.yaml/badge.svg?branch=main)](https://github.com/sandialabs/TChem-atm/actions)
[![GitHub Pages](https://img.shields.io/badge/GitHub%20Pages-121013?logo=github&logoColor=white)](https://sandialabs.github.io/TChem-atm/)
## Installation 

The [installation](docs/installation.md) guide demonstrates how to obtain, build, and install TChem-atm along with the requisite third-party libraries.

## Citations
* [TChem-atm ](https://www.osti.gov/biblio/2472634)
```bibtex
@misc{osti_2472634,
title = {TChem-atm v1.0, Version 1.0},
author = {Safta, Cosmin and Diaz-Ibarra, Oscar and Schmidt, Michael and USDOE},
abstractNote = {SAND2024-11300O TChem-atm is a software library that was developed to solve complex kinetic models for atmospheric chemistry applications. TChem-atm interface employs a hierarchical parallelism design to exploit the massive parallelism available from modern computing platforms. It also supports gas atmospheric chemistry applications, e.g., the energy exascale earth system model. TChem can be used as a box model or coupled with a climate model to compute the time evolution of gas tracer species. Sandia National Laboratories is a multimission laboratory managed and operated by National Technology & Engineering Solutions of Sandia, LLC, a wholly owned subsidiary of Honeywell International Inc., for the U.S. Department of Energy’s National Nuclear Security Administration under contract DE-NA0003525.},
url = {https://www.osti.gov//servlets/purl/2472634},
doi = {10.11578/dc.20240925.6},
url = {https://www.osti.gov/biblio/2472634}, year = {2024},
month = {6},
note =
}
```

* [TChem: A performance portable parallel software
toolkit for complex kinetic mechanisms.](https://www.sciencedirect.com/science/article/pii/S0010465522003472)

```bibtex
@article{tchem:Kim:2022,
  title    = {{TChem: A performance portable parallel software toolkit for complex kinetic mechanisms}},
  journal  = {Computer Physics Communications},
  volume   = {285},
  pages    = {108628},
  year     = {2023},
  issn     = {0010-4655},
  author   = {Kyungjoo Kim and Oscar H. Díaz-Ibarra and Habib N. Najm and Judit Zádor and Cosmin Safta},
  keywords = {TChem, Kokkos, Performance portability, GPU, Flow chemistry}
}
```

* ["Benchmarking TChem for Potential Incorporation into E3SM as a Replacement Chemical Kinetics Solver"](sand_report/QTI_tchemV1.pdf)
```bibtex
@techreport{Diaz-Ibarra:2024:tchem,
  author      = {Diaz-Ibarra, Oscar and Schmidt, Michael J.  and Safta, Cosmin },
  title       = {{Benchmarking TChem for Potential Incorporation into E3SM as a Replacement Chemical Kinetics Solver}},
  institution = {Sandia National Laboratories},
  year        = {2024},
  number      = {SAND2024-01807R}
}
```
