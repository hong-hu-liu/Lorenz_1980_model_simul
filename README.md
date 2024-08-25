# Simulation of the Lorenz 1980 model
Matlab codes for integrating the Lorenz 80 (L80) model

The model equations are those given by Eqns (33)-(35) in

[L80] E. N. Lorenz (1980): Attractor sets and quasi-geostrophic equilibrium. J Atmos Sci 37, 1685-1699.
We refer to this model as the L80 model.

This folder contains:

-A script that defines the coefficients required to form the RHS of the L80 model  ("get_par_Lorenz9D.m")

-A script that forms the linear, nonlinear and forcing terms (RHS of the L80 model) necessary to integrate the model ("get_VF.m")  

-A script that actually integrates the L80 model using a standard RK4 scheme ("int_Lorenz9D.m")

- "run_Lorenz80_model.m" is the main script that calls the above functions to integrate the L80 model. It produces a figure showing time series of the model and a projection of the L80 model's attractor 

In "run_Lorenz80_model.m" 
par_id = 1 corresponds to the High-low Frequency (HLF) regime analyzed in [CLSM24] and [CLM21] 
par_id = 2 corresponds to the slow chaos regime analyzed in [CLSM24] and [CLM21].


# References

[CLSM24] M. D. Chekroun, H. Liu, K. Srinivasan, and J. C. McWilliams (2024): The high-frequency and rare events barriers to neural closures of atmospheric dynamics. J. Physics Complexity, 5, 025004. https://iopscience.iop.org/article/10.1088/2632-072X/ad3e59

```
@article{CLSM24,
  title={The high-frequency and rare events barriers to neural closures of atmospheric dynamics},
  author={Chekroun, M. D. and Liu, H. and Srinivasan, K. and McWilliams, J. C.},
  journal={Journal of Physics: Complexity},
  volume={5},
  pages={025004},
  year={2024},
  doi={10.1088/2632-072X/ad3e59},	
  publisher={IOP Publishing}
}
```


[CLM21] M.D. Chekroun, H. Liu, and J. McWilliams (2021): Stochastic rectification of fast oscillations on slow manifold closures. 
Proc. Natl. Acad. Sci. (PNAS), 118(48):e2113650118, 2021. 
https://doi.org/10.1073/pnas.2113650118

```
@article{CLM21,
  title={Stochastic rectification of fast oscillations on slow manifold closures},
  author={Chekroun, M. D. and Liu, H. and McWilliams, J. C.},
  journal={Proceedings of the National Academy of Sciences},
  volume={118},
  number={48},
  pages={e2113650118},
  year={2021},
  doi={10.1073/pnas.2113650118},
  publisher={National Acad Sciences}
}
```

[CLM17] M. D. Chekroun, H. Liu, and J. C. McWilliams (2017): The emergence of fast oscillations in a reduced primitive equation model and its implications for closure theories. Comput. Fluids 151, 3–22. https://doi.org/10.1016/j.compfluid.2016.07.005

```
@article{CLM17,
  title={The emergence of fast oscillations in a reduced primitive equation model and its implications for closure theories},
  author={Chekroun, M. D. and Liu, H. and McWilliams, J. C.},
  journal={Computers \& Fluids},
  volume={151},
  pages={3--22},
  year={2017},
  doi={10.1016/j.compfluid.2016.07.005},
  publisher={Elsevier}
}
```

[L80] E. N. Lorenz (1980): Attractor sets and quasi-geostrophic equilibrium. J Atmos Sci 37, 1685–1699. https://doi.org/10.1175/1520-0469(1980)037<1685:ASAQGE>2.0.CO;2

```
@article{L80,
  title={Attractor sets and quasi-geostrophic equilibrium},
  author={Lorenz, E. N.},
  journal={Journal of Atmospheric Sciences},
  volume={37},
  number={8},
  pages={1685--1699},
  year={1980},
  doi={10.1175/1520-0469(1980)037<1685:ASAQGE>2.0.CO;2}
}
```


