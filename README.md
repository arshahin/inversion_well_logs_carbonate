# inversion_well_logs_carbonate
This is a package for the inversion of well logs associated with carbonates with two pore systems.

Very Fast Simulated Annealing (VFSA) is the global and stochastic optimization algorithm and here we customize it for well log inversion.  

To run the codes, please run "vfsa_density_sonic_resisit_synthetic_sw_equal_v07.m" code. It will call all the depencencies for the project.

Please cite these papers if you use these codes.

https://www.researchgate.net/publication/356863718_Multi-physics_rock_templates_to_seismically_characterise_complex_carbonates
Multi-physics rock templates to seismically characterise complex carbonates
December 2021Exploration Geophysics
DOI:10.1080/08123985.2021.2010500



https://www.researchgate.net/publication/361250177_MULTI-SCALE_INVERSION_OF_SUBSURFACE_DATA_AIMED_AT_CHARACTERIZING_HETEROGENEOUS_CARBONATE_RESERVOIRS



Contrubutions are welcome to update or convert to other programming languages.


#####################################################################################################################
Single-physics rock templates have been widely addressed in the literature especially for sandstone reservoirs. Nevertheless, multi-physics rock templates (MPRT) have not been broadly studied to characterise complex carbonates. Multi-physics measurements lead to generate MPRTs which provide visual aid to petrophysicists and seismic rock physicists to classify facies and determine reservoir rock and fluids. Our research is oriented around two sequential stages. In the first stage, we make three independent porosity measurements (Archimedes, µCT and NMR) on core carbonate plugs from northern Niagaran reef. Resistivity, P&S-wave ultrasonic measurement and joint modelling of the same brine saturated plugs help us to fine-tune the model parameters through a global optimisation algorithm. Optimisation algorithm provides vuggy and micro-porosities close to independently measured porosities using NMR and µCT. In the second stage, we extend the technique from core data to well logs. We integrate mass balance equations to model bulk density and staged differential effective medium (SDEM) theory to model elastic and electrical resistivity of dual-porosity carbonates. Then, we design a stochastic global algorithm to simultaneously invert petrophysical properties. The inversion algorithm iteratively recovers the petrophysical properties including intergranular porosity, vuggy porosity, water saturation, salinity and matrix properties. Critical porosity, resistivity lithology exponents and sonic length scales for different pore systems are also estimated with meaningful accuracy.
