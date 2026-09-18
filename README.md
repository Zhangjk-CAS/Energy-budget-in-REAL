# Code Description

This archive contains the calculation and plotting programs associated with the figures in *Lateral Open Boundaries in the Energy Budget of a Regional Ocean Model for the Western Pacific and Adjacent Seas*.

## Organization and use

Each subdirectory contains the calculation programs and plotting programs for the corresponding diagnostics or figures. Run the calculation program to produce the required diagnostic data, then run the plotting program to generate the figure.

A single program may be used for multiple subplots. To generate a specific subplot, change the relevant key fields in the program, such as the selected experiment, variable, or region, and rerun it. Check the input and output paths before running the programs in a new environment.

## `model and analysis scripts/real/`

- `real/Code/` contains the complete LICOM-REAL source code and the model configuration used for the simulations.
- `real/fort22toNC/` contains the utility used to convert the relevant instantaneous model output into NetCDF format.
