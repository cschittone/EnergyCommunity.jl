![CI](https://github.com/SPSUnipi/EnergyCommunity.jl/actions/workflows/CI.yml/badge.svg)
[![Docs dev](https://img.shields.io/badge/docs-latest-blue.svg)](https://spsunipi.github.io/EnergyCommunity.jl/dev/)
[![Docs stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://spsunipi.github.io/EnergyCommunity.jl/stable/)

# EnergyCommunity.jl
Optimization of Energy Communities becomes easy with EnergyCommunity.jl!

A simple optimization of the model can be performed with


```julia
using EnergyCommunity, JuMP
using HiGHS, Plots

# create a sample Energy Community model input files in folder "data"
create_example_data("data")

# define input configuration (available in the package)
input_file = "./data/energy_community_model.yml"

# create the Energy Community model in Cooperation mode GroupCO()
ECModel = ModelEC(input_file, EnergyCommunity.GroupCO(), HiGHS.Optimizer)

# build the model
build_model!(ECModel)

# optimize the model
optimize!(ECModel)

# create some plots
plot(ECModel, output_plot_combined)
```
