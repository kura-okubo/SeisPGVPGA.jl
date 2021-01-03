module SeisPGVPGA

# module to be used over modules
using  SeisIO, FileIO, Dates, DSP, Statistics, DataFrames, CSV, Trapz, Distributed, Plots, PlotSeis

export

    map_seisPGVPGA,
    map_readPGVPGA,
    get_requeststations,
    get_starttimelist

# NOTE: Do not change the order of include.
include("func_SeisPGVPGA.jl")
include("map_SeisPGVPGA.jl")

end
