using SeisIO, FileIO, Dates, DSP, Statistics, DataFrames, CSV, Trapz, Distributed, Plots, PlotSeis
@everywhere using SeisPGVPGA


#---output project directory (all outputs are dumped into here)---#
fodir = "output"

#---parameters for request stations---#
fi_stations = "BP_gmap-stations.txt" #"BP_gmap-stations.txt"
locchan = Dict(
            "BP" => [(*, "BP*"), (*, "SP*")]
            )

#---parameters for processsing---#
starttime = DateTime(2004, 10, 1)
endtime  = DateTime(2004, 10, 2)
download_time_unit = 86400.0 #[s]
sampling_frequency = 20.0 #[Hz]
get_data_opt = [true, true, true, true, true]#: [unscale, demean, detrend, taper, ungap]
download_margin = 300.0 # [s] download margin
save_seisdata = true # True if save seisdata, otherwise save only pgvpga output.
bp_fmin = 0.1;
bp_fmax = 5.0;
drms_fmin = 0.9;
drms_fmax = 2.0;

#---make output directories---#
InputDict=Dict()
mkpath(fodir)
InputDict["pgvpga_dir"] = joinpath(fodir, "pgvpga")
InputDict["seisdata_dir"] = joinpath(fodir, "seisdata")
mkpath(InputDict["pgvpga_dir"])
mkpath(InputDict["seisdata_dir"])

#---make stationlist---#
stationlist_fopath = joinpath(fodir, "request_stations")
InputDict["requeststations"] = get_requeststations(fi_stations, locchan, stationlist_fopath);

starttimelist = get_starttimelist(starttime, endtime, download_time_unit)
InputDict["starttimelist"] = starttimelist
InputDict["download_time_unit"] = download_time_unit
InputDict["fodir"] = fodir
InputDict["sampling_frequency"] = sampling_frequency
InputDict["get_data_opt"] = get_data_opt
InputDict["download_margin"] = download_margin
InputDict["save_seisdata"] = save_seisdata

if (drms_fmin < bp_fmin) || (bp_fmax<drms_fmax) || (drms_fmax < drms_fmin)
    error("The order of frequency band is incorrect. Please order like bp_fmin < drms_fmin< drms_fmax< bp_fmax.")
end

InputDict["bp_fmin"] = bp_fmin
InputDict["bp_fmax"] = bp_fmax
InputDict["drms_fmin"] = drms_fmin
InputDict["drms_fmax"] = drms_fmax

#----Restrict number of processors------#
#NEVER CHANGE THIS THRESHOLD OTHERWISE IT OVERLOADS THE DATA SERVER
np = nprocs()
if np > 64 throw(DomainError(np, "np must be smaller than 64.")) end
#---------------------------------------#

# Start downloading data
println("-------START DOWNLOADING AND PROCESSING PGVPGA--------")
t_download = 0
t_download = @elapsed pmap(x -> map_seisPGVPGA(x, InputDict), 1:length(starttimelist))

println("-------START CONVERTING OUTPUT INTO CSV--------")

files_all = String[]

for (root, dirs, files) in walkdir(InputDict["pgvpga_dir"])
   for file in files
       fi = joinpath(root, file)
        (split(fi, ".")[end] == "dat") && push!(files_all, fi)# filter if it is .seisio
   end
end
t_convert = 0
t_convert = @elapsed df_mapped = pmap(x -> map_readPGVPGA(x), files_all)

df_all = DataFrame(id=String[], starttime=DateTime[], endtime=DateTime[], midtime=DateTime[], PGV=Float64[],
                    PGA=Float64[], drms_mean=Float64[])

for df in df_mapped
    (!isnothing(df) && !isnan(df.PGV[1])) && append!(df_all, df)
end

csvname = joinpath(InputDict["fodir"], "pgv_pga_all.csv")
CSV.write(csvname, df_all)

printstyled("---Summary---\n"; color = :cyan, bold = true)
println("Total time for download and process data = $(t_download)[s]")
println("Total time for gather data = $(t_convert)[s]")
printstyled("SeisPGVPGA has been successfully done.\n"; color = :cyan, bold = true)
