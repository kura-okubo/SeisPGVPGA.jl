"""
	map function for parallelized processing of PGV and PGA

2021.1.4 Kurama Okubo
"""

using CSV, DataFrames

"""
	map_seisPGVPGA(startid, InputDict::Dict)

Download 3 components seismogram and compute PGV, PGA and noise strength d_rms.
"""
function map_seisPGVPGA(startid, InputDict::Dict)

	download_time_unit = InputDict["download_time_unit"]
	fodir              = InputDict["fodir"]
	starttimelist = InputDict["starttimelist"]

	#SeisIO getdata option
	!haskey(InputDict, "get_data_opt") && (InputDict["get_data_opt"] = [true, true, true, true, true])#: [unscale, demean, detrend, taper, ungap]
	!haskey(InputDict, "sampling_frequency") && (InputDict["sampling_frequency"] = false) # default margin: 5 minutes
	!haskey(InputDict, "download_margin") && (InputDict["download_margin"] = 5 * 60) # default margin: 5 minutes

	#show progress
	if starttimelist[startid][end-8:end] == "T00:00:00"
	    println("start downloading $(starttimelist[startid])")
	end

	request_src_chanks = InputDict["requeststations"]

	for src in keys(request_src_chanks)
	    # replace source name
	    src=="IRISDMC" && (reqsrc="IRIS")
	    
	    #---download data---#
	    requeststrs = get_requeststr(request_src_chanks[src], 1) #  request single station by station

	    for requeststr in requeststrs

    		#---Skip if this request station has been already performed due to restart processing etc.---#
    		s_str = string(DateTime(starttimelist[startid]))[1:19]
	        # compute end time:
	        et = DateTime(starttimelist[startid]) + Second(download_time_unit)
	        e_str=string(et)[1:19]
	        # split network, station, location and channel
	        net, sta, loc, cha = split(requeststr[1], ".")
	        varname = join([net, sta], ".")*"__"*s_str*"__"*e_str

	        # make hierarchical directory
	        dir1 = joinpath(InputDict["pgvpga_dir"], join([net, sta], ".")) #ex. BP.LCCB..BP1
	        ts_year = split(split(varname, "__")[2], "-")[1]# start year
	        dir2 = joinpath(dir1, ts_year) #ex. 2014
	        mkpath(dir2)
	        # dump file
	        if ispath(joinpath(dir2, varname*".dat"))
	        	println("skip $(varname*".dat").")
	        	continue;
	        end

	        #NOTE: method is currently fixed "FDSN", which downloads waveforms.
	        method = "FDSN"

	        # including download margin
	        starttime = string(DateTime(starttimelist[startid]) - Second(InputDict["download_margin"]))
	        dltime = download_time_unit + 2 * InputDict["download_margin"]

	        ex = :(get_data($(method), $(requeststr), s=$(starttime), t=$(dltime),
	         v=$(0), src=$(reqsrc), unscale=$(InputDict["get_data_opt"][1]),
	          demean=$(InputDict["get_data_opt"][2]), detrend=$(InputDict["get_data_opt"][3]),taper=$(InputDict["get_data_opt"][4]),
	          ungap=$(InputDict["get_data_opt"][5]), rr=true))
	        
	        t_dl = @elapsed Stemp = check_and_get_data(ex)
	        # If download is failed, skip
	        Stemp == 1 && continue
	        # If Stemp does not have 3 components, skip
	        Stemp.n != 3 && continue
	        
	        # manipulate download_margin
	        # NOTE: using SeisIO.sync()
	        stsync = DateTime(starttimelist[startid])
	        etsync = u2d(d2u(stsync) + download_time_unit)
	        SeisIO.sync!(Stemp, s=stsync, t=etsync, v=0)

	        #---downsampling---#
	        if InputDict["sampling_frequency"] isa Number
	            for j = 1:Stemp.n
	                if Stemp.misc[j]["dlerror"] == 0 && !isempty(Stemp[j].t)
	                    # downsample
	                    if Stemp.fs[j] > InputDict["sampling_frequency"]
	                        SeisIO.resample!(Stemp, chans=j, fs=float(InputDict["sampling_frequency"]))
	                    end
	                end
	            end
	        end
	        
	        #---compute pgv, pga and d_rms---#
	        pgv, pga, drms_mean = compute_pgvpga_drms(Stemp, InputDict["bp_fmin"], InputDict["bp_fmax"],
	                                                    InputDict["drms_fmin"], InputDict["drms_fmax"])
	        
	        #---output result into dat file---#
	        # dump file
	        open(joinpath(dir2, varname*".dat"), "w") do io
	            write(io, "id,starttime,endtime,PGV,PGA,drms_mean\n")
	            write(io, "$(join([net, sta], ".")), $(s_str), $(e_str), $(pgv), $(pga), $(drms_mean)")
	        end;        
	        
	        #---save seisdata---#
	        if InputDict["save_seisdata"] 
	            # make hierarchical directory
	            dir1 = joinpath(InputDict["seisdata_dir"], join([net, sta], ".")) #ex. BP.LCCB..BP1
	            ts_year = split(split(varname, "__")[2], "-")[1]# start year
	            dir2 = joinpath(dir1, ts_year) #ex. 2014
	            mkpath(dir2)
	            # dump file
	            wseis(joinpath(dir2, varname*".seisio"), Stemp)
	        end
	    end
	end

end


"""
    map_readPGVPGA(fi_csv::String)

read pgv, pga and d_rms from dat file.
"""
function map_readPGVPGA(fi_csv::String)
    df = CSV.File(fi_csv,  header=1, delim=',') |> DataFrame
    df["midtime"] = u2d((d2u(df["starttime"][1]) + d2u(df["endtime"][1]))/2)
    return df
end
