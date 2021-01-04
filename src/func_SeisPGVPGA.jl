"""
	Functions used in SeisPGVPGA.jl

Modified from SeisMonitoring.jl(https://github.com/kura-okubo/SeisMonitoring.jl)
2021.1.3 Kurama Okubo
"""

"""
	get_requeststations(fi_stations::String, locchan::Dict, stationlist_fopath::String)

return dictionary of request stations.

# Arguments
- 'fi_stations::String': path to the station list downloaded from IRIS GMAP2 text format.
- 'locchan::Dict': Dictionary of channel and location.

e.g. locchan=Dict("BP" => [(*, "BP*")]) requests loc=* and chan=BP* for BP network

- 'stationlist_fopath::String': Output file name for backup.

"""
function get_requeststations(fi_stations::String, locchan::Dict, stationlist_fopath::String)

	# read input file
	lines = readlines(fi_stations)

	#---make stationlist---#
	# parse algorithm:
	# 1. search until encounting #DATACENTER
	# 2. append stations to the chunc

	requeststations = Dict()
	i = 1
	while i <= length(lines)
	    if occursin("#DATACENTER", lines[i])
	        #start appending this srcchunk
	        src = string(split(split(lines[i], "=")[2], ",")[1])
	        # #manipulate IRIS variation
	        # occursin("IRIS", src) && (src = "IRIS")
	        srcchunk = DataFrame(network=String[],
	                               station=String[],
	                               location=String[],
	                               channel=String[],
	                               latitude=String[],
	                               longitude=String[],
	                               elevation=String[],
	                               sitename=String[],
	                               starttime=DateTime[],
	                               endtime=DateTime[])

	        i += 1
	        while true
	            # append request stations
	            lc = split(lines[i],"|")
	            net=lc[1]
	            sta=lc[2]
	            lat=lc[3]
	            lon=lc[4]
	            ele=lc[5]
	            sitename=lc[6]
	            starttime=DateTime(lc[7])
	            endtime=DateTime(lc[8])

	            # choose request channel if locchan dict contains the network; otherwise request all (e.g. *, *)
	            if net ∈ keys(locchan)
	                reqs = locchan[net]
	                for req in reqs
	                    loc, chan = string.(req)

	                    dftemp = DataFrame(network=net, station=sta,location=loc,channel=chan,
	                    latitude=lat, longitude=lon,elevation=ele, sitename=sitename, starttime=starttime,endtime=endtime)
	                    append!(srcchunk,dftemp)
	                end

	            else
	                dftemp = DataFrame(network=net, station=sta,location="*",channel="*",
	                latitude=lat, longitude=lon,elevation=ele, sitename=sitename, starttime=starttime,endtime=endtime)
	                append!(srcchunk,dftemp)
	            end

	            i += 1
	            if (i > length(lines) || isempty(lines[i]))
	                # push srcchunk to master Dictionary
	                requeststations[src] = srcchunk

	                # output csv
	                CSV.write(stationlist_fopath[1:end-5]*"_$(src).csv", srcchunk)
	                # next src chunk or end of file
	                break
	            end
	        end

	    else
	        i += 1
	        continue
	    end
	end
	return requeststations
end


"""

get_starttimelist(st::DateTime, et::DateTime, unittime::Float64)
calculate start time list for parallel downloading

    st: start time
    et: end time
    unittime: unit time in Second

    this function returns
    stlist: list of start time

    e.g.
    st = DateTime(2019,1,1,0,0,0)
    et = DateTime(2019,1,1,12,0,0)
    unittime = 3600

    stlist = get_starttimelist(st, et, unittime)

"""
function get_starttimelist(st::DateTime, et::DateTime, unittime::Real)

    reftime = st
    stlist = []

    while reftime < et
        push!(stlist, string(reftime))
        reftime += Dates.Second(float(unittime))
    end

    return stlist
end


"""
	get_requeststr(df::DataFrame, numstationperrequest::Int))

return request str following web_chanspec of SeisIO.get_data.

- 'numstationperrequest::Int': number of stations per one request to avoid too many request statinos at one time.
"""
function get_requeststr(df::DataFrame, numstationperrequest::Int)

	reqstrs = Array{Array{String,1},1}(undef, 0) # array of request stations per HTTP request
	for inds = Iterators.partition(1:size(df)[1], numstationperrequest)
		reqstr  = String[]
		for i in inds
			rst = join([df.network[i], df.station[i], df.location[i], df.channel[i]], ".")
			push!(reqstr, rst)
		end
		push!(reqstrs, reqstr)
	end
	#
	# for i = 1:size(df)[1]
	# 	rst = join([df.network[i], df.station[i], df.location[i], df.channel[i]], ".")
	# end
	return reqstrs
end


"""
    check_and_get_data(ex::Expr, requeststr::String)

Download seismic data, removing instrumental response and saving into JLD2 file.

# Arguments
- `ex::Expr`        : expression of get data includin all request information

# Output
- `S::SeisData`     : downloaded SeisData
"""
function check_and_get_data(ex::Expr)

	# we try the same download request upto 3 times because it sometimes fails
	# due to network error.
	download_itr = 1

	S = SeisData()

	while download_itr <= 3
		S = try
			#remove comment out below if you want to print contents of get_data()
			# println(ex)
			eval(ex);
		catch
		end

		if !isnothing(S) && !isempty(S)
			break;
		else
			download_itr += 1
			#println("retry downloading $(download_itr).")
		end
	end

	if isnothing(S) || isempty(S)
		return 1
	end

	for j = 1:S.n
		!isnothing(S[j]) && !isempty(S[j]) && !isempty(S[j].t) ? S.misc[j]["dlerror"] = 0 : S.misc[j]["dlerror"] = 1
	end

	return S
end


"""
    norm2(A; dims)

Compute norm along dimension of array.
https://github.com/JuliaLang/julia/issues/34830
"""
function norm2(A; dims)
    B = sum(x -> x^2, A; dims=dims)
    B .= sqrt.(B)
end


"""
    df_rms(P_vel::DSP.Periodograms.Periodogram,  fmin::Real, fmax::Real)

compute rms displacement from velocity power spectrum density.
Modified from SeismoRMS (https://github.com/ThomasLecocq/SeismoRMS) written by Thomas Lecocq.
https://github.com/ThomasLecocq/SeismoRMS/blob/master/seismosocialdistancing.py l.256
"""
function df_rms(P_vel::DSP.Periodograms.Periodogram, fmin::Real, fmax::Real)

    ix = findall( x-> (fmin<x<fmax), P_vel.freq)
    spec = P_vel.power[ix]
    f = P_vel.freq[ix]
    w2f = (2.0 * pi .* f)
    # NOTE: welch_pgram() returns power spectrum density, not in dB. so we don't compute 10^(S/10)
    vamp = spec
    damp = vamp ./ (w2f.^2)
    d_rms = sqrt(trapz(f, damp))

    return d_rms
end

"""
    compute_pgvpga_drms(S::SeisData, bp_fmin::Real, bp_fmax::Real, drms_fmin::Real, drms_fmax::Real)

compute pgv, pga and rms displacement from 3 components SeisData.

# Arguments
- `S::SeisData`     : 3 components (i.e. S.n==3) SeisData.
- `bp_fmin::Real`   : minimum frequency for preprocessing.
- `bp_fmax::Real`   : maximum frequency for preprocessing.
- `drms_fmin::Real`   : minimum frequency for d_rms integration
- `drms_fmax::Real`   : maximum frequency for d_rms integration

# Output
- `pgv::Float64`     : Peak ground velocity magnitude (norm2(vx, vy, vz))
- `pga::Float64`     : Peak ground acceleration magnitude (norm2(ax, ay, az))
- `d_rms_mean::Float64`     : mean strength of noise displacement between drms_fmin and drms_fmax
"""
function compute_pgvpga_drms(S::SeisData, bp_fmin::Real, bp_fmax::Real, drms_fmin::Real, drms_fmax::Real)

    if S.n != 3
        @warn "$(S.id) does not have 3 components. skip."
        return (NaN, NaN, NaN)
    end

	#assert data length
	if !(length(S[1].x[:]) == length(S[2].x[:]) == length(S[3].x[:]))
		@warn "$(S.id) data length is different between components. skip."
		return (NaN, NaN, NaN)
	end

    #---Preprocessing---#
    detrend!(S)
    demean!(S)
    taper!(S, t_max=180, α=0.1)
    filtfilt!(S, fl=bp_fmin, fh=bp_fmax, rt="Bandpass", dm="Butterworth")

    #---Convert trace---#
    S_acc = convert_seis(S, units_out="m/s2", v=0);
    S_vel = convert_seis(S, units_out="m/s", v=0);

    #---Preprocessing again for velocity and acceleration---#
    detrend!(S_acc)
    demean!(S_acc)
    taper!(S_acc, t_max=180, α=0.1)
    filtfilt!(S_acc, fl=bp_fmin, fh=bp_fmax, rt="Bandpass", dm="Butterworth")

    detrend!(S_vel)
    demean!(S_vel)
    taper!(S_vel, t_max=180, α=0.1)
    filtfilt!(S_vel, fl=bp_fmin, fh=bp_fmax, rt="Bandpass", dm="Butterworth")

    #---Compute acceleration and velocity magnitude---#
    acc_arr=zeros(size(S_acc[1].x, 1), 3)
    vel_arr=zeros(size(S_vel[1].x, 1), 3)

    for i = 1:3
        acc_arr[:,i]=S_acc[i].x[:]
        vel_arr[:,i]=S_vel[i].x[:]
    end

    S_accmag = deepcopy(S_acc[1])
    S_accmag.id = S_accmag.id[1:end-3]*"MAG"
    S_accmag.x = norm2(acc_arr, dims=2)[:, 1]

    S_velmag = deepcopy(S_vel[1])
    S_velmag.id = S_velmag.id[1:end-3]*"MAG"
    S_velmag.x = norm2(vel_arr, dims=2)[:, 1];

    #---Pick up PGV and PGA---#
    pgv = maximum(S_velmag.x)
    pga = maximum(S_accmag.x)

    #---Compute d_rms---#

    d_rms_3comp=zeros(3)
    for i = 1:3
        N = length(S_vel[i].x)
        welch_n = div(N, 8)
        welch_noverlap = div(welch_n, 2)
        P_vel = welch_pgram(S_vel[i].x, welch_n, welch_noverlap, onesided=true, nfft=nextfastfft(welch_n), fs=S_acc[1].fs, window=DSP.Windows.hanning)
        P_vel.power[P_vel.power.<=1e-20] .= 1e-20 # avoid zero in power
        d_rms_3comp[i] = df_rms(P_vel, drms_fmin, drms_fmax)
    end

    # take mean rms of 3 components
    d_rms_mean = mean(d_rms_3comp)

    return (pgv, pga, d_rms_mean)

end
