export hanncoeffs, hannsum, hannsum!, hannfilter

"""
    hanncoeffs(tin,tout,L)
# Arguments
- `tin`: input dataset time (vector valued)
- `tout`: time where filtered output will be created
- `L`: filter length 
# Output
- `h`: normalized coefficients of Hann(ing) filter
"""
function hanncoeffs(tin,tout,L)
    ctin = collect(tin)
    δ = abs.(ctin .- tout)./(2.0.*L)
    h = 0 .*ctin
    #ih = findall(δ.< 0.5)
    for hh = 1:length(h)
        if δ[hh] < 0.5
            h[hh] = (cos.(π.*δ[hh])).^2
        end
    end
    h ./= sum(h)  # normalize
    return h
end

"""
    hannsum(θin,tin,tout,L)
    use this when a timeseries is available
# Arguments
- `θ`: input value at all times t (e.g., θ(t))
- `tout`: time where filtered output will be created
- `tin`: input dataset time (vector valued)
- `L`: filter length 
# Output
- `θout`: output value at time th (e.g., θh(th))
"""
function hannsum(θin,tin,tout,L)
    # function hannsum(θin,tin,tout,L)

    h = hanncoeffs(tin,tout,L)
    
    # iterate over all filter coeffs
    i = 1
    nh = length(tin)
    θout = 0.0
    for tt = 1:nh
        θout += h[tt]*θin[tt]
    end
    return θout
end

"""
    hannsum!(θout,θin,tin,tout,L)
    use this when a timeseries is available
# Arguments
- `θout`: output value at time th (e.g., θh(th))
- `θin`: input value at all times t (e.g., θ(t))
- `tin`: input dataset time (vector valued)
- `tout`: time where filtered output will be created
- `L`: filter length 
"""
function hannsum!(θout,θin,tin,tout,L)
    # function hannsum(θin,tin,tout,L)
    # use this when a timeseries is available
    #
    # θ  = input value at all times t (e.g., θ(t))
    # th = time where filtered output will be created
    # t  = input dataset time (vector valued)
    # L = filter length 
    #
    # θh = output value at time th (e.g., θh(th))

    h = hanncoeffs(tin,tout,L)
    
    # iterate over all filter coeffs
    i = 1
    nh = length(tin)
    #θout = 0.0
    for tt = 1:nh
        θout += h[tt]*θin[tt]
    end
end

"""
    function hannsum!(θout,θin,tin,tout,toutindex,L::Float64)
    
     Call this version for mesharrays
# Arguments
- `θout::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}`: output value at time tout (e.g., θout(tout))
- `θin::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}`: input value at all times t (e.g., θ(t))
- `tin::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}`: input dataset time (vector valued)
- `tout::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}`: time where filtered output will be created
- `toutindex::Integer`: output times?
- `L::Float64`: filter length
"""
function hannsum!(θout::MeshArrays.gcmarray{Float32,2,Array{Float32,2}},θin::MeshArrays.gcmarray{Float32,2,Array{Float32,2}},tin::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},tout::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},toutindex::Integer,L::Float64)

    h = hanncoeffs(tin,tout[toutindex],L)
    
    # iterate over all filter coeffs
    i = 1
    nin = length(tin)
    for tt = 1:nin
        if h[tt] > 0
            θout[:,toutindex] += h[tt]*θin[:,tt]
        end
    end
end

"""
    function hannfilter(θin,tin,tout,L)
# Arguments
- `θin::Array{Float32,1}`: input
- `tin`
- `tout`
- `L`
# Output
- `θout`: output
"""
function hannfilter(θin::Array{Float32,1},tin,tout,L)
    # assumes all values available in timeseries
    # initialize θh

    nt = length(tin)
    nh = length(tout)
    θout = zeros(Float32,nh)
    println("calling simple hanningfilter")
    for tt= 1:nh
        hannsum!(θout[tt],θin,tin,tout[tt],L)
    end
    return θout
end

"""
    function hannfilter(θin,tin,tout,L,γ)
# Arguments
- `θin::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}`
- `tin::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}`
- `tout::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}`
- `L::Float64`
- `γ::gcmgrid`
# Output
- `θout`
"""
function hannfilter(θin::MeshArrays.gcmarray{Float32,2,Array{Float32,2}},tin::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},tout::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}},L::Float64,γ::gcmgrid)
    # use this for gcmarray
    # initialize θout

    nin = length(tin)
    nout = length(tout)
    #println(nin,nout)
    #θout = zeros(Float32,nh)
    θout = MeshArray(γ,Float32,nout) # some nans here
    tmp1=zeros(Float32,Tuple(γ.ioSize))

    for tt= 1:nout
        println(tt)
        # initialize on the fly.
        θout[:,tt]=γ.read(tmp1,MeshArray(γ,Float32))
        #println("call hanningsum!")
        hannsum!(θout,θin,tin,tout,tt,L)
    end
    return θout
end
