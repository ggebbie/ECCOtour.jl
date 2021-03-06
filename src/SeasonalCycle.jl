"""
    function seasonal_matrices(fcycle,t,overtones=1)
# Arguments
- `fcycle`: frequency of seasonal cycle
- `t`: time
- `overtones=1`: optional argument for number of overtones
# Output
- `E`: matrix that solves E*parameters= seasonal cycle
- `F=E†`: generalized inverse of `E`
"""
function seasonal_matrices(fcycle,t,overtones=1)

    ω = 2π * fcycle # day^{-1} , seasonal frequency
    nt = length(t)
    freqs = overtones+1
        
    # times 2 for sines and cosines 
    nβ = freqs*2 + 1 # plus one to remove mean

    # would be neat to return whatever type goes in
    E = zeros(Float32,nt,nβ)
    println(size(E))
    
    # mean value: make it an option
    E[:,1] = ones(Float32,nt,1)
    println(size(E))

    for i = 1:freqs # change to 3 or 4 to get higher harmonics
        println(length(sin.(i*ω*t)))
        E[:,i+1] = sin.(i*ω*t)
        E[:,i+freqs+1] = cos.(i*ω*t)
    end

    F = (E'*E)\E' # least squares estimator

    return E,F
end

"""
    function trend_matrices(t)
# Arguments
- `t`: time
# Output
- `E`: matrix that solves E*parameters= timeseries
- `F=E†`: generalized inverse of `E`
"""
function trend_matrices(t)

    nt = length(t)
    tbar = StatsBase.mean(t)

    # 2 = y-intercept and trend
    nβ = 2

    # would be neat to return whatever type goes in
    E = zeros(Float32,nt,nβ)

    # parameter 1: value at mean t
    # mean value: make it an option
    E[:,1] = ones(Float32,nt,1)

    E[:,2] = t .- tbar

    F = (E'*E)\E' # least squares estimator

    return E,F
end


