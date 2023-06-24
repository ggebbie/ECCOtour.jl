"""
    function get_filtermatrix(tin,tout)

# Arguments
- `tin`: input times
- `tout`: output times
# Output
- `Eout2in`: matrix that maps from tout to tin
- `Fin2out`: matrix that maps from tin to tout
"""
function get_filtermatrix(tin,tout)
    # function get_matrixfilter(tin,tout)

    ntin  = length(tin)
    ntout = length(tout)

    # assume evenly spaced in time
    Δtout = tout[2]-tout[1]

    # find all hi-res values bounded by tiepoint 1 and 2
    inrange = x -> ( tout[1] < x <  tout[2]) # found on plot
    iblock = findall(inrange,collect(tin))

    nblock = length(iblock)
    Eblock = zeros(nblock,2)
    for ii = 1:nblock
        Eblock[ii,1] = -(tin[iblock[ii]] - tout[2])/Δtout
        Eblock[ii,2] = (tin[iblock[ii]] - tout[1])/Δtout
    end

    # put the blocks together in a big sparse matrix.
    Eout2in = zeros(ntin+56,ntout) # make ntin too long and then trim later

    # for a real function, would want to fix this.
    ival = 3 # first two 6-hourly values off range
    jval = 1

    for ii = 1:ntout-1
        Eout2in[ival:ival+nblock-1,jval:jval+1] = Eblock
        ival += nblock
        jval += 1
    end

    # some trimming must be done.
    Eout2in = Eout2in[1:ntin,:]
    ETE = Eout2in'*Eout2in
    Fin2out = ETE\(Eout2in')

    Eout2in = convert(Array{Float32,2},Eout2in)
    Fin2out = convert(Array{Float32,2},Fin2out)

    return Eout2in, Fin2out
end

"""
    function get_filtermatrixfull(tin,tout)

    By Anthony Meza

# Arguments
- `tin`: input times
- `tout`: output times
# Output
- `Eout2in`: matrix that maps from tout to tin
- `Fin2out`: matrix that maps from tin to tout
"""
function get_filtermatrixfull(tin,tout)
    # function get_matrixfilter(tin,tout)

    ntin  = length(tin)
    ntout = length(tout)

    # assume evenly spaced in time
    Δtout = tout[2]-tout[1]

    # find all hi-res values bounded by tiepoint 1 and 2
    inrange = x -> ( tout[1] < x <  tout[2]) # found on plot
    iblock = findall(inrange,collect(tin))

    nblock = length(iblock)
    Eblock = zeros(nblock,2)
    for ii = 1:nblock
        Eblock[ii,1] = -(tin[iblock[ii]] - tout[2])/Δtout
        Eblock[ii,2] = (tin[iblock[ii]] - tout[1])/Δtout
    end

    # put the blocks together in a big sparse matrix.
    Eout2in = zeros(ntin+56,ntout) # make ntin too long and then trim later

    # for a real function, would want to fix this.
    ival = 1 # first two 6-hourly values off range
    jval = 1

    for ii = 1:ntout-1
        Eout2in[ival:ival+nblock-1,jval:jval+1] = Eblock
        ival += nblock
        jval += 1
    end

    # some trimming must be done.
    Eout2in = Eout2in[1:ntin,:]
    ETE = Eout2in'*Eout2in
    Fin2out = ETE\(Eout2in')

    Eout2in = convert(Array{Float32,2},Eout2in)
    Fin2out = convert(Array{Float32,2},Fin2out)

    return Eout2in, Fin2out
end

"""
    matrixfilter(F,froot,years,γ)
    writing it in a funny way to save computation
    issue with timeseries being read in different files
# Arguments
- `F`: filter in matrix form
- `froot`: filename root
- `years`: iterator for multiple files
- `γ`: GCM grid (meshArray type)
"""
function matrixfilter(F,froot,years,γ)

    # initialize known variables

    # use size of F to determine size of input/output
    nout = size(F,1)
    nin  = size(F,2)

    nyr = length(years)

    fname = froot*string(years[1])
    println("initialize θout")
    field = read_bin(fname,Float32,γ)

    # initialize θh
    # θout = 0.0f0 .*field[:,1:nout] # prone to have NaNs

    θout = MeshArray(γ,Float32,nout) # some nans here
    fill!(θout,0.0)
    #for tt= 1:nout
        # initialize a sub-optimal way
    #    θout[:,tt]=γ.read(tmp1,MeshArray(γ,Float32))
    #end
    println("initialization finished")

    istart = 1
    for tt = 1:nyr
        fname = froot*string(years[tt])
        println(fname)
        field = read_bin(fname,Float32,γ)

        # may need to keep track of indices
        nt = size(field,2)
        iend = istart+nt-1
        @time annualcontrib!(θout,F[:,istart:iend],field)
        istart = iend+1
    end

    return θout
end

"""
    function matrixfilter_annualcontrib!(θout,F,θin)
# Arguments
- `θout`: filtered values
- `F`: filter function as a matrix
- `θin`: input values at all times t (e.g., θ(t))
"""
function annualcontrib!(θout,F,θin)
    nt = size(θin,2)
    for tt = 1:nt
         # igood = findall(x -> abs(x)>1e-8,F[:,tt])
         # if !isempty(igood)
         #     matmul2!(θout[:,igood],F[igood,tt],θin[:,tt])
         # end
        columnscale!(θout,F[:,tt],θin[:,tt])
    end
end

"""
    function columnscale!(product,M,flux)
- `product::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}`
- `M::Array{Float32,1}`
- `flux::MeshArrays.gcmarray{Float32,1,Array{Float32,2}}`
"""
function columnscale!(product::MeshArrays.gcmarray{Float32,2,Array{Float32,2}},M::Array{Float32,1},flux::MeshArrays.gcmarray{Float32,1,Array{Float32,2}})

    nM = size(M,1)  # matrix size M
    nQ  = size(flux,1) # repeat matmul Q times

     for qq = 1:nQ # inner product over faces
         for mm = 1:nM
             if abs(M[mm]) > 1e-8
                 # igood = findall(x -> abs(x)>1e-8,F[:,tt])
                 # if !isempty(igood)
                 #     matmul2!(θout[:,igood],F[igood,tt],θin[:,tt])
                 # end
                 product[qq,mm] += flux[qq] * M[mm]
             end
        end
    end
end

"""
    function columnscale!(product,M,flux)
- `product::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}`
- `M::Array{Float32,1}`
- `flux::MeshArrays.gcmarray{Float32,2,Array{Float32,2}}`
"""
function columnscale!(product::MeshArrays.gcmarray{Float32,2,Array{Float32,2}},M::Array{Float32,2},flux::MeshArrays.gcmarray{Float32,2,Array{Float32,2}})

    nM = size(M,1)  # matrix size M
    nQ  = size(flux,1) # repeat matmul Q times
    nN  = size(M,2) # repeat matmul Q times

    for mm = 1:nM
        for qq = 1:nQ # inner product over faces
            for nn = 1:nN # inner product over faces
                product[qq,mm] += flux[qq,nn] * M[mm,nn]
            end
        end
    end
end

"""
    matrixspray(F,rmfield,frootin,frootout,years,γ)

    writing it in a funny way to save computation
    issue with timeseries being read in different files
# Arguments
- `F`: spray operator in matrix form
- `rmfield`: field to be removed from files
- `frootin`: filename root of input field
- `frootout`: filename root of output field
- `years`: iterator for multiple files
- `γ`: GCM grid (meshArray type)
"""
function matrixspray(F,rmfield,frootin,frootout,years,γ)

    # initialize known variables

    # use size of F to determine size of input/output
    nout = size(F,1)
    nin  = size(F,2)

    nyr = length(years)

    fnamein = frootin*string(years[1])
    println("initialize θout")
    field = read_bin(fnamein,Float32,γ)

    istart = 1
    nseries = []
    for tt = 1:nyr
        fnamein = frootin*string(years[tt])
        println("reading file "*fnamein)
        field = read_bin(fnamein,Float32,γ)

        # may need to keep track of indices
        push!(nseries,size(field,2))
        iend = istart+nseries[end]-1
        @time annualcontrib!(field,F[istart:iend,:],rmfield)
        istart = iend+1

        fnameout = frootout*string(years[tt])
        println("saved file "*fnameout)
        write(fnameout,field)
    end
end

"""
    function matrixsaveinterpolation(E,savefield,frootin,frootout,years,γ)

    writing it in a funny way to save computation
    issue with timeseries being read in different files. 
    Interpolates "savefield" using the E matrix onto a monthly grid. 
    Then saves the interpolated field "savefield".

    By Anthony Meza

# Arguments
- `E`: interpolating spray operator in matrix form
- `savefield`: field to be saved to multiple files
- `frootin`: filename root of input field
- `frootout`: filename root of output field
- `years`: iterator for multiple files
- `γ`: GCM grid (meshArray type)
"""
function matrixsaveinterpolation(E,savefield,frootin,frootout,years,γ)

    nyr = length(years)

    fnamein = frootin*string(years[1])
    istart = 1
    nseries = []

    for tt = 1:nyr
        fnamein = frootin*string(years[tt])
        println("reading file "*fnamein)
        field = read_bin(fnamein,Float32,γ); fill!(field,0.0)

        push!(nseries,size(field,2))

        # may need to keep track of indices
        iend = istart+nseries[end]-1
        # iend = istart+nt-1
        @time annualcontrib!(field,E[istart:iend,:],savefield)
        istart = iend+1
        fnameout = frootout*string(years[tt])
        println("saved file "*fnameout)
        write(fnameout,field)
    end
end

