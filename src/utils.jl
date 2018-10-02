"""
    cleanWrite()

Clears all timestamp directories in the current folder
"""
function cleanWrite()
    dirvec = readdir()
    dirresults = map(x->(v = tryparse(Float64,x); typeof(v) == Nothing ? 0.0 : v),dirvec)
    for i =1:length(dirresults)
        rm("$(dirresults[i])", force=true, recursive=true)
    end
    rm("*~", force=true)
end


"""
    simpleInterp(x1, x2, y1, y2, x)

Performs a linear interpolation between points (x1,y1) and (x2,y2) to
find the value of `y` at `x`.

"""
function simpleInterp(x1 ::Float64, x2 :: Float64, y1 :: Float64, y2 :: Float64, x::Float64)
    y = y1 + (y2 - y1)*(x - x1)/(x2 - x1)
    return y
end

"""
    find_tstep([kin])

kin is a definition of either Sin / Cos / Ramp kinematics

A suitable time-step is calculated based on the values of amplitude and frequency

For Sin / Cos kinematics, time step for amp*k = 0.2 is 0.015

For Eldredge ramp kinematics, time step for K = 0.2 is 0.015

For Constant kinematics, time step is 0.015

"""
function find_tstep(kin:: CosDef)
    dtstar = 15.
    dt_tmp = 0.015*0.2/(kin.k*kin.amp)
    dtstar = minimum([dt_tmp dtstar])
    return dtstar
end

function find_tstep(kin:: SinDef)
    dtstar = 15.
    dt_tmp = 0.015*0.2/(kin.k*kin.amp)
    dtstar = minimum([dt_tmp dtstar])
    return dtstar
end

function find_tstep(kin:: EldUpDef)
    dtstar = 1.
    dtstar = minimum([0.015*0.2/kin.K 0.015])
    return dtstar
end

function find_tstep(kin:: EldUptstartDef)
    dtstar = 1.
    dtstar = minimum([0.015*0.2/kin.K 0.015])
    return dtstar
end

function find_tstep(kin:: EldRampReturnDef)
    dtstar = 1.
    dtstar = minimum([0.015*0.2/kin.K 0.015])
    return dtstar
end

function find_tstep(kin :: BendingDef)
    dtstar = 15
    amp = evaluate(kin.spl, kin.spl.t[end])
    dt_tmp = 0.015*0.2/(kin.k*amp)
    dtstar = minimum([dt_tmp dtstar])
    return dtstar
end

function find_tstep(kin :: ConstDef)
    dtstar = 0.015
    return dtstar
end

# Numerical integration method: Trapezoidal
function simpleTrapz(y::Vector{T}, x::Vector{T}) where {T<:Real}
    local len = length(y)
    if (len != length(x))
        error("Vectors must be of same length")
    end
    r = 0.0
    for i in 2:len
        r += (x[i] - x[i-1]) * (y[i] + y[i-1])
    end
    r/2.0
end

# Aerofoil camber calculation from coordinate file
function camber_calc(x::Vector,airfoil::String)
    #Determine camber and camber slope on airfoil from airfoil input file

    ndiv = length(x);
    c = x[ndiv];

    cam = Array{Float64}(undef, ndiv)
    cam_slope = Array{Float64}(undef, ndiv)
    in_air = DelimitedFiles.readdlm(airfoil, Float64);
    xcoord = in_air[:,1];
    ycoord = in_air[:,2];
    ncoord = length(xcoord);
    xcoord_sum = zeros(ncoord);
    xcoord_sum[1] = 0;
    for i = 1:ncoord-1
        xcoord_sum[i+1] = xcoord_sum[i] + abs(xcoord[i+1]-xcoord[i]);
    end
    y_spl = Spline1D(xcoord_sum,ycoord);
    y_ans = Array{Float64}(undef, 2*ndiv);

    for i=1:ndiv
        y_ans[i] = evaluate(y_spl,x[i]/c);
    end
    for i=ndiv+1:2*ndiv
        y_ans[i] = evaluate(y_spl,(x[ndiv]/c) + (x[i-ndiv]/c));
    end
    cam[1:ndiv] = [(y_ans[i] + y_ans[(2*ndiv) + 1 - i])*c/2 for i = ndiv:-1:1];
    cam[1] = 0;
    cam_spl = Spline1D(x,cam);
    cam_slope[1:ndiv] = derivative(cam_spl,x);
    return cam, cam_slope
end
