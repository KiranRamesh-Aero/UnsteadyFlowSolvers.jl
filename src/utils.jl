"""
    cleanWrite()

Clears all timestamp directories in the current folder
"""
function cleanWrite()

    dirvec = readdir()
    if "Step Files" in dirvec
        try
            rm("Step Files", recursive=true)
        catch
            println(" ERROR: Unable to reset 'Step Files' directory")
        end
    else
        dirvec = readdir()
        dirresults = map(x->(v = tryparse(Float64,x); typeof(v) == Nothing ? 0.0 : v),dirvec)
        for i =1:length(dirresults)
            rm("$(dirresults[i])", force=true, recursive=true)
        end
        # rm("*~", force=true)
    end
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

# added EldRampReturntstartDef by Joe
function find_tstep(kin:: EldRampReturntstartDef)
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
    if len != length(x)
        error("Vectors must be of same length")
    end
    r = 0.0
    for i in 2:len
        r += (x[i] - x[i-1]) * (y[i] + y[i-1])
    end
    r/2.0
end

# Aerofoil camber calculation from coordinate file
function camber_calc(x::Vector,airfoil::Union{String,SubString{String}})
    #Determine camber and camber slope on airfoil from airfoil input file

    ndiv = length(x);
    c = x[ndiv];

    cam = zeros(ndiv)
    cam_slope = zeros(ndiv)
    in_air = DelimitedFiles.readdlm(airfoil, Float64);
    xcoord = in_air[:,1];
    ycoord = in_air[:,2];
    ncoord = length(xcoord);
    xcoord_sum = zeros(ncoord);

    for i = 1:ncoord-1
        xcoord_sum[i+1] = xcoord_sum[i] + abs(xcoord[i+1]-xcoord[i]);
    end
    y_spl = Spline1D(xcoord_sum,ycoord);

    y_ans = zeros(2*ndiv);

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

function camber_thick_calc(x::Vector,coord_file::String)
    #Determine camber and camber slope on airfoil from airfoil input file

    ndiv = length(x);
    c = x[ndiv];

    cam = zeros(ndiv)
    cam_slope = zeros(ndiv)
    thick = zeros(ndiv)
    thick_slope = zeros(ndiv)

    if coord_file[1:6] == "NACA00"
        m = parse(Int, coord_file[5])/100.
        p = parse(Int, coord_file[6])/10.
        th = parse(Int, coord_file[7:8])/100.

        b1 = 0.2969; b2 = -0.1260; b3 = -0.3516; b4 = 0.2843; b5 = -0.1015
        for i = 2:ndiv
            thick[i] = 5*th*(b1*sqrt(x[i]) + b2*x[i] + b3*x[i]^2 + b4*x[i]^3 + b5*x[i]^4)
            thick_slope[i] = 5*th*(b1/(2*sqrt(x[i])) + b2 + 2*b3*x[i] + 3*b4*x[i]^2 + 4*b5*x[i]^3)
        end
        thick[1] = 5*th*(b1*sqrt(x[1]) + b2*x[1] + b3*x[1]^2 + b4*x[1]^3 + b5*x[1]^4)
        thick_slope[1] = 2*thick_slope[2] - thick_slope[3]
        rho = 1.1019*th*th*c
        cam[1:ndiv] .= 0.
        cam_slope[1:ndiv] .= 0.

    elseif coord_file[1:8] == "Cylinder"
        for i = 1:ndiv
            theta = acos(1. - 2*x[i]/c)
            thick[i] = 0.5*c*sin(theta)
        end
        thick_spl = Spline1D(x,thick)
        thick_slope[1:ndiv] = derivative(thick_spl,x)
        cam[1:ndiv] .= 0.
        cam_slope[1:ndiv] .= 0.
        rho = 0.5*c
    elseif coord_file[1:9] == "FlatPlate"
        th = parse(Int, coord_file[10:13])/10000.
        r = th*c/2
        for i = 2:ndiv-1
            if x[i] <= r
                thick[i] = sqrt(r^2 - (x[i] - r)^2)
                thick_slope[i] = -(x[i] - r)/(sqrt(r^2 - (x[i] - r)^2))
            elseif  x[i] >= c-r
                thick[i] = sqrt(r^2 - (x[i] - c + r)^2)
                thick_slope[i] = -(x[i] - c + r)/sqrt(r^2 - (x[i] - c + r)^2)
            else
                thick[i] = r
            end
        end
        thick[1] = sqrt(r^2 - (x[1] - r)^2)
        thick_slope[1] = 2*thick_slope[2] - thick_slope[3]
        thick[ndiv] = sqrt(r^2 - (x[ndiv] - c + r)^2)
        thick_slope[ndiv] = 2*thick_slope[ndiv-1] - thick_slope[ndiv-2]

        rho = r
        for i = 1:ndiv
            cam[i] = 0.
            cam_slope[i] = 0.
        end
    else
        coord = readdlm(coord_file)
        ncoord = length(coord[:,1])
        if (0. in coord[:,1]) == false
            error("Airfoil file must contain leading edge coordinate (0,0)")
        else
            nle = find(x->x==0, coord[:,1])[1]
            if coord[nle,2] != 0.
                error("Airfoil leading edge must be at (0,0)")
            end
        end

        zu_spl = Spline1D(reverse(coord[1:nle,1]), reverse(coord[1:nle,2]),k=1)
        zl_spl = Spline1D(coord[nle:ncoord,1], coord[nle:ncoord,2],k=1)

        zu = Array{Float64}(ndiv)
        zl = Array{Float64}(ndiv)

        for i=1:ndiv
            zu[i] = evaluate(zu_spl,x[i]/c)
            zl[i] = evaluate(zl_spl,x[i]/c)
        end

        cam[1:ndiv] = [(zu[i] + zl[i])*c/2 for i = 1:ndiv]
        thick[1:ndiv] = [(zu[i] - zl[i])*c/2 for i = 1:ndiv]
        cam_spl = Spline1D(x,cam)
        thick_spl = Spline1D(x,thick)
        cam_slope[1:ndiv] = derivative(cam_spl,x)
        thick_slope[1:ndiv] = derivative(thick_spl,x)
        rho = readdlm("rho")[1]
    end
    return thick, thick_slope,rho, cam, cam_slope
end
