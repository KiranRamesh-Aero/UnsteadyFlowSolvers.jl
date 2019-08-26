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


    cam = zeros(ndiv)
    cam_slope = zeros(ndiv)
    in_air = readdlm(airfoil, Float64);
    xcoord = in_air[:,1];
    ycoord = in_air[:,2];
    ncoord = length(xcoord);
    xcoord_sum = zeros(ncoord);

    xcoord_sum[1] = 0;
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

        function nacaTh(x)
            5*th*(b1*sqrt(x) + b2*x + b3*x^2 + b4*x^3 + b5*x^4)
        end

        for i = 1:ndiv
            thick[i] = nacaTh(x[i])
            thick_slope[i] = ForwardDiff.derivative(nacaTh, x[i])
        end
        rho = 1.1019*th*th*c
        cam[1:ndiv] .= 0.
        cam_slope[1:ndiv] .= 0.

    elseif coord_file[1:7] == "NACAcut"
        m = parse(Int, coord_file[8])/100.
        p = parse(Int, coord_file[9])/10.
        th = parse(Int, coord_file[10:11])/100.
        cut = parse(Int, coord_file[12:13])/100.
        
        b1 = 0.2969; b2 = -0.1260; b3 = -0.3516; b4 = 0.2843; b5 = -0.1015
        for i = 2:ndiv
            x_mod = x[i]*cut
            thick[i] = 5*th*(b1*sqrt(x_mod) + b2*x_mod + b3*x_mod^2 + b4*x_mod^3 + b5*x_mod^4)
            thick_slope[i] = 5*th*(b1/(2*sqrt(x_mod)) + b2 + 2*b3*x_mod + 3*b4*x_mod^2 + 4*b5*x_mod^3)
        end
        x1mod = x[1]*cut
        thick[1] = 5*th*(b1*sqrt(x1mod) + b2*x1mod + b3*x1mod^2 + b4*x1mod^3 + b5*x1mod^4)
        thick_slope[1] = 2*thick_slope[2] - thick_slope[3]
        rho = 1.1019*th*th*c
        cam[1:ndiv] .= 0.
        cam_slope[1:ndiv] .= 0.

    elseif coord_file[1:7] == "NACAwak"
        m = parse(Int, coord_file[8])/100.
        p = parse(Int, coord_file[9])/10.
        th = parse(Int, coord_file[10:11])/100.
        cut = parse(Int, coord_file[12:13])/100.
        
        b1 = 0.2969; b2 = -0.1260; b3 = -0.3516; b4 = 0.2843; b5 = -0.1015
        for i = 2:Int(ndiv/2)
            x_mod = x[i]/cut
            thick[i] = 5*th*(b1*sqrt(x_mod) + b2*x_mod + b3*x_mod^2 + b4*x_mod^3 + b5*x_mod^4)
            thick_slope[i] = 5*th*(b1/(2*sqrt(x_mod)) + b2 + 2*b3*x_mod + 3*b4*x_mod^2 + 4*b5*x_mod^3)
        end
        for i = Int(ndiv/2)+1:ndiv
            thick[i] = 0.
            thick_slope[i] = 0.#(thick[i] - thick[i-1])/(x[i] - x[i-1])
        end
        x1mod = x[1]*cut
        thick[1] = 5*th*(b1*sqrt(x1mod) + b2*x1mod + b3*x1mod^2 + b4*x1mod^3 + b5*x1mod^4)
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

function secant_method(f, x0, x1, ftol, max_iter)

    iter = 2

    f0 = f(x0)
    f1 = f(x1)

    m = length(f0)
    n = length(x0)

    J = zeros(m, n)
    
    while iter <= max_iter

        for i = 1:m
            for j = 1:n
                J[i,j] = (f1[i] - f0[i])/(x1[j] - x0[j])
            end
        end

        #println(det(transpose(J)*J))
        #println(J[1,:])
        x = x1 .- inv(transpose(J)*J)*transpose(J)*f1
        fx = f(x)
        println(iter, "  ", sqrt(sum(fx.^2)))

        if sqrt(sum(fx.^2)) < ftol
            break
        else
            iter += 1
            x0 = x1
            f0 = f1
            x1 = x
            f1 = fx
        end
    end
    return x
end

function getEndCycle(mat, k, str="resultsSummaryEndCycle")

    T = pi/k
    
    end_cycle = mat[end,1]
    
    #Find number of cycles
    ncyc = 0
    for i = 1:1000
        t_l = real(i)*T
        if t_l > end_cycle
            break
        end
        ncyc = ncyc + 1
    end
    ncol = size(mat)[2]
    
    println("Cycles    ", ncyc)

    start_t = real(ncyc-1)*T
    end_t = real(ncyc)*T
    start_ind = argmin(abs.(mat[:,1] .- start_t))
    end_ind = argmin(abs.(mat[:,1] .- end_t))

    nlast = end_ind - start_ind + 1

    newmat = zeros(nlast, ncol)
    newmat[:,1] = (mat[start_ind:end_ind,1] .- start_t)/T
    for i = 2:ncol
        newmat[:,i] = mat[start_ind:end_ind,i]
    end

    f = open(str, "w")
    Serialization.serialize(f, ["#t/T \t", "alpha (deg) \t", "h/c \t", "u/uref \t", "A0 \t", "Cl \t", "Cd \t", "Cm \n"])
    writedlm(f, newmat)
    close(f)

    #Print convergence statistics
    println("Convergence statistics")
    for icyc = 1:ncyc-1
        println("Cycle $(icyc) compared to $ncyc")
        for i = 2:ncol
            x_base = newmat[:,1]
            y_base = newmat[:,i]

            st = real(icyc-1)*T
            et = real(icyc)*T
            si = argmin(abs.(mat[:,1] .- st))
            ei = argmin(abs.(mat[:,1] .- et))

            x_q = (mat[si:ei,1] .- st)/T
            y_q = mat[si:ei,i]
            
            er = errorCalc(x_base, y_base, x_q, y_q)
            er *= 100
            println("Col $i  , $er %")
        end
    end
    
    return newmat
end


function errorCalc(base_x, base_y, quant_x, quant_y; ndiv=1000)
    
    xs = maximum([base_x[1]; quant_x[1]])
    xe = minimum([base_x[end]; quant_x[end]])
    
    x = collect(xs:(xe-xs)/(ndiv-1):xe)
    
    basespl = Spline1D(base_x, base_y)
    quantspl = Spline1D(quant_x, quant_y)
    
    base_o = evaluate(basespl, x)
    quant_o = evaluate(quantspl, x)
    
    rmserror = sqrt(mean(((quant_o .- base_o)./base_o).^2))
end       

function smoothEdges!(q::Array{Float64}, nsm::Int)
    ndiv = length(q)
    
    for i = 1:nsm
        q[nsm-i+1] = 2*q[nsm+2-i] - q[nsm+3-i]
        q[ndiv-nsm+i] = 2*q[ndiv-nsm+i-1] - q[ndiv-nsm+i-2]
    end

    return q
end

function smoothEnd!(q::Array{Float64}, nsm::Int)
    ndiv = length(q)
    
    for i = 1:nsm
        q[ndiv-nsm+i] = 2*q[ndiv-nsm+i-1] - q[ndiv-nsm+i-2]
    end

    return q
end

function smoothStart!(q::Array{Float64}, nsm::Int)
    ndiv = length(q)
    
    for i = 1:nsm
        q[nsm-i+1] = 2*q[nsm+2-i] - q[nsm+3-i]
      end

    return q
end


function smoothScaledEnd!(x::Array{Float64}, q::Array{Float64}, nsm::Int)
    ndiv = length(q)
    
    for i = 1:nsm
        q[ndiv-nsm+i] = q[ndiv-nsm+i-1] + (q[ndiv-nsm+i-1] - q[ndiv-nsm+i-2])*(x[ndiv-nsm+i] - x[ndiv-nsm+i-1])/(x[ndiv-nsm+i-1] - x[ndiv-nsm+i-2])
    end

    return q
end
