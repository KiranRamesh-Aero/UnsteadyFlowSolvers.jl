function calc_y(x::Array{Float64}, eps::Array{Float64}, eta::Array{Float64}, del_eps::Float64, del_eta::Float64, del_t::Float64, m::Int64, n::Int64)
    y = zeros(m+1, n+1)

    for i = 3:m-1
        for j = 3:n-1
            if x[i,j] < 3.1
                epscalc = eps[i]
                etacalc = eta[j]
                y[i,j] = 0
                while true
                    #Find out the location of eps and eta for the interpolation
                    i_l = 1; i_u = m+1; j_l = 1; j_u = n+1
                    for k = 2:m
                        if epscalc >= eps[k] && epscalc < eps[k+1]
                            i_l = k
                            i_u = k+1
                        end
                    end
                    for k = 2:n
                        if etacalc >= eta[k] && etacalc < eta[k+1]
                            j_l = k
                            j_u = k+1
                        end
                    end

                    p_ycal = (epscalc - eps[i_l])/del_eps
                    q_ycal = (etacalc - eta[j_l])/del_eta

                    if j_l == 1
                        #!i_l,j_l
                        a_ycal_p1 = -(x[i_l,j_l+1] - x[i_l,j_l])/del_eta
                        #!i_u,j_l
                        a_ycal_p2 = -(x[i_l+1,j_l+1] - x[i_l+1,j_l])/del_eta
                        #i_u, j_u
                        a_ycal_p3 = -0.5*(x[i_l+1,j_l+2] - x[i_l+1,j_l])/del_eta
                        #i_l, j_u
                        a_ycal_p4 = -0.5*(x[i_l,j_l+2] - x[i_l,j_l])/del_eta
                    elseif j_l == n
                        #i_l, j_l
                        a_ycal_p1 = -0.5*(x[i_l,j_l+1] - x[i_l,j_l-1])/del_eta
                        #i_u, j_l
                        a_ycal_p2 = -0.5*(x[i_l+1,j_l+1] - x[i_l+1,j_l-1])/del_eta
                        #i_u, j_u
                        a_ycal_p3 = -(x[i_l+1,j_l+1] - x[i_l+1,j_l])/del_eta
                        #i_l, j_u
                        a_ycal_p4 = -(x[i_l,j_l+1] - x[i_l,j_l])/del_eta
                    else
                        #i_l, j_l
                        a_ycal_p1 = -0.5*(x[i_l,j_l+1] - x[i_l,j_l-1])/del_eta
                        #i_u, j_l
                        a_ycal_p2 = -0.5*(x[i_l+1,j_l+1] - x[i_l+1,j_l-1])/del_eta
                        #i_u, j_u
                        a_ycal_p3 = -0.5*(x[i_l+1,j_l+2] - x[i_l+1,j_l])/del_eta
                        #i_l, j_u
                        a_ycal_p4 = -0.5*(x[i_l,j_l+2] - x[i_l,j_l])/del_eta
                        end

                        if i_l == 1
                            b_ycal_p1 = (x[i_l+1,j_l] - x[i_l,j_l])/del_eps
                            b_ycal_p2 = 0.5*(x[i_l+2,j_l] - x[i_l,j_l])/del_eps
                            b_ycal_p3 = 0.5*(x[i_l+2,j_l+1] - x[i_l,j_l+1])/del_eps
                            b_ycal_p4 = (x[i_l+1,j_l+1] - x[i_l,j_l+1])/del_eps
                        elseif i_l == m
                            b_ycal_p1 = 0.5*(x[i_l+1,j_l] - x[i_l-1,j_l])/del_eps
                            b_ycal_p2 = (x[i_l+1,j_l] - x[i_l,j_l])/del_eps
                            b_ycal_p3 = (x[i_l+1,j_l+1] - x[i_l,j_l+1])/del_eps
                            b_ycal_p4 = 0.5*(x[i_l+1,j_l+1] - x[i_l-1,j_l+1])/del_eps
                        else
                            b_ycal_p1 = 0.5*(x[i_l+1,j_l] - x[i_l-1,j_l])/del_eps
                            b_ycal_p2 = 0.5*(x[i_l+2,j_l] - x[i_l,j_l])/del_eps
                            b_ycal_p3 = 0.5*(x[i_l+2,j_l+1] - x[i_l,j_l+1])/del_eps
                            b_ycal_p4 = 0.5*(x[i_l+1,j_l+1] - x[i_l-1,j_l+1])/del_eps
                        end

                    ansy2 = (1 - p_ycal)*(1 - q_ycal)*a_ycal_p1 + p_ycal*(1 - q_ycal)*a_ycal_p2 + q_ycal*(1 - p_ycal)*a_ycal_p3 + p_ycal*q_ycal*a_ycal_p4
                    ansy1 = (1-p_ycal)*(1-q_ycal)*b_ycal_p1 + p_ycal*(1 - q_ycal)*b_ycal_p2 + q_ycal*(1-p_ycal)*b_ycal_p3 + p_ycal*q_ycal*b_ycal_p4

                    del_s = del_eps/(sqrt(ansy1^2 + ansy2^2))
                    zet = 2*cos(pi*etacalc/2)*cos(pi*etacalc/2)/pi
                    y[i,j] = y[i,j] + del_s/zet
                    epscalc = epscalc - ansy2*del_s
                    etacalc = etacalc - ansy1*del_s

                    if etacalc <10e-3
                        break
                    end
                end
            end
        end
    end
    return y
end

function calc_shear(u::Array{Float64}, del_eta::Float64, m::Int64)
    tau = zeros(m)

    for i = 2:m
        tau[i] = (2/pi)*((18*u[i,2] - 9*u[i,3] + 2*u[i,4])/(6*del_eta))
    end

    return tau
end

function calc_minnorm(x::Array{Float64}, del_eps::Float64, del_eta::Float64, m::Int64, n::Int64)
    norm = zeros(m, n)

    for i = 2:m
        for j = 2:n
            del_x_eps=0.5*(x[i+1,j] - x[i-1,j])/del_eps
            del_x_eta=0.5*(x[i,j+1] - x[i,j-1])/del_eta
            norm[i,j] = sqrt(del_x_eps^2+del_x_eta^2)
        end
    end
    minnorm = minimum(norm)

    return norm, minnorm
end

function calc_vort(x::Array{Float64}, u::Array{Float64}, del_eps::Float64, del_eta::Float64, m::Int64, n::Int64)
    vort = zeros(m,n)

    for i = 2:m
        for j = 2:n
            del_x_eps = 0.5*(x[i+1,j] - x[i-1,j])/del_eps
            del_x_eta = 0.5*(x[i,j+1] - x[i,j-1])/del_eta
            del_u_eps = 0.5*(u[i+1,j] - u[i-1,j])/del_eps
            del_u_eta = 0.5*(u[i,j+1] - u[i,j-1])/del_eta
            vort[i,j] = -z[j]*(del_x_eps*del_u_eta - del_x_eta*del_u_eps)
        end
    end

    return vort
end

function writeStampBl(dirname::String, x::Array{Float64}, u::Array{Float64}, y::Array{Float64}, tau::Array{Float64}, norm::Array{Float64}, vort::Array{Float64})
    dirvec = readdir()
    if dirname in dirvec
        rm(dirname, recursive=true)
    end
    mkdir(dirname)
    cd(dirname)

    f = open("x", "w")
    writedlm(f, x[2:m, 2:n])
    close(f)

    f = open("u", "w")
    writedlm(f, u[2:m, 2:n])
    close(f)

    f = open("y", "w")
    writedlm(f, y[2:m, 2:n])
    close(f)

    f = open("tau", "w")
    writedlm(f, tau[2:m])
    close(f)

    f = open("vort", "w")
    writedlm(f, vort[2:m, 2:n])
    close(f)

    f = open("norm", "w")
    writedlm(f, norm[2:m, 2:n])
    close(f)
    cd("..")
end


function thomas_alg(a, b, c, d, n)

    a = zeros(n)
    b = zeros(n)
    c = zeros(n)
    d = zeros(n)
    x = zeros(n)
    cc = zeros(n)
    bb = zeros(n)
    dd = zeros(n)

    #Forward sweep
    bb[1] = b[1]
    dd[1] = d[1]
    cc[1:n] = c[1:n]

    for i = 2:n
    bb[i] = b[i] - cc(i-1)*a[i]/bb(i-1)
    end

    for i = 2:n
        dd[i] = d[i] - dd(i-1)*a[i]/bb(i-1)
    end

    x[n] = dd[n]/bb[n]
    for i = n-1:-1:1
    x[i] = (dd[i]-cc[i]*x(i+1))/bb[i]
    end
    return x
end

function tridiag(a,b,c,r,n)
    nmax = 500

    u = zeros(n)

    #Solves for a vector u(1:n) of length n the tridiagonal linear set given by equation (2.4.1) in numerical recipes.
    #a(1:n), b(1:n), c(1:n), and r(1:n) are input vectors and are not modified.

    #Parameter: NMAX is the maximum expected value of n.

    gam = zeros(nmax)

    if b[1] == 0.
        println("rewrite equations")
        #If this happens then you should rewrite your equations as a set of order N, with u2 trivially eliminated.
    end
    bet = b[1]
    u[1] = r[1]/bet

    #Decomposition and forward substitution.
    for j = 2:n
        gam[j] = c[j-1]/bet
        bet = b[j] - a[j]*gam[j]
        if bet == 0.
            println("tridag failed'")
        end
        u[j] = (r[j] - a[j]*u[j-1])/bet
    end

    #Backsubstitution
    for j = n-1:-1:1
        u[j] = u[j] - gam[j+1]*u[j+1]
    end
    return u
end

function calc_eqncoeff(x::Array{Float64}, x_prev::Array{Float64}, z::Array{Float64}, del_eps::Float64, del_eta::Float64, m::Int64, n::Int64)
    p = zeros(m+1, n+1)
    q = zeros(m+1, n+1)
    r = zeros(m+1, n+1)
    s = zeros(m+1, n+1)
    t = zeros(m+1, n+1)

    for i = 2:m
        for j = 2:n
            #Calculate all required partial derivatives of xbar
            d_xbar_eps = 0.25*(x[i+1,j] - x[i-1,j] + x_prev[i+1,j] - x_prev[i-1,j])/del_eps
            d_xbar_eta = 0.25*(x[i,j+1] - x[i,j-1] + x_prev[i,j+1] - x_prev[i,j-1])/del_eta
            d2_xbar_eps = 0.5*(x[i+1,j] - 2*x[i,j] + x[i-1,j] + x_prev[i+1,j] - 2*x_prev[i,j] + x_prev[i-1,j])/del_eps^2
            d2_xbar_eta = 0.5*(x[i,j+1] - 2*x[i,j] + x[i,j-1] + x_prev[i,j+1] - 2*x_prev[i,j] + x_prev[i,j-1])/del_eta^2
            d2_xbar_epseta = 0.125*(x[i+1,j+1]  -x[i+1,j-1] - x[i-1,j+1] + x[i-1,j-1] + x_prev[i+1,j+1] - x_prev[i+1,j-1]
                                    - x_prev[i-1,j+1] + x_prev[i-1,j-1])/(del_eta*del_eps)

            #Calculate the equation coefficients (P, Q ....)
            p[i,j] = z[j]*z[j]*(d_xbar_eps*(gam[j]*d_xbar_eps +d2_xbar_epseta) - d_xbar_eta*d2_xbar_eps)
            t[i,j] = z[j]*z[j]*d_xbar_eps^2
	    s[i,j] = z[j]*z[j]*(-2*d_xbar_eps*d_xbar_eta)
            r[i,j] = z[j]*z[j]*d_xbar_eta^2
            q[i,j] = z[j]*z[j]*(d_xbar_eta*d2_xbar_epseta - d_xbar_eps*(gam[j]*d_xbar_eta + d2_xbar_eta))
        end
    end

    i = 1
    for j = 2:n
        d_xbar_eps = 0.5*(x[i+1,j] - x[i,j] + x_prev[i+1,j] - x_prev[i,j])/del_eps
        d_xbar_eta = 0.25*(x[i,j+1] - x[i,j-1] + x_prev[i,j+1] - x_prev[i,j-1])/del_eta
        d2_xbar_eps = 0.5*(x[i+2,j] - 2*x[i+1,j] + x[i,j] + x_prev[i+2,j] -2 *x_prev[i+1,j] + x_prev[i,j])/del_eps^2
        d2_xbar_eta = 0.5*(x[i,j+1] - 2*x[i,j] + x[i,j-1] + x_prev[i,j+1] - 2*x_prev[i,j] + x_prev[i,j-1])/del_eta^2
        d2_xbar_epseta = 0.25*(x[i+1,j+1] - x[i+1,j-1] - x[i,j+1]+x[i,j-1] + x_prev[i+1,j+1] - x_prev[i+1,j-1] - x_prev[i,j+1]+x_prev[i,j-1])/(del_eta*del_eps)

        #Calculate the equation coefficients (P, Q ....)
        p[i,j] = z[j]*z[j]*(d_xbar_eps*(gam[j]*d_xbar_eps + d2_xbar_epseta) - d_xbar_eta*d2_xbar_eps)
        t[i,j] = z[j]*z[j]*d_xbar_eps^2
    end

    i = m+1
    for j = 2:n
        d_xbar_eps = 0.5*(x[i,j] - x[i-1,j] + x_prev[i,j] - x_prev[i-1,j])/del_eps
        d_xbar_eta = 0.25*(x[i,j+1] - x[i,j-1] + x_prev[i,j+1] - x_prev[i,j-1])/del_eta
        d2_xbar_eps = 0.5*(x[i,j] - 2*x[i-1,j] + x[i-1,j] + x_prev[i,j] - 2*x_prev[i-1,j] + x_prev[i-2,j])/del_eps^2
        d2_xbar_eta = 0.5*(x[i,j+1] - 2*x[i,j] + x[i,j-1] + x_prev[i,j+1] - 2*x_prev[i,j] + x_prev[i,j-1])/del_eta^2
        d2_xbar_epseta = 0.25*(x[i,j+1] - x[i,j-1] - x[i-1,j+1] + x[i-1,j-1] + x_prev[i,j+1] - x_prev[i,j-1] - x_prev[i-1,j+1] + x_prev[i-1,j-1])/(del_eta*del_eps)

        #Calculate the equation coefficients (P, Q ....)
        p[i,j] = z[j]*z[j]*(d_xbar_eps*(gam[j]*d_xbar_eps + d2_xbar_epseta) - d_xbar_eta*d2_xbar_eps)
        t[i,j] = z[j]*z[j]*d_xbar_eps^2
    end

    return p, q, r, s, t
end

function calc_rhs(u::Array{Float64}, u_prev::Array{Float64}, p::Array{Float64}, q::Array{Float64}, r::Array{Float64}, s::Array{Float64}, t::Array{Float64}, del_eps::Float64, del_eta::Float64, m::Int64, n::Int64)
    d = zeros(m+1, n+1)

    for i = 2:m
        for j = 2:n
            del2_uprev_eps = u_prev[i+1,j] -2 *u_prev[i,j] + u_prev[i-1,j]
            del2_uprev_eta = u_prev[i,j+1] - 2*u_prev[i,j] + u_prev[i,j-1]
            del2_ubar_epseta = 0.125*(u[i+1,j+1] - u[i+1,j-1] - u[i-1,j+1] + u[i-1,j-1] + u_prev[i+1,j+1] - u_prev[i-1,j+1] -u_prev[i+1,j-1] + u_prev[i-1,j-1])/(del_eps*del_eta)

            if q[i,j] >= 0
                ximns_uprev_eps = u_prev[i,j] - u_prev[i-1,j]
            else
                ximns_uprev_eps = u_prev[i+1,j] - u_prev[i,j]
            end

            if p[i,j] >= 0
                ximns_uprev_eta = u_prev[i,j] - u_prev[i,j-1]
            else
                ximns_uprev_eta = u_prev[i,j+1] - u_prev[i,j]
            end

            #RHS specific for cylinder (includes terms from Uinf)
            d[i,j] = u_prev[i,j] + del_t*(0.5*((r[i,j]*del2_uprev_eps/del_eps^2)
             + (q[i,j]*ximns_uprev_eps/del_eps) + (t[i,j]*del2_uprev_eta/del_eta^2)
             + (p[i,j]*ximns_uprev_eta/del_eta)) + sin(0.5*(x[i,j] + x_prev[i,j]))*cos(0.5*(x[i,j] + x_prev[i,j]))
             + s[i,j]*del2_ubar_epseta)
        end
    end

    i = 1
    for j = 2:n
        del2_uprev_eta = u_prev[i,j+1] - 2*u_prev[i,j] + u_prev[i,j-1]

        if p[i,j] >= 0
            ximns_uprev_eta = u_prev[i,j] - u_prev[i,j-1]
        else
            ximns_uprev_eta = u_prev[i,j+1] - u_prev[i,j]
        end

        d[i,j] = u_prev[i,j] + del_t*((p[i,j]*ximns_uprev_eta/del_eta) + (t[i,j]*del2_uprev_eta/del_eta^2))/2
    end

    i = m+1
    for j = 2:n
        del2_uprev_eta = u_prev[i,j+1] - 2*u_prev[i,j] + u_prev[i,j-1]

        if p[i,j] >= 0
            ximns_uprev_eta = u_prev[i,j]-u_prev[i,j-1]
        else
            ximns_uprev_eta = u_prev[i,j+1]-u_prev[i,j]
        end

        d[i,j] = u_prev[i,j] + del_t*((p[i,j]*ximns_uprev_eta/del_eta) + (t[i,j]*del2_uprev_eta/del_eta^2))/2
    end

    return d
end

function calc_triag_util(r::Array{Float64}, q::Array{Float64}, d::Array{Float64}, util::Array{Float64}, del_eps::Float64, del_eta::Float64, del_t::Float64, m::Int64, n::Int64)
    inf_util_a = zeros(m+1)
    inf_util_b = zeros(m+1)
    inf_util_c = zeros(m+1)

    for j = 2:n
        mul_util = del_t/(2*del_eps*del_eps)
        for i = 2:m
            if q[i,j] >= 0
                inf_util_b[i] = 1 - mul_util*(-2*r[i,j] - del_eps*q[i,j])
                inf_util_c[i] = -mul_util*(r[i,j] + del_eps*q[i,j])
                inf_util_a[i] = -mul_util*r[i,j]
            else
                inf_util_b[i] = 1 - mul_util*(-2*r[i,j] + del_eps*q[i,j])
                inf_util_c[i] = -mul_util*r[i,j]
                inf_util_a[i] = -mul_util*(r[i,j] - del_eps*q[i,j])
            end
        end

        #Use boundary conditions
        d[2,j] = d[2,j] - util[1,j]*inf_util_a[2]
        d[m,j] = d[m,j] - inf_util_c[m]*util[m+1,j]

        util[2:m,j] = tridiag(inf_util_a[2:m], inf_util_b[2:m], inf_util_c[2:m], d[2:m,j],  m-1)
    end

    return util
end

function calc_triag_u(p::Array{Float64}, t::Array{Float64}, d::Array{Float64}, util::Array{Float64}, u::Array{Float64}, del_eps::Float64, del_eta::Float64, del_t::Float64, m::Int64, n::Int64)
    inf_u_a = zeros(n+1)
    inf_u_b = zeros(n+1)
    inf_u_c = zeros(n+1)

    for i = 1:m+1
        mul_u = del_t/(2*del_eta*del_eta)
        for j = 2:n
            if p[i,j] >= 0
                inf_u_b[j] = 1 - mul_u*(-2*t[i,j] - del_eta*p[i,j])
                inf_u_c[j] = -mul_u*(t[i,j] + del_eta*p[i,j])
                inf_u_a[j] = -mul_u*t[i,j]
            else
                inf_u_b[j] = 1 - mul_u*(-2*t[i,j] + del_eta*p[i,j])
                inf_u_c[j] = -mul_u*t[i,j]
                inf_u_a[j] = -mul_u*(t[i,j] - del_eta*p[i,j])
            end
        end

        #Use BCs
        util[i,2] = util[i,2] - inf_u_a[2]*u[i,1]
        util[i,n] = util[i,n] - inf_u_c[n]*u[i,n+1]

        u[i,2:n] = tridiag(inf_u_a[2:n], inf_u_b[2:n], inf_u_c[2:n], util[i,2:n], n-1)
    end

    return u
end


function blLagrangian(dtstar = 0.001, tTot = 1.0, m=200, n=60, tol=1e-6,
     dels=1e-3, writeflag = 0, writeInterval = 0.5; maxwrite=100, nround=6)

## Needs to be de

    # if writeflag is on, determine the timesteps to write at
    if writeflag == 1
        writeArray = Int64[]

        for i = 1:maxwrite
            tcur = writeInterval*real(i)
            if tcur > tTot
                break
            else
                push!(writeArray, Int(round(tcur/dtstar)))
            end
        end
    end

    #Define domain in eta and eps
    del_eps = pi/m
    del_eta = 1./n
    eps = [linspace(0., pi, m+1);]
    eta = [linspace(0, 1., n+1);]

    x = zeros(m+1, n+1)
    u = zeros(m+1, n+1)
    util = zeros(m+1, n+1)
    x_prev = zeros(m+1, n+1)
    u_prev = zeros(m+1, n+1)
    u_previter = zeros(m+1, n+1)
    newt = zeros(m+1, n+1)
    z = zeros(n+1)
    gam = zeros(n+1)

    #Assign initial conditions
    del_t = dtstar

    for i = 1:m+1
        x[i,1:n+1] = eps[i]
    end

    for i = 1:m+1
        u[i,1:n+1] = sin(eps[i])
    end

    #Precalculations
    for j = 1:n+1
        z[j] = 2*cos(pi*eta[j]/2)*cos(pi*eta[j]/2)/pi
        gam[j] = -pi*tan(pi*eta[j]/2)
    end

    tim = 0.

    p = q = r = s = t = d = zeros(m+1, n+1)

    #Start time loop
    for istep = 1:5#3200
        tim = tim + del_t

        #Store values of u and x from previous time step
        u_prev = u
        x_prev = x

        #Start iteration loop
        iter=0
        iter_res = 1.
        while iter_res > tol
            iter +=  1
            println(iter, iter_res)

            p, q, r, s, t = calc_eqncoeff(x, x_prev, z, del_eps, del_eta, m, n)

            d = calc_rhs(u, u_prev, p, q, r, s, t, del_eps, del_eta, m, n)

            #Define BCs for util and u
            #Edge velocity is specific for cylinder
            util[1,2:n] = d[1,2:n]
            util[m+1,2:n] = d[m+1,2:n]
            u[1:m+1,1] = 0.
            u[1:m+1,n+1] = sin.(x[1:m+1,n+1])

            util = calc_triag_util(r, q, d, util, del_eps, del_eta, del_t, m, n)

            #Store values of u from previous iteration
            u_previter=u

            u = calc_triag_u(p, t, d, util, u, del_eps, del_eta, del_t, m, n)

            #Update x
            for i = 1:m+1
                for j = 1:n+1
                    x[i,j] = x_prev[i,j] + del_t*(u[i,j] + u_prev[i,j])/2
                end
            end

            #Check for convergence
            newt[:,:] = 0
            for i = 2:m
                for j = 2:n
                    if u[i,j] != 0
                        newt[i,j] = abs((u[i,j] - u_previter[i,j])/u[i,j])
                    end
                end
            end
            iter_res =  maximum(newt)
        end

        #Write flow details if required
        if writeflag == 1
            if istep in writeArray
                dirname = "$(round(tim,nround))"

                norm = vort = y = zeros(m,n)
                tau = zeros(m)

                y = calc_y(x, eps, eta, del_eps, del_eta, del_t, m, n)
                tau = calc_shear(u, del_eta, m)
                norm, minnorm = calc_minnorm(c, del_eps, del_eta, m, n)
                println(tim, minnorm)
                vort = calc_vort(x, u, del_eps, del_eta, m, n)

                writeStampBl(dirname, x, u, y, tau, norm, vort)
            else
                _,minnorm = calc_minnorm(c, del_eps, del_eta, m, n)
                println(tim, "    ", minnorm)
            end
        end

        #Check for separation (both gradients equal zero)
        for i = 2:m
            for j = 2:n
                del_x_eps=0.5*(x[i+1,j] - x[i-1,j])/del_eps
                del_x_eta=0.5*(x[i,j+1] - x[i,j-1])/del_eta
                if del_x_eps<tol && del_x_eta<tol
                    println("singularity (separation) detected. time = $tim,
                        i=$i, j=$j, eps=$(eps[i]), eta=$(eta[j]), x=$(x[i,j])")
                    break
                end
            end
        end
    end
end
