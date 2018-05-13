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

function calc_eqncoeff(x::Array{Float64}, x_prev::Array{Float64}, z::Array{Float64}, del_eps::Float64, del_eta::Float64, m::Inte64, n::Int64)
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

function calc_rhs(u::Arrat{Float64}, u_prev::Array{Float64}, p::Array{Float64}, q::Array{Float64}, r::Array{Float64}, s::Array{Float64}, t::Array{Float64}, del_eps::Float64, del_eta::Float64, m::Int64, n::Int64)
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
            d[i,j] = u_prev[i,j] + del_t*(0.5*((r[i,j]*del2_uprev_eps/del_eps^2) + (q[i,j]*ximns_uprev_eps/del_eps) + (t[i,j]*del2_uprev_eta/del_eta^2) + (p[i,j]*ximns_uprev_eta/del_eta)) + 0.5*(sin(x[i,j])*cos(x[i,j]) + sin(x_prev[i,j])*cos(x_prev[i,j])) + s[i,j]*del2_ubar_epseta)
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

function calc_triag_util(p::Array{Float64}, t::Array{Float64}, d::Array{Float64}, util::Array{Float64}, u::Array{Float64}, del_eps::Float64, del_eta::Float64, del_t::Float64, m::Int64, n::Int64)
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
