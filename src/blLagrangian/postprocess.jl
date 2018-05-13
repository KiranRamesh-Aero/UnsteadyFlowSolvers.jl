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
            norm[i,j] = sqrt(del_x_eps**2+del_x_eta**2)
        end
    end
    minnorm = mimimum(norm)

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
            vort[i,j] = -z(j)*(del_x_eps*del_u_eta-del_x_eta*del_u_eps)
        end
    end

    return vort
end

function writeStamp(dirname::String, x::Array{Float64}, u::Array{Float64}, y::Array{Float64}, tau::Array{Float64}, norm::Array{Float64}, vort::Array{Float64})
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
