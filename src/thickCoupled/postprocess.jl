function writeStamp(dirname::String, t::Float64, surf::TwoDSurfThick, curfield::TwoDFlowField, qu::Vector{Float64}, ql::Vector{Float64}, cpu::Vector{Float64}, cpl::Vector{Float64}, suc::Vector{Float64}, del::Vector{Float64}, E::Vector{Float64}, thick_orig::Vector{Float64}, quc::Vector{Float64}, qucx::Vector{Float64}, quct::Vector{Float64})
    dirvec = readdir()
    if dirname in dirvec
        rm(dirname, recursive=true)
    end
    mkdir(dirname)
    cd(dirname)

    f = open("timeKinem", "w")
    Serialization.serialize(f, ["#time \t", "alpha (deg) \t", "h/c \t", "u/uref \t"])
    writedlm(f, [t, surf.kinem.alpha, surf.kinem.h, surf.kinem.u]')
    close(f)

    f = open("FourierCoeffsA", "w")
    Serialization.serialize(f, ["#Fourier coeffs (0-n) \t", "d/dt (Fourier coeffs) \n"])
    mat = zeros(surf.naterm+1, 2)
    mat[:,1] = [surf.a0[1];surf.aterm[:]]
    mat[:,2] = [surf.a0dot[1];surf.adot[:]]
    writedlm(f, mat)
    close(f)

    f = open("FourierCoeffsB", "w")
    Serialization.serialize(f, ["#Fourier coeffs (0-n) \t", "d/dt (Fourier coeffs) \n"])
    mat = zeros(surf.naterm, 1)
    mat[:,1] = surf.bterm[:]
    writedlm(f, mat)
    close(f)

    # f = open("forces", "w")
    # write(f, ["#cl \t", "cd \t", "cm \t", "Gamma \t", "cn \t", "cs \t", "cnc \t",
    #           "cnnc \t", "nonl \t", "cm_n \t", "cm_pvt \t", "nonl_m \n"])
    # cl, cd, cm, gamma, cn, cs, cnc, cncc, nonl, cm_n, cm_pvt, nonl_m = calc_forces_more(surf)
    # writedlm(f, [cl, cd, cm, gamma, cn, cs, cnc, cnnc, nonl, cm_n, cm_pvt, nonl_m])
    # close(f)

    f = open("tev", "w")
    Serialization.serialize(f, ["#strength \t", "x-position \t", "z-position \n"])
    tevmat = zeros(length(curfield.tev), 3)
    for i = 1:length(curfield.tev)
        tevmat[i,:] = [curfield.tev[i].s curfield.tev[i].x curfield.tev[i].z]
    end
    writedlm(f, tevmat)
    close(f)

    f = open("lev", "w")
    Serialization.serialize(f, ["#strength \t", "x-position \t", "z-position \n"])
    levmat = zeros(length(curfield.lev), 3)
    for i = 1:length(curfield.lev)
        levmat[i,:] = [curfield.lev[i].s curfield.lev[i].x curfield.lev[i].z]
    end
    writedlm(f, levmat)
    close(f)

    f = open("boundv", "w")
    Serialization.serialize(f, ["#strength \t", "x-position \t", "z-position \n"])
    bvmat = zeros(length(surf.bv), 3)
    for i = 1:length(surf.bv)
        bvmat[i,:] = [surf.bv[i].s surf.bv[i].x surf.bv[i].z]
    end
    writedlm(f, bvmat)
    close(f)

    f = open("src", "w")
    Serialization.serialize(f, ["#strength \t", "x-position \t", "z-position \n"])
    bvmat = zeros(length(surf.src), 3)
    for i = 1:length(surf.src)
        bvmat[i,:] = [surf.src[i].s surf.src[i].x surf.src[i].z]
    end
    writedlm(f, bvmat)
    close(f)

    f = open("cp_edgevel", "w")
    Serialization.serialize(f, ["#x \t", "cp_us \t", "cp_ls \t", "qu \t", "ql \n"])
    mat = zeros(length(surf.x), 5)
    mat[:,1] = surf.x[:]
    mat[:,2] = cpu[:]
    mat[:,3] = cpl[:]
    mat[:,4] = qu[:]
    mat[:,5] = ql[:]
    writedlm(f, mat)
    close(f)

    f = open("bl_properties", "w")
    Serialization.serialize(f, ["#x \t", "sc \t", "del \t", "E \t", "thick_orig \t", "thick \t", "q_c \t", "q_c_x \t", "q_c_t \n"])
    mat = zeros(length(surf.x)-1, 9)
    mat[:,1] = surf.x[2:end]
    mat[:,2] = suc[:]
    mat[:,3] = del[:]
    mat[:,4] = E[:]
    mat[:,5] = thick_orig[2:end]
    mat[:,6] = surf.thick[2:end]
    mat[:,7] = quc[:]
    mat[:,8] = qucx[:]
    mat[:,9] = quct[:]
    writedlm(f, mat)
    close(f)

    cd("..")
end


