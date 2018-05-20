function writeStamp(dirname::String, t::Float64, surf::ThreeDSurfSimple, curfield::ThreeDFieldSimple)
    dirvec = readdir()
    if dirname in dirvec
        rm(dirname, recursive=true)
    end
    mkdir(dirname)
    cd(dirname)

    f = open("time", "w")
    write(f, ["#time \n"])
    writedlm(f, [t]')
    close(f)

    f = open("spanVariations", "w")
    write(f, ["#y location \t", "alpha (def) \t", "h \t", "u \t", "scaled bound circ \t",
    "a0 \t", "a03d \t", "tevstr \t", "cl \t", "cd \t", "cm \n"])

    mat = zeros(surf.nspan,11)
    for i = 1:surf.nspan
        mat[i,1] = surf.yle[i]
        mat[i,2] = surf.s2d[i].kinem.alpha*180/pi
        mat[i,3] = surf.s2d[i].kinem.h
        mat[i,4] = surf.s2d[i].kinem.u
        mat[i,5] = surf.bc[i]
        mat[i,6] = surf.s2d[i].a0[1]
        mat[i,7] = surf.a03d[i]
        mat[i,8] = curfield.f2d[i].tev[end].s
        mat[i,9] = surf.fc[i,1]
        mat[i,10] = surf.fc[i,2]
        mat[i,11] = surf.fc[i,3]
    end
    writedlm(f, mat)
    close(f)

    f = open("FourierCoeffsRoot", "w")
    write(f, ["#y=$(surf.yle[surf.nspan]) \n"])
    write(f, ["#Fourier coeffs at root (0-n) \t", "d/dt (Fourier coeffs at root) \n"])
    mat = zeros(surf.naterm+1, 2)
    mat[:,1] = [surf.s2d[surf.nspan].a0[1]; surf.s2d[surf.nspan].aterm[:]]
    mat[1:4,2] = [surf.s2d[surf.nspan].a0dot[1]; surf.s2d[surf.nspan].adot[:]]
    writedlm(f, mat)
    close(f)

    f = open("FourierCoeffsMid", "w")
    nmid = Int(round(surf.nspan/2))
    write(f, ["#y=$(surf.yle[nmid]) \n"])
    write(f, ["#Fourier coeffs at root (0-n) \t", "d/dt (Fourier coeffs at root) \n"])
    mat = zeros(surf.naterm+1, 2)
    mat[:,1] = [surf.s2d[nmid].a0[1]; surf.s2d[nmid].aterm[:]]
    mat[1:4,2] = [surf.s2d[nmid].a0dot[1];surf.s2d[nmid].adot[:]]
    writedlm(f, mat)
    close(f)

    f = open("FourierCoeffsTip", "w")
    write(f, ["#y=$(surf.yle[1]) \n"])
    write(f, ["#Fourier coeffs at root (0-n) \t", "d/dt (Fourier coeffs at root) \n"])
    mat = zeros(surf.naterm+1, 2)
    mat[:,1] = [surf.s2d[1].a0[1]; surf.s2d[1].aterm[:]]
    mat[1:4,2] = [surf.s2d[1].a0dot[1];surf.s2d[1].adot[:]]
    writedlm(f, mat)
    close(f)

    # f = open("forces", "w")
    # write(f, ["#cl \t", "cd \t", "cm \t", "Gamma \t", "cn \t", "cs \t", "cnc \t",
    #           "cnnc \t", "nonl \t", "cm_n \t", "cm_pvt \t", "nonl_m \n"])
    # cl, cd, cm, gamma, cn, cs, cnc, cncc, nonl, cm_n, cm_pvt, nonl_m = calc_forces_more(surf)
    # writedlm(f, [cl, cd, cm, gamma, cn, cs, cnc, cnnc, nonl, cm_n, cm_pvt, nonl_m])
    # close(f)

    f = open("tevRoot", "w")
    write(f, ["#y=$(surf.yle[surf.nspan]) \n"])
    write(f, ["#strength \t", "x-position \t", "z-position \n"])
    ntev = length(curfield.f2d[surf.nspan].tev)
    mat = zeros(ntev, 3)
    for i = 1:ntev
        mat[i,:] = [curfield.f2d[surf.nspan].tev[i].s curfield.f2d[surf.nspan].tev[i].x curfield.f2d[surf.nspan].tev[i].z]
    end
    writedlm(f, mat)
    close(f)

    f = open("tevMid", "w")
    write(f, ["#y=$(surf.yle[nmid]) \n"])
    write(f, ["#strength \t", "x-position \t", "z-position \n"])
    ntev = length(curfield.f2d[nmid].tev)
    mat = zeros(ntev, 3)
    for i = 1:ntev
        mat[i,:] = [curfield.f2d[nmid].tev[i].s curfield.f2d[nmid].tev[i].x curfield.f2d[nmid].tev[i].z]
    end
    writedlm(f, mat)
    close(f)

    f = open("tevTip", "w")
    write(f, ["#y=$(surf.yle[1]) \n"])
    write(f, ["#strength \t", "x-position \t", "z-position \n"])
    ntev = length(curfield.f2d[1].tev)
    mat = zeros(ntev, 3)
    for i = 1:ntev
        mat[i,:] = [curfield.f2d[1].tev[i].s curfield.f2d[1].tev[i].x curfield.f2d[1].tev[i].z]
    end
    writedlm(f, mat)
    close(f)

    f = open("levRoot", "w")
    write(f, ["#y=$(surf.yle[surf.nspan]) \n"])
    write(f, ["#strength \t", "x-position \t", "z-position \n"])
    nlev = length(curfield.f2d[surf.nspan].lev)
    mat = zeros(nlev, 3)
    for i = 1:nlev
        mat[i,:] = [curfield.f2d[surf.nspan].lev[i].s curfield.f2d[surf.nspan].lev[i].x curfield.f2d[surf.nspan].lev[i].z]
    end
    writedlm(f, mat)
    close(f)

    f = open("levMid", "w")
    write(f, ["#y=$(surf.yle[nmid]) \n"])
    write(f, ["#strength \t", "x-position \t", "z-position \n"])
    nlev = length(curfield.f2d[nmid].lev)
    mat = zeros(nlev, 3)
    for i = 1:nlev
        mat[i,:] = [curfield.f2d[nmid].lev[i].s curfield.f2d[nmid].lev[i].x curfield.f2d[nmid].lev[i].z]
    end
    writedlm(f, mat)
    close(f)

    f = open("levTip", "w")
    write(f, ["#y=$(surf.yle[1]) \n"])
    write(f, ["#strength \t", "x-position \t", "z-position \n"])
    nlev = length(curfield.f2d[1].lev)
    mat = zeros(nlev, 3)
    for i = 1:nlev
        mat[i,:] = [curfield.f2d[1].lev[i].s curfield.f2d[1].lev[i].x curfield.f2d[1].lev[i].z]
    end
    writedlm(f, mat)
    close(f)

    f = open("boundvRoot", "w")
    write(f, ["#y=$(surf.yle[surf.nspan]) \n"])
    write(f, ["#strength \t", "x-position \t", "z-position \n"])
    nbv = length(surf.s2d[surf.nspan].bv)
    mat = zeros(nbv, 3)
    for i = 1:nbv
        mat[i,:] = [surf.s2d[surf.nspan].bv[i].s surf.s2d[surf.nspan].bv[i].x surf.s2d[surf.nspan].bv[i].z]
    end
    writedlm(f, mat)
    close(f)

    f = open("boundvMid", "w")
    write(f, ["#y=$(surf.yle[nmid]) \n"])
    write(f, ["#strength \t", "x-position \t", "z-position \n"])
    nbv = length(surf.s2d[nmid].bv)
    mat = zeros(nbv, 3)
    for i = 1:nbv
        mat[i,:] = [surf.s2d[nmid].bv[i].s surf.s2d[nmid].bv[i].x surf.s2d[nmid].bv[i].z]
    end
    writedlm(f, mat)
    close(f)

    f = open("boundvTip", "w")
    write(f, ["#y=$(surf.yle[1]) \n"])
    write(f, ["#strength \t", "x-position \t", "z-position \n"])
    nbv = length(surf.s2d[1].bv)
    mat = zeros(nbv, 3)
    for i = 1:nbv
        mat[i,:] = [surf.s2d[1].bv[i].s surf.s2d[1].bv[i].x surf.s2d[1].bv[i].z]
    end
    writedlm(f, mat)
    close(f)

    cd("..")
end

function viewVort3Dstrip(ztot::Float64, tevR::Array{Float64}, levR::Array{Float64}, bvR::Array{Float64},
    tevM::Array{Float64}, levM::Array{Float64}, bvM::Array{Float64},
    tevT::Array{Float64}, levT::Array{Float64}, bvT::Array{Float64})

    scatter(tevT[:,2], tevT[:,3] + ztot/9., s=5, c=tevT[:,1], edgecolors="none")
    sc = scatter(levT[:,2], levT[:,3] + ztot/9., s=5, c=levT[:,1], edgecolors="none")
    sc2 = scatter(bvT[:,2], bvT[:,3] + ztot/9., s=5, c=bvT[:,1], edgecolors="none")
    plot(bvT[:,2], bvT[:,3] + ztot/9., color = "black", linewidth=1.0)

    scatter(tevM[:,2], tevM[:,3] + 4.5*ztot/9., s=5, c=tevM[:,1], edgecolors="none")
    sc = scatter(levM[:,2], levM[:,3] + 4.5*ztot/9., s=5, c=levM[:,1], edgecolors="none")
    sc2 = scatter(bvM[:,2], bvM[:,3] + 4.5*ztot/9., s=5, c=bvM[:,1], edgecolors="none")
    plot(bvT[:,2], bvM[:,3] + ztot/9., color = "black", linewidth=1.0)

    scatter(tevR[:,2], tevR[:,3] + 8.*ztot/9., s=5, c=tevR[:,1], edgecolors="none")
    sc = scatter(levR[:,2], levR[:,3] + 8.*ztot/9., s=5, c=levR[:,1], edgecolors="none")
    sc2 = scatter(bvR[:,2], bvR[:,3] + 8.*ztot/9., s=5, c=bvR[:,1], edgecolors="none")
    plot(bvT[:,2], bvR[:,3] + ztot/9., color = "black", linewidth=1.0)
end

function makeVortPlots3Dstrip()
    dirvec = readdir()
    dirresults = map(x->(v = tryparse(Float64,x); isnull(v) ? 0.0 : get(v)),dirvec)
    #Determine axis limits
    dirmax = maximum(dirresults)
    cd("$(dirmax)")
    tev = readdlm("tevRoot")
    lev = try
        readdlm("levRoot")
    catch
        Array{Float64}(0,3)
    end
    bv = readdlm("boundvRoot")

    xmin = minimum([tev[:,2];lev[:,2];bv[:,2];])
    zmin = minimum([tev[:,3];lev[:,3];bv[:,3];])
    xmax = maximum([tev[:,2];lev[:,2];])
    zmax = maximum([bv[:,2];tev[:,3];lev[:,3];bv[:,3];])

    zdist = zmax - zmin
    ztot = 9*zdist

    cd("..")

    if "vortPlots" in dirvec
        rm("vortPlots", recursive=true)
    end
    mkdir("vortPlots")
    for i=1:length(dirresults)
        if dirresults[i] != 0
            dirstr="$(dirresults[i])"
            cd(dirstr)
            tevR = readdlm("tevRoot")
            levR = try
                readdlm("levRoot")
            catch
                Array{Float64}(0,3)
            end
            bvR = readdlm("boundvRoot")

            tevM = readdlm("tevMid")
            levM = try
                readdlm("levMid")
            catch
                Array{Float64}(0,3)
            end
            bvM = readdlm("boundvMid")

            tevT = readdlm("tevTip")
            levT = try
                readdlm("levTip")
            catch
                Array{Float64}(0,3)
            end
            bvT = readdlm("boundvTip")

            viewVort3Dstrip(ztot, tevR, levR, bvR, tevM, levM, bvM, tevT, levT, bvT)
            axis([xmin-1, xmax+1, -1, ztot+1])
            savefig("../vortPlots/$(dirresults[i]).png")
            close()
            cd("..")
        end
    end
end

function makeForcePlots3Dstrip()
    dirvec = readdir()
    if "forcePlots" in dirvec
        rm("forcePlots", recursive=true)
    end
    mkdir("forcePlots")

    mat = readdlm("resultsSummary")

    t = mat[:,1]
    cl = mat[:,3]
    cd = mat[:,4]
    cm = mat[:,5]

    len = length(t)

    plot(t, cl)
    range = Int(round(0.05*len)):Int(round(0.95*len))
    xmin = t[1]
    xmax = t[end]
    zmin = minimum(cl[range]) - 0.1*abs(minimum(cl[range]))
    zmax = maximum(cl[range]) + 0.1*abs(maximum(cl[range]))
    axis([xmin, xmax, zmin, zmax])
    xlabel(L"$t^*$")
    ylabel(L"$C_L$")
    savefig("forcePlots/cl.png")
    close()

    plot(t, cd)
    range = Int(round(0.05*len)):Int(round(0.95*len))
    xmin = t[1]
    xmax = t[end]
    zmin = minimum(cd[range]) - 0.1*abs(minimum(cd[range]))
    zmax = maximum(cd[range]) + 0.1*abs(maximum(cd[range]))
    axis([xmin, xmax, zmin, zmax])
    xlabel(L"$t^*$")
    ylabel(L"$C_D$")
    savefig("forcePlots/cd.png")
    close()

    plot(t, cm)
    range = Int(round(0.05*len)):Int(round(0.95*len))
    xmin = t[1]
    xmax = t[end]
    zmin = minimum(cm[range]) - 0.1*abs(minimum(cm[range]))
    zmax = maximum(cm[range]) + 0.1*abs(maximum(cm[range]))
    axis([xmin, xmax, zmin, zmax])
    xlabel(L"$t^*$")
    ylabel(L"$C_M$")
    savefig("forcePlots/cm.png")
    close()
end

function makeTevstrPlots3Dstrip()
    dirvec = readdir()
    dirresults = map(x->(v = tryparse(Float64,x); isnull(v) ? 0.0 : get(v)),dirvec)

    if "TevstrPlots" in dirvec
        rm("TevstrPlots", recursive=true)
    end
    mkdir("TevstrPlots")

    #Determine axis range over all the time steps
    strmax = 0; strmin = 0
    for i=1:length(dirresults)
        if dirresults[i] != 0
            dirstr="$(dirresults[i])"
            cd(dirstr)
            mat = readdlm("spanVariations")
            y = mat[:,1]
            tevstr = mat[:,8]
            strmax = maximum([strmax; maximum(tevstr)])
            strmin = minimum([strmin; minimum(tevstr)])
            cd("..")
        end
    end

    for i=1:length(dirresults)
        if dirresults[i] != 0
            dirstr="$(dirresults[i])"
            cd(dirstr)
            mat = readdlm("spanVariations")
            y = mat[:,1]
            tevstr = mat[:,8]
            plot(y,tevstr,"k-")
            axis([y[1],y[end],strmin, strmax])
            xlabel(L"$y_{LE}$")
            ylabel(L"$\Gamma_{TEV}$")
            savefig("../TevstrPlots/$(dirresults[i]).png")
            close()
            cd("..")
        end
    end
end

function calc_forces(surf::ThreeDSurfSimple)
    for i = 1:surf.nspan
        # First term in eqn (2.30) Ramesh et al. in coefficient form
        wi = 0
        for n = 1:surf.nspan
            nn = 2*n -2
            wi += real(nn)*surf.bcoeff[n]*sin(nn*surf.psi[i])/sin(surf.psi[i])
        end

        cnc = 2*pi*(surf.s2d[i].kinem.u*cos(surf.s2d[i].kinem.alpha)/surf.s2d[i].uref + surf.s2d[i].kinem.hdot*sin(surf.s2d[i].kinem.alpha)/surf.s2d[i].uref + wi*sin(surf.s2d[i].kinem.alpha))*(surf.s2d[i].a0[1] + surf.s2d[i].aterm[1]/2.)

        # Second term in eqn (2.30) Ramesh et al. in coefficient form
        cnnc = 2*pi*(3*surf.s2d[i].c*surf.s2d[i].a0dot[1]/(4*surf.s2d[i].uref) + surf.s2d[i].c*surf.s2d[i].adot[1]/(4*surf.s2d[i].uref) + surf.s2d[i].c*surf.s2d[i].adot[2]/(8*surf.s2d[i].uref))

        # Suction force given in eqn (2.31) Ramesh et al.
        cs = 2*pi*surf.s2d[i].a0[1]*surf.s2d[i].a0[1]

        #The components of normal force and moment from induced velocities are calulcated in dimensional units and nondimensionalized later
        nonl=0
        nonl_m=0
        for ib = 1:surf.s2d[i].ndiv-1
            nonl = nonl + (surf.s2d[i].uind[ib]*cos(surf.s2d[i].kinem.alpha) - surf.s2d[i].wind[ib]*sin(surf.s2d[i].kinem.alpha))*surf.s2d[i].bv[ib].s
            nonl_m = nonl_m + (surf.s2d[i].uind[ib]*cos(surf.s2d[i].kinem.alpha) - surf.s2d[i].wind[ib]*sin(surf.s2d[i].kinem.alpha))*surf.s2d[i].x[ib]*surf.s2d[i].bv[ib].s
        end
        nonl = nonl*2./(surf.s2d[i].uref*surf.s2d[i].uref*surf.s2d[i].c)
        nonl_m = nonl_m*2./(surf.s2d[i].uref*surf.s2d[i].uref*surf.s2d[i].c*surf.s2d[i].c)

        # Normal force coefficient
        cn = cnc + cnnc + nonl

        # Lift and drag coefficients
        surf.fc[i,1] = cn*cos(surf.s2d[i].kinem.alpha) + cs*sin(surf.s2d[i].kinem.alpha)
        surf.fc[i,2] = cn*sin(surf.s2d[i].kinem.alpha)-cs*cos(surf.s2d[i].kinem.alpha)

        #Pitching moment is clockwise or nose up positive
        surf.fc[i,3] = cn*surf.s2d[i].pvt - 2*pi*((surf.s2d[i].kinem.u*cos(surf.s2d[i].kinem.alpha)/surf.s2d[i].uref + surf.s2d[i].kinem.hdot*sin(surf.s2d[i].kinem.alpha)/surf.s2d[i].uref)*(surf.s2d[i].a0[1]/4. + surf.s2d[i].aterm[1]/4. - surf.s2d[i].aterm[2]/8.) + (surf.s2d[i].c/surf.s2d[i].uref)*(7.*surf.s2d[i].a0dot[1]/16. + 3.*surf.s2d[i].adot[1]/16. + surf.s2d[i].adot[2]/16. - surf.s2d[i].adot[3]/64.)) - nonl_m
    end
    cl3d = 0
    cd3d = 0
    cm3d = 0

    for i = 1:surf.nspan-1
        cl3d += (surf.fc[i,1] + surf.fc[i+1,1])*sin(0.5*(surf.psi[i]
        + surf.psi[i+1]))*(surf.psi[i+1] - surf.psi[i])/2

        cd3d +=  (surf.fc[i,2] + surf.fc[i+1,2])*sin(0.5*(surf.psi[i]
        + surf.psi[i+1]))*(surf.psi[i+1] - surf.psi[i])/2

        cm3d += (surf.fc[i,3] + surf.fc[i+1,3])*sin(0.5*(surf.psi[i]
        + surf.psi[i+1]))*(surf.psi[i+1] - surf.psi[i])/2

    end

    return cl3d, cd3d, cm3d
end
