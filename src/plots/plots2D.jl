function viewVort2D(tev::Matrix{Float64}, lev::Matrix{Float64}, bv::Matrix{Float64})
    scatter(tev[:,2], tev[:,3], s=5, c=tev[:,1], edgecolors="none")
    sc = scatter(lev[:,2], lev[:,3], s=5, c=lev[:,1], edgecolors="none")
    sc2 = scatter(bv[:,2], bv[:,3], s=5, c=bv[:,1], edgecolors="none")
    plot(bv[:,2], bv[:,3], color = "black", linewidth=1.0)
end

function viewVortConnect2D(tev::Matrix{Float64}, lev::Matrix{Float64}, bv::Matrix{Float64})
    scatter(tev[:,2], tev[:,3], s=5, c=tev[:,1], edgecolors="none")
    sc = scatter(lev[:,2], lev[:,3], s=5, c=lev[:,1], edgecolors="none")
    sc2 = scatter(bv[:,2], bv[:,3], s=5, c=bv[:,1], edgecolors="none")

    plot(tev[:,2], tev[:,3], color = "red", linewidth=1.0)
    plot(lev[:,2], lev[:,3], color = "blue", linewidth=1.0)
    plot(bv[:,2], bv[:,3], color = "black", linewidth=1.0)
end

function makeVortPlots2D()
    dirvec = readdir()
    dirresults = map(x->(v = tryparse(Float64,x); typeof(v) == Nothing ? 0.0 : v),dirvec)
    #Determine axis limits
    dirmax = maximum(dirresults)
    cd("$(dirmax)")
    tev, _ = DelimitedFiles.readdlm("tev", '\t', Float64, header=true)
    lev =   try
                DelimitedFiles.readdlm("lev", '\t', Float64, header=true)[1]
            catch
                zeros(0,3)
            end
    bv, _ = DelimitedFiles.readdlm("boundv", '\t', Float64, header=true)

    xmin = minimum([tev[:,2];lev[:,2];bv[:,2];])
    zmin = minimum([tev[:,3];lev[:,3];bv[:,3];])
    xmax = maximum([tev[:,2];lev[:,2];])
    zmax = maximum([bv[:,2];tev[:,3];lev[:,3];bv[:,3];])
    cd("..")

    if "vortPlots" in dirvec
        rm("vortPlots", recursive=true)
    end
    mkdir("vortPlots")
    for i=1:length(dirresults)
        if dirresults[i] != 0
            dirstr="$(dirresults[i])"
            cd(dirstr)
            tev, _ = DelimitedFiles.readdlm("tev", '\t', Float64, header=true)
            lev =   try
                        DelimitedFiles.readdlm("lev", '\t', Float64, header=true)[1]
                    catch
                        zeros(0,3)
                    end
            bv, _ = DelimitedFiles.readdlm("boundv", '\t', Float64, header=true)

            viewVort2D(tev, lev, bv)
            axis([xmin-1, xmax+1, zmin-1, zmax+1])
            savefig("../vortPlots/$(dirresults[i]).png")
            close()
            cd("..")
        end
    end
end

function makeForcePlots()
    dirvec = readdir()
    if "forcePlots" in dirvec
        rm("forcePlots", recursive=true)
    end
    mkdir("forcePlots")

    mat, _ = DelimitedFiles.readdlm("resultsSummary", '\t', Float64, header=true)

    t = mat[:,1]
    alpha = mat[:,2]*180/pi
    h = mat[:,3]
    cl = mat[:,6]
    cd = mat[:,7]
    cm = mat[:,8]

    plot(t, alpha)
    len = length(t)
    xmin = t[1]
    xmax = t[end]
    zmin = minimum(alpha) - 0.1*abs(minimum(alpha))
    zmax = maximum(alpha) + 0.1*abs(maximum(alpha))
    axis([xmin, xmax, zmin, zmax])
    xlabel(L"$t^*$")
    ylabel(L"$\alpha$ (deg)")
    savefig("forcePlots/pitch.png")
    close()

    plot(t, h)
    xmin = t[1]
    xmax = t[end]
    zmin = minimum(h) - 0.1*abs(minimum(h))
    zmax = maximum(h) + 0.1*abs(maximum(h))
    axis([xmin, xmax, zmin, zmax])
    xlabel(L"$t^*$")
    ylabel(L"$h/c$")
    savefig("forcePlots/plunge.png")
    close()

    plot(t, cl)
    range = Int(round(0.05*len)):Int(round(0.95*len))
    xmin = t[1]
    xmax = t[end]
    zmin = minimum(cl[range]) - 0.1*abs(minimum(cl[range]))
    zmax = maximum(cl[range]) + 0.1*abs(maximum(cl[range]))
    axis([xmin, xmax, zmin, zmax])
    xlabel(L"$t^*$")
    ylabel(L"$C_l$")
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
    ylabel(L"$C_d$")
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
    ylabel(L"$C_m$")
    savefig("forcePlots/cm.png")
    close()
end
