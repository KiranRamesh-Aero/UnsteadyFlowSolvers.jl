function viewVort2D(tev::Matrix{Float64}, lev::Matrix{Float64}, bv::Matrix{Float64})
    scatter(tev[:,2], tev[:,3], s=5, c=tev[:,1], edgecolors="none")
    sc = scatter(lev[:,2], lev[:,3], s=5, c=lev[:,1], edgecolors="none")
    sc2 = scatter(bv[:,2], bv[:,3], s=5, c=bv[:,1], edgecolors="none")
    #plot(bv[:,2], bv[:,3], color = "black", linewidth=1.0)
end

#For quick views while debugging
function viewVort2D(tev::Vector{TwoDVort}, lev::Vector{TwoDVort}, bv::Vector{TwoDVort})
    tev = hcat(map(q->q.s, tev), map(q->q.x, tev), map(q->q.z, tev))
    lev = hcat(map(q->q.s, lev), map(q->q.x, lev), map(q->q.z, lev))
    bv = hcat(map(q->q.s, bv), map(q->q.x, bv), map(q->q.z, bv))

    scatter(tev[:,2], tev[:,3], s=5, c=tev[:,1], edgecolors="none")
    sc = scatter(lev[:,2], lev[:,3], s=5, c=lev[:,1], edgecolors="none")
    sc2 = scatter(bv[:,2], bv[:,3], s=5, c=bv[:,1], edgecolors="none")
    #plot(bv[:,2], bv[:,3], color = "black", linewidth=1.0)
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

    multsurfflag = 0
    if isfile("boundv-1") == true  #only 1 surface
        multsurfflag = 1
    end

    if multsurfflag == 0 #single surface

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

    else #multiple surfaces

        tev, _ = DelimitedFiles.readdlm("tev", '\t', Float64, header=true)
        lev =   try
            DelimitedFiles.readdlm("lev", '\t', Float64, header=true)[1]
        catch
            zeros(0,3)
        end

        nsurf = 0
        for i = 1:1000
            if isfile("boundv-$i") == true
                nsurf += 1

            end
        end
        bv, _ = DelimitedFiles.readdlm("boundv-1", '\t', Float64, header=true)
        for i = 2:nsurf
            bvt, _ = DelimitedFiles.readdlm("boundv-$i", '\t', Float64, header=true)
            bv = [bv;bvt;]
        end

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
                bv, _ = DelimitedFiles.readdlm("boundv-1", '\t', Float64, header=true)
                for i = 2:nsurf
                    bvt, _ = DelimitedFiles.readdlm("boundv-$i", '\t', Float64, header=true)
                    bv = [bv;bvt;]
                end
                viewVort2D(tev, lev, bv)
                axis([xmin-1, xmax+1, zmin-1, zmax+1])
                savefig("../vortPlots/$(dirresults[i]).png")
                close()
                cd("..")
            end
        end
    end
end

function makeForcePlots2D()
    dirvec = readdir()
    if "forcePlots" in dirvec
        rm("forcePlots", recursive=true)
    end
    mkdir("forcePlots")

    mat, _ = DelimitedFiles.readdlm("resultsSummary", '\t', Float64, header=true)

    if length(mat[1,:]) == 8 #only 1 surface
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

    else #multiple surfaces

        nsurf = (length(mat[1,:]) - 1)/7
        t = mat[:,1]

        for i = 1:nsurf
            alpha = mat[:,Int((i-1)*7+2)]*180/pi
            h = mat[:,Int((i-1)*7+3)]
            cl = mat[:,Int((i-1)*7+6)]
            cd = mat[:,Int((i-1)*7+7)]
            cm = mat[:,Int((i-1)*7+8)]

            plot(t, alpha)
            len = length(t)
            xmin = t[1]
            xmax = t[end]
            zmin = minimum(alpha) - 0.1*abs(minimum(alpha))
            zmax = maximum(alpha) + 0.1*abs(maximum(alpha))
            axis([xmin, xmax, zmin, zmax])
            xlabel(L"$t^*$")
            ylabel(L"$\alpha$ (deg)")
            savefig("forcePlots/pitch-$i.png")
            close()

            plot(t, h)
            xmin = t[1]
            xmax = t[end]
            zmin = minimum(h) - 0.1*abs(minimum(h))
            zmax = maximum(h) + 0.1*abs(maximum(h))
            axis([xmin, xmax, zmin, zmax])
            xlabel(L"$t^*$")
            ylabel(L"$h/c$")
            savefig("forcePlots/plunge-$i.png")
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
            savefig("forcePlots/cl-$i.png")
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
            savefig("forcePlots/cd-$i.png")
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
            savefig("forcePlots/cm-$i.png")
            close()
        end

    end

end

function checkConverge(k::Float64)

    mat, _ = DelimitedFiles.readdlm("resultsSummary", '\t', Float64, header=true)

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

    #Lift convergence
    figure
    for i = 1:ncyc
        start_t = real(i-1)*T
        end_t = real(i)*T
        start_ind = argmin(abs.(mat[:,1] .- start_t))
        end_ind = argmin(abs.(mat[:,1] .- end_t))
        plot((mat[start_ind:end_ind,1] .- start_t)/T, mat[start_ind:end_ind,6])
    end
    xmin = 0.
    xmax = 1.
    zmin = minimum(mat[:,6])
    zmax = maximum(mat[:,6])
    axis([xmin, xmax, zmin, zmax])

    xlabel(L"$t^*$")
    ylabel(L"$C_l$")
    savefig("forcePlots/cl-convergence.png")
    close()

    #Drag convergence
    figure
    for i = 1:ncyc
        start_t = real(i-1)*T
        end_t = real(i)*T
        start_ind = argmin(abs.(mat[:,1] .- start_t))
        end_ind = argmin(abs.(mat[:,1] .- end_t))
        plot((mat[start_ind:end_ind,1] .- start_t)/T, mat[start_ind:end_ind,7])
    end
    xmin = 0.
    xmax = 1.
    zmin = minimum(mat[:,7])
    zmax = maximum(mat[:,7])
    axis([xmin, xmax, zmin, zmax])

    xlabel(L"$t^*$")
    ylabel(L"$C_d$")
    savefig("forcePlots/cd-convergence.png")
    close()

    #Pitching moment convergence
    figure
    for i = 1:ncyc
        start_t = real(i-1)*T
        end_t = real(i)*T
        start_ind = argmin(abs.(mat[:,1] .- start_t))
        end_ind = argmin(abs.(mat[:,1] .- end_t))
        plot((mat[start_ind:end_ind,1] .- start_t)/T, mat[start_ind:end_ind,8])
    end
    xmin = 0.
    xmax = 1.
    zmin = minimum(mat[:,8])
    zmax = maximum(mat[:,8])
    axis([xmin, xmax, zmin, zmax])

    xlabel(L"$t^*$")
    ylabel(L"$C_m$")
    savefig("forcePlots/cm-convergence.png")
    close()
end

function makeKinemClVortPlots2D()

    mat, _ = DelimitedFiles.readdlm("resultsSummary", '\t', Float64, header=true)

    t = mat[:,1]
    len = length(t)
    cl = mat[:,6]
    alpha = mat[:,2]*180/pi
    h = mat[:,3]
    u = mat[:,4]

    dirvec = readdir()
    dirresults = map(x->(v = tryparse(Float64,x); typeof(v) == Nothing ? 0.0 : v),dirvec)

    #Determine axis limits
    dirmax = maximum(dirresults)

    cd("$(dirmax)")

    multsurfflag = 0
    if isfile("boundv-1") == true  #only 1 surface
        error("this plot function is only written for single surface")
    end

    tev, _ = DelimitedFiles.readdlm("tev", '\t', Float64, header=true)
    lev =   try
        DelimitedFiles.readdlm("lev", '\t', Float64, header=true)[1]
    catch
        zeros(0,3)
    end
    bv, _ = DelimitedFiles.readdlm("boundv", '\t', Float64, header=true)

    xminv = minimum([tev[:,2];lev[:,2];bv[:,2];])
    zminv = minimum([tev[:,3];lev[:,3];bv[:,3];])
    xmaxv = maximum([tev[:,2];lev[:,2];])
    zmaxv = maximum([bv[:,2];tev[:,3];lev[:,3];bv[:,3];])

    cd("..")

    if "infoPlots" in dirvec
        rm("infoPlots", recursive=true)
    end
    mkdir("infoPlots")

    for i=1:length(dirresults)
        if dirresults[i] != 0
            dirstr="$(dirresults[i])"
            cd(dirstr)

            t_cur = dirresults[i]

            figure(figsize=(8,8))

            subplot2grid((6, 2), (0, 0))
            plot(t, alpha)
            plot([t_cur; t_cur], [-10000; 10000], "k-")
            range = Int(round(0.05*len)):Int(round(0.95*len))
            xmin = t[1]
            xmax = t[end]
            zmin = minimum(alpha[range]) - 0.1*abs(minimum(alpha[range]))
            zmax = maximum(alpha[range]) + 0.1*abs(maximum(alpha[range]))
            axis([xmin, xmax, zmin, zmax])
            #xlabel(L"$t^*$")
            xticks([])
            ylabel(L"$\alpha (deg)$")

            subplot2grid((6, 2), (1, 0))
            plot(t, h)
            plot([t_cur; t_cur], [-10000; 10000], "k-")
            range = Int(round(0.05*len)):Int(round(0.95*len))
            xmin = t[1]
            xmax = t[end]
            zmin = minimum(h[range]) - 0.1*abs(minimum(h[range]))
            zmax = maximum(h[range]) + 0.1*abs(maximum(h[range]))
            axis([xmin, xmax, zmin, zmax])
            #xlabel(L"$t^*$")
            xticks([])
            ylabel(L"$h$")

            subplot2grid((6, 2), (2, 0))
            plot(t, u)
            plot([t_cur; t_cur], [-10000; 10000], "k-")
            range = Int(round(0.05*len)):Int(round(0.95*len))
            xmin = t[1]
            xmax = t[end]
            zmin = minimum(u[range]) - 0.1*abs(minimum(u[range]))
            zmax = maximum(u[range]) + 0.1*abs(maximum(u[range]))
            axis([xmin, xmax, zmin, zmax])
            xlabel(L"$t^*$")
            ylabel(L"$u$")

            subplot2grid((6, 2), (0, 1), rowspan=3)
            plot(t, cl)
            plot([t_cur; t_cur], [-10000; 10000], "k-")
            range = Int(round(0.05*len)):Int(round(0.95*len))
            xmin = t[1]
            xmax = t[end]
            zmin = minimum(cl[range]) - 0.1*abs(minimum(cl[range]))
            zmax = maximum(cl[range]) + 0.1*abs(maximum(cl[range]))
            axis([xmin, xmax, zmin, zmax])
            xlabel(L"$t^*$")
            ylabel(L"$C_l$")

            subplot2grid((6, 2), (3, 0), rowspan=3, colspan=2)
            tev, _ = DelimitedFiles.readdlm("tev", '\t', Float64, header=true)
            lev =   try
                DelimitedFiles.readdlm("lev", '\t', Float64, header=true)[1]
            catch
                zeros(0,3)
            end
            bv, _ = DelimitedFiles.readdlm("boundv", '\t', Float64, header=true)
            viewVort2D(tev, lev, bv)
            axis([xminv-1, xmaxv+1, zminv-1, zmaxv+1])

            tight_layout()

            savefig("../infoPlots/$(dirresults[i]).png")
            close()
            cd("..")
        end
    end
end

function makeKinemVelVortPlots2D()

    mat, _ = DelimitedFiles.readdlm("resultsSummary", '\t', Float64, header=true)

    t = mat[:,1]
    len = length(t)
    cl = mat[:,6]
    alpha = mat[:,2]*180/pi
    h = mat[:,3]
    u = mat[:,4]

    dirvec = readdir()
    dirresults = map(x->(v = tryparse(Float64,x); typeof(v) == Nothing ? 0.0 : v),dirvec)

    #Determine axis limits
    dirmax = maximum(dirresults)

    cd("$(dirmax)")

    multsurfflag = 0
    if isfile("boundv-1") == true  #only 1 surface
        error("this plot function is only written for single surface")
    end

    tev, _ = DelimitedFiles.readdlm("tev", '\t', Float64, header=true)
    lev =   try
        DelimitedFiles.readdlm("lev", '\t', Float64, header=true)[1]
    catch
        zeros(0,3)
    end
    bv, _ = DelimitedFiles.readdlm("boundv", '\t', Float64, header=true)

    xminv = minimum([tev[:,2];lev[:,2];bv[:,2];])
    zminv = minimum([tev[:,3];lev[:,3];bv[:,3];])
    xmaxv = maximum([tev[:,2];lev[:,2];])
    zmaxv = maximum([bv[:,2];tev[:,3];lev[:,3];bv[:,3];])

    cd("..")

    if "infoPlots" in dirvec
        rm("infoPlots", recursive=true)
    end
    mkdir("infoPlots")

    for i=1:length(dirresults)
        if dirresults[i] != 0
            dirstr="$(dirresults[i])"
            cd(dirstr)

            t_cur = dirresults[i]

            figure(figsize=(8,8))

            subplot2grid((6, 2), (0, 0))
            plot(t, alpha)
            plot([t_cur; t_cur], [-10000; 10000], "k-")
            range = Int(round(0.05*len)):Int(round(0.95*len))
            xmin = t[1]
            xmax = t[end]
            zmin = minimum(alpha[range]) - 0.1*abs(minimum(alpha[range]))
            zmax = maximum(alpha[range]) + 0.1*abs(maximum(alpha[range]))
            axis([xmin, xmax, zmin, zmax])
            #xlabel(L"$t^*$")
            xticks([])
            ylabel(L"$\alpha (deg)$")

            subplot2grid((6, 2), (1, 0))
            plot(t, h)
            plot([t_cur; t_cur], [-10000; 10000], "k-")
            range = Int(round(0.05*len)):Int(round(0.95*len))
            xmin = t[1]
            xmax = t[end]
            zmin = minimum(h[range]) - 0.1*abs(minimum(h[range]))
            zmax = maximum(h[range]) + 0.1*abs(maximum(h[range]))
            axis([xmin, xmax, zmin, zmax])
            #xlabel(L"$t^*$")
            xticks([])
            ylabel(L"$h$")

            subplot2grid((6, 2), (2, 0))
            plot(t, u)
            plot([t_cur; t_cur], [-10000; 10000], "k-")
            range = Int(round(0.05*len)):Int(round(0.95*len))
            xmin = t[1]
            xmax = t[end]
            zmin = minimum(u[range]) - 0.1*abs(minimum(u[range]))
            zmax = maximum(u[range]) + 0.1*abs(maximum(u[range]))
            axis([xmin, xmax, zmin, zmax])
            xlabel(L"$t^*$")
            ylabel(L"$u$")

            subplot2grid((6, 2), (0, 1), rowspan=3)
            plot(t, cl)
            plot([t_cur; t_cur], [-10000; 10000], "k-")
            range = Int(round(0.05*len)):Int(round(0.95*len))
            xmin = t[1]
            xmax = t[end]
            zmin = minimum(cl[range]) - 0.1*abs(minimum(cl[range]))
            zmax = maximum(cl[range]) + 0.1*abs(maximum(cl[range]))
            axis([xmin, xmax, zmin, zmax])
            xlabel(L"$t^*$")
            ylabel(L"$C_l$")

            subplot2grid((6, 2), (3, 0), rowspan=3, colspan=2)
            tev, _ = DelimitedFiles.readdlm("tev", '\t', Float64, header=true)
            lev =   try
                DelimitedFiles.readdlm("lev", '\t', Float64, header=true)[1]
            catch
                zeros(0,3)
            end
            bv, _ = DelimitedFiles.readdlm("boundv", '\t', Float64, header=true)
            viewVort2D(tev, lev, bv)
            axis([xminv-1, xmaxv+1, zminv-1, zmaxv+1])

            tight_layout()

            savefig("../infoPlots/$(dirresults[i]).png")
            close()
            cd("..")
        end
    end
end

function makeVelContourPlots2D()

    dirvec = readdir()
    dirresults = map(x->(v = tryparse(Float64,x); typeof(v) == Nothing ? 0.0 : v),dirvec)
    #Determine axis limits
    dirmax = maximum(dirresults)

    #Add code for multiple surfaces later

    field = jldopen("$(dirmax).jld", "r") do file
        read(file, "field")
    end
    surf = jldopen("$(dirmax).jld", "r") do file
        read(file, "surf")
    end

    stev = field.tev
    slev = field.lev
    sbv = surf.bv

    tev = zeros(length(stev), 3)
    for i = 1:length(stev)
        tev[i,:] = [stev[i].s stev[i].x stev[i].z]
    end
    lev = zeros(length(slev), 3)
    for i = 1:length(slev)
        lev[i,:] = [slev[i].s slev[i].x slev[i].z]
    end
    lev = zeros(length(slev), 3)
    for i = 1:length(slev)
        lev[i,:] = [slev[i].s slev[i].x slev[i].z]
    end
    bv = zeros(length(sbv), 3)
    for i = 1:length(sbv)
        bv[i,:] = [sbv[i].s sbv[i].x sbv[i].z]
    end


    xmin = maximum([minimum([tev[:,2];lev[:,2];bv[:,2];]) .- surf.c; surf.bnd_x[1] - surf.c])
    zmin = maximum([minimum([tev[:,3];lev[:,3];bv[:,3];]) .- surf.c; surf.bnd_z[1] - 2*surf.c])
    xmax = minimum([maximum([tev[:,2];lev[:,2];]) + surf.c; surf.bnd_x[surf.ndiv] + 5*surf.c])
    zmax = minimum([maximum([bv[:,2];tev[:,3];lev[:,3];bv[:,3];]) + surf.c; surf.bnd_z[1] + 2*surf.c])

    if "velContourPlots" in dirvec
        rm("velContourPlots", recursive=true)
    end
    mkdir("velContourPlots")
    for i=1:length(dirresults)
        if dirresults[i] != 0
            dirstr="$(dirresults[i])"

            field = jldopen("$(dirstr).jld", "r") do file
                read(file, "field")
            end
            surf = jldopen("$(dirstr).jld", "r") do file
                read(file, "surf")
            end

            #Transform body to inertial axis
            #Translation
            xb = surf.bnd_x[:]
            zb = surf.bnd_z[:]
            pvt_ind = argmin(abs.(surf.x .- surf.pvt*surf.c))
            dist = surf.pvt*surf.c - surf.bnd_x[pvt_ind]
            xt1 = xb .+ dist
            zt1 = zb

            #Rotation
            theta = surf.kinem.alpha
            trans = [cos(theta) -sin(theta); sin(theta) cos(theta)]
            xt = zeros(surf.ndiv)
            zt = zeros(surf.ndiv)
            for i = 1:surf.ndiv
                xt[i], zt[i] = trans*[xt1[i]; zt1[i]]
            end

            #Create grid of points around body
            #Below body
            dx = surf.c/surf.ndiv
            dz = dx
            n = Int(ceil((zmax-zmin)/dz/2))
            println(n)

            xmint = xmin + dist
            xmaxt = xmax + dist

            x = [collect(xmint:dx:xt[1]); collect(xt[2:surf.ndiv-1]); collect(xt[surf.ndiv]:dx:xmaxt)];
            l1 = length(collect(xmint:dx:xt[1]))
            l2 = length(collect(xt[surf.ndiv]:dx:xmaxt))
            zu = [ones(l1)*dz; collect(zt[2:surf.ndiv-1]).+dz; ones(l2)*dz]
            zl = [-ones(l1)*dz; collect(zt[2:surf.ndiv-1]).-dz; -ones(l2)*dz]


            xmat = repeat(x', 2*n, 1)
            zmat = zeros(size(xmat))
            zmat[n+1,:] = zu[:]
            for i = 1:n-1
                zmat[n+1+i,:] = zu[:] .+ i*dz
            end
            zmat[n,:] = zl[:]
            for i = 1:n-1
                zmat[n-i,:] = zu[:] .- i*dz
            end

            #Transform grid back to body frame
            xtr = zeros(size(xmat))
            ztr = zeros(size(xmat))
            #Rotation
            trans = [cos(theta) sin(theta); -sin(theta) cos(theta)]
            for i = 1:size(xmat,1)
                for j = 1:size(xmat,2)
                    xr, zr = trans*[xmat[i,j]; zmat[i,j]]
                    xtr[i,j] = xr - dist
                    ztr[i,j] = zr
                end
            end

            uvel_i = zeros(size(xmat))
            wvel_i = zeros(size(xmat))
            vel_i = zeros(size(xmat))
            #Find velocity field over grid
            for i = 1:size(xtr,1)
                for j = 1:size(xtr,2)
                    utemp, wtemp = UnsteadyFlowSolvers.ind_vel([field.tev; field.lev; surf.bv], xtr[i,j], ztr[i,j])
                    uvel_i[i,j] += utemp[1]
                    wvel_i[i,j] += wtemp[1]
                    uvel_i[i,j] += field.u[1]
                    wvel_i[i,j] += field.w[1]
                    vel_i[i,j] = sqrt(uvel_i[i,j]^2 + wvel_i[i,j]^2)
                end
            end


            viewVort2D(tev, lev, bv)
            axis([xmin-1, xmax+1, zmin-1, zmax+1])
            savefig("velContourPlots/$(dirresults[i]).png")
            close()
        end
    end
end


