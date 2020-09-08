
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

function makeVortPlots2D(tt = "all")
 # tt is a target time to make the plot for. It will choose the closest time available
    try
        cd("Step Files")
    catch
        println(" Error: No 'Step Files' directory. Make sure 'wflag' is enabled before running simulation.")
        return
    end

    dirvec = readdir()
    dirresults = map(x->(v = tryparse(Float64,x); typeof(v) == Nothing ? 0.0 : v),dirvec)
    if tt != "all"
        figure()
        ttidx = findmin(abs.(dirresults.-tt))[2]
    end
    #Determine axis limits
    dirmax = maximum(dirresults)
    cd("$(dirmax)")

    multsurfflag = 0
    if isfile("boundv-1") == true  #only 1 surface
        multsurfflag = 1
    end

    if multsurfflag == 0 #single surface

        tev = DelimitedFiles.readdlm("tev", '\t', Float64,  skipstart = 1)
        lev =   try
            DelimitedFiles.readdlm("lev", '\t', Float64,  skipstart = 1)
        catch
            zeros(0,3)
        end
        bv = DelimitedFiles.readdlm("boundv", '\t', Float64,  skipstart = 1)
        #tupletest = minimum([UnsteadyFlowSolvers.tev[:,2];UnsteadyFlowSolvers.lev[:,2];UnsteadyFlowSolvers.bv[:,2]])
        xmin = minimum([tev[:,2];lev[:,2];bv[:,2]])
        zmin = minimum([tev[:,3];lev[:,3];bv[:,3]])
        xmax = maximum([tev[:,2];lev[:,2]])
        zmax = maximum([bv[:,2];tev[:,3];lev[:,3];bv[:,3]])

        cd("../..")

        if tt == "all"
            dirvec = readdir()
            if "vortPlots" in dirvec
                rm("vortPlots", recursive=true)
            end
            mkdir("vortPlots")
        end
        cd("Step Files")

        for i=1:length(dirresults)
            if tt != "all" && ttidx != i
                continue
            end
            if dirresults[i] != 0
                dirstr="$(dirresults[i])"
                cd(dirstr)
                tev= DelimitedFiles.readdlm("tev", '\t', Float64,  skipstart = 1)
                lev =   try
                    DelimitedFiles.readdlm("lev", '\t', Float64,  skipstart = 1)
                catch
                    zeros(0,3)
                end
                bv= DelimitedFiles.readdlm("boundv", '\t', Float64,  skipstart = 1)

                viewVort2D(tev, lev, bv)
                axis([xmin-1, xmax+1, zmin-1, zmax+1])

                if tt == "all"
                    savefig("../../vortPlots/$(dirresults[i]).png")
                    close()
                end
                cd("..")
            end
        end
        if tt == "all"
            println(" All plots saved under the 'vortPlots' directory.")
        else
            pause(0.0001)
            show(block = false)
        end

    else #multiple surfaces

        tev, _ = DelimitedFiles.readdlm("tev", '\t', Float64,  skipstart = 1)
        lev =   try
            DelimitedFiles.readdlm("lev", '\t', Float64,  skipstart = 1)
        catch
            zeros(0,3)
        end

        nsurf = 0
        for i = 1:1000
            if isfile("boundv-$i") == true
                nsurf += 1

            end
        end
        bv, _ = DelimitedFiles.readdlm("boundv-1", '\t', Float64,  skipstart = 1)
        for i = 2:nsurf
            bvt, _ = DelimitedFiles.readdlm("boundv-$i", '\t', Float64,  skipstart = 1)
            bv = [bv;bvt;]
        end

        xmin = minimum([tev[:,2];lev[:,2];bv[:,2]])
        zmin = minimum([tev[:,3];lev[:,3];bv[:,3]])
        xmax = maximum([tev[:,2];lev[:,2]])
        zmax = maximum([bv[:,2];tev[:,3];lev[:,3];bv[:,3]])

        cd("..")

        if "vortPlots" in dirvec
            rm("vortPlots", recursive=true)
        end
        mkdir("vortPlots")
        for i=1:length(dirresults)
            if dirresults[i] != 0
                dirstr="$(dirresults[i])"
                cd(dirstr)
                tev, _ = DelimitedFiles.readdlm("tev", '\t', Float64,  skipstart = 1)
                lev =   try
                    DelimitedFiles.readdlm("lev", '\t', Float64,  skipstart = 1)
                catch
                    zeros(0,3)
                end
                bv, _ = DelimitedFiles.readdlm("boundv-1", '\t', Float64,  skipstart = 1)
                for i = 2:nsurf
                    bvt, _ = DelimitedFiles.readdlm("boundv-$i", '\t', Float64,  skipstart = 1)
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
    cd("..")
end

function makeForcePlots2D(mat = "")

    dirvec = readdir()
    if "forcePlots" in dirvec
        rm("forcePlots", recursive=true)
    end
    mkdir("forcePlots")

    if mat == ""
       mat = DelimitedFiles.readdlm("resultsSummary", Float64,  skipstart = 1)
    end

    if length(mat[1,:])==8 #only 1 surface
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
        xlabel("t*")
        ylabel(L"$\alpha$ (deg)")
        savefig("forcePlots/pitch.png")
        close()

        plot(t, h)
        xmin = t[1]
        xmax = t[end]
        zmin = minimum(h) - 0.1*abs(minimum(h))
        zmax = maximum(h) + 0.1*abs(maximum(h))
        axis([xmin, xmax, zmin, zmax])
        xlabel("t*")
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
        xlabel("t*")
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
        xlabel("t*")
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
        xlabel("t*")
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
            xlabel("t*")
            ylabel(L"$\alpha$ (deg)")
            savefig("forcePlots/pitch-$i.png")
            close()

            plot(t, h)
            xmin = t[1]
            xmax = t[end]
            zmin = minimum(h) - 0.1*abs(minimum(h))
            zmax = maximum(h) + 0.1*abs(maximum(h))
            axis([xmin, xmax, zmin, zmax])
            xlabel("t*")
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
            xlabel("t*")
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
            xlabel("t*")
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
            xlabel("t*")
            ylabel(L"$C_m$")
            savefig("forcePlots/cm-$i.png")
            close()
        end

    end

end

function checkConverge(k::Float64)

    mat = DelimitedFiles.readdlm("resultsSummary", '\t', Float64,  skipstart = 1)

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

    xlabel("t*")
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

    xlabel("t*")
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

    xlabel("t*")
    ylabel(L"$C_m$")
    savefig("forcePlots/cm-convergence.png")
    close()
end

function makeKinemClVortPlots2D(mat = "", tt = "all")
# tt is a target time to make the plot for. It will choose the closest time available

    if mat == ""
        mat = DelimitedFiles.readdlm("resultsSummary", '\t', Float64,  skipstart = 1)
    end

    t = mat[:,1]
    len = length(t)
    cl = mat[:,6]
    alpha = mat[:,2]*180/pi
    h = mat[:,3]
    u = mat[:,4]

    try
        cd("Step Files")
    catch
        println(" Error: No 'Step Files' directory. Make sure wflag is enabled.")
        return
    end

    dirvec = readdir()
    dirresults = map(x->(v = tryparse(Float64,x); typeof(v) == Nothing ? 0.0 : v),dirvec)
    if tt != "all"
        ttidx = findmin(abs.(dirresults.-tt))[2]
    end

    #Determine axis limits
    dirmax = maximum(dirresults)

    cd("$(dirmax)")

    multsurfflag = 0
    if isfile("boundv-1") == true  #only 1 surface
        error("this plot function is only written for single surface")
    end

    tev = DelimitedFiles.readdlm("tev", '\t', Float64,  skipstart = 1)
    lev =   try
        DelimitedFiles.readdlm("lev", '\t', Float64,  skipstart = 1)[1]
    catch
        zeros(0,3)
    end
    bv = DelimitedFiles.readdlm("boundv", '\t', Float64,  skipstart = 1)

    xminv = minimum([tev[:,2];lev[:,2];bv[:,2]])
    zminv = minimum([tev[:,3];lev[:,3];bv[:,3]])
    xmaxv = maximum([tev[:,2];lev[:,2]])
    zmaxv = maximum([bv[:,2];tev[:,3];lev[:,3];bv[:,3]])

    cd("../..")
    if tt == "all"
        dirvec = readdir()
        if "infoPlots" in dirvec
            rm("infoPlots", recursive=true)
        end
        mkdir("infoPlots")
    end

    cd("Step Files")
    for i=1:length(dirresults)
        if tt != "all" && ttidx != i
            continue
        end
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
            #xlabel("t*")
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
            #xlabel("t*")
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
            xlabel("t*")
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
            xlabel("t*")
            ylabel(L"$C_l$")

            subplot2grid((6, 2), (3, 0), rowspan=3, colspan=2)
            tev= DelimitedFiles.readdlm("tev", '\t', Float64,  skipstart = 1)
            lev =   try
                DelimitedFiles.readdlm("lev", '\t', Float64,  skipstart = 1)
            catch
                zeros(0,3)
            end
            bv = DelimitedFiles.readdlm("boundv", '\t', Float64,  skipstart = 1)
            viewVort2D(tev, lev, bv)
            axis([xminv-1, xmaxv+1, zminv-1, zmaxv+1])

            tight_layout()

            if tt == "all"
                savefig("../../infoPlots/$(dirresults[i]).png")
                close()
            end
            cd("..")
        end
    end
    cd("..")
    if tt == "all"
        println(" All plots saved under the 'infoPlots' directory.")
    else
        pause(0.0001)
        show(block = false)
    end
end

function makeKinemVelVortPlots2D(mat = "",tt = "all")

    if mat == ""
        mat = DelimitedFiles.readdlm("resultsSummary", '\t', Float64,  skipstart = 1)
    end

    t = mat[:,1]
    len = length(t)
    cl = mat[:,6]
    alpha = mat[:,2]*180/pi
    h = mat[:,3]
    u = mat[:,4]

    try
        cd("Step Files")
    catch
        println(" Error: No 'Step File' directory. Make sure wflag is enabled.")
        return
    end

    dirvec = readdir()
    dirresults = map(x->(v = tryparse(Float64,x); typeof(v) == Nothing ? 0.0 : v),dirvec)
    if tt != "all"
        ttidx = findmin(abs(dirresults-tt))[2]
    end


    #Determine axis limits
    dirmax = maximum(dirresults)

    cd("$(dirmax)")

    multsurfflag = 0
    if isfile("boundv-1") == true  #only 1 surface
        error("this plot function is only written for single surface")
    end

    tev = DelimitedFiles.readdlm("tev", '\t', Float64,  skipstart = 1)
    lev =   try
        DelimitedFiles.readdlm("lev", '\t', Float64,  skipstart = 1)
    catch
        zeros(0,3)
    end
    bv = DelimitedFiles.readdlm("boundv", '\t', Float64,  skipstart = 1)

    xminv = minimum([tev[:,2];lev[:,2];bv[:,2]])
    zminv = minimum([tev[:,3];lev[:,3];bv[:,3]])
    xmaxv = maximum([tev[:,2];lev[:,2]])
    zmaxv = maximum([bv[:,2];tev[:,3];lev[:,3];bv[:,3]])

    cd("../..")

    dirvec = readdir()
    if "infoplots" in dirvec
        rm("infoplots", recursive=true)
    end
    mkdir("infoplots")

    cd("Step Files")
    for i=1:length(dirresults)
        if tt != "all" && ttidx != i
            continue
        end
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
            #xlabel("t*")
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
            #xlabel("t*")
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
            xlabel("t*")
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
            xlabel("t*")
            ylabel(L"$C_l$")

            subplot2grid((6, 2), (3, 0), rowspan=3, colspan=2)
            tev = DelimitedFiles.readdlm("tev", '\t', Float64,  skipstart = 1)
            lev =   try
                DelimitedFiles.readdlm("lev", '\t', Float64,  skipstart = 1)
            catch
                zeros(0,3)
            end
            bv = DelimitedFiles.readdlm("boundv", '\t', Float64,  skipstart = 1)
            viewVort2D(tev, lev, bv)
            axis([xminv-1, xmaxv+1, zminv-1, zmaxv+1])

            tight_layout()

            savefig("../../infoplots/$(dirresults[i]).png")
            close()
            cd("..")
        end
    end
    cd("..")
end

function subPlotForces(mat, xVar :: Int64, xVarLabel :: String, saveName :: String)
    # Create Cl,Cd,Cm,LESP vs xVar plot
    # mat = [t , alpha , h , u , LESP , cl , cd , cm] for each timestep
    xVar = mat[:,xVar] ;
    cl = mat[:,6] ;
    cd = mat[:,7] ;
    cm = mat[:,8] ;
    LESP = mat[:,5] ;

    subplot(221)
    plot(xVar,cl) ;
    xlabel(xVarLabel)
    ylabel("Cl")

    subplot(222)
    plot(xVar,cd) ;
    xlabel(xVarLabel)
    ylabel("Cd")

    subplot(223)
    plot(xVar,cm) ;
    xlabel(xVarLabel)
    ylabel("Cm")

    subplot(224)
    plot(xVar,LESP) ;
    xlabel(xVarLabel)
    ylabel("LESP")

    tight_layout()
    savefig(saveName)
    close()
end

## UI Plots
function initPlot(A,C,t,a,h,u)
    # Create intitilization plots that shows airfoil and camber and motion parameters

    # Camber Plot
    if A != 0 # airfoil defined
        subplot(121,aspect = 1)
        plot(A[:,1],A[:,2],color = "black") # full airfoil
        plot(C[:,1],C[:,2], color = "blue", marker = "|") # Camber
        ylim(-.5,.5)
        xlabel("x")
        ylabel("z")
        title("Geometry")
    end

    # Motion plots
    if sum(a) != 0 # alphadef defined
        subplot(122)
        plot(t,h, color = "black") # height
        plot(t,u, color = "blue") # velocity
        plot(t,a, color = "red") # alpha
        xlabel("t*")
        legend(["Height","Velocity","Alpha [deg]"])
        title("Motion")
    end
    pause(0.0001)
    show(block = false)
end

function plotMat(mat,pltCnt,plotx,ploty,matIdx)
    # plot subplots (number determined by plotCnt) with an x-axis of plotx
    #   and a y axis of ploty. plotx/y are indices of the columns of mat.

    figure()
    # find subplot input dimensions
    if pltCnt <= 2 # plot on one row
        subRow = 1
        subCol = pltCnt
    else
        subRow = 2
        subCol = ceil(pltCnt/2)
    end
    subDim = subRow*100 + subCol*10

    for i = 1:pltCnt
        subplot(subDim + i)
        x = matIdx[lowercase(plotx)]
        y = matIdx[lowercase(ploty[i])]
        plot(mat[:,x],mat[:,y],"b-")
        xlabel(plotx)
        ylabel(ploty[i])
        title("$(ploty[i]) vs $plotx")
    end
    tight_layout()
    pause(0.0001)
    show(block = false)
end
