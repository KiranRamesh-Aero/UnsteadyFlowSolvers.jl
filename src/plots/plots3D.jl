function makeForcePlots3D()
    
    dirvec = readdir()
    if "forcePlots" in dirvec
        rm("forcePlots", recursive=true)
    end
    mkdir("forcePlots")

    mat, _ = DelimitedFiles.readdlm("resultsSummary", '\t', Float64, header=true)

    nspan = (length(mat[1,:]) - 4)/8
    
    t = mat[:,1]
    len = length(t)
    cl = mat[:,2]
    cd = mat[:,3]
    cm = mat[:,4]

    #Plot wing coefficients
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

    #Co-plot kinematics for the various sections
    alpha_all = Float64[]
    h_all = Float64[]
    for i = 1:nspan
        alpha = mat[:,Int((i-1)*8+5)]*180/pi
        alpha_all = [alpha_all; alpha]
        h = mat[:,Int((i-1)*8+6)]
        h_all = [h_all; h]
    end
    for i = 1:nspan
        alpha = mat[:,Int((i-1)*8+5)]*180/pi
        plot(t, alpha)
    end
    xmin = t[1]
    xmax = t[end]
    zmin = minimum(alpha_all) - 0.1*abs(minimum(alpha_all))
    zmax = maximum(alpha_all) + 0.1*abs(maximum(alpha_all))
    axis([xmin, xmax, zmin, zmax])
    xlabel(L"$t^*$")
    ylabel(L"$\alpha$ (deg)")
    savefig("forcePlots/pitch-spanwise.png")
    close()

    for i = 1:nspan
        h = mat[:,Int((i-1)*8+6)]
        plot(t, h)
    end
    xmin = t[1]
    xmax = t[end]
    zmin = minimum(h_all) - 0.1*abs(minimum(h_all))
    zmax = maximum(h_all) + 0.1*abs(maximum(h_all))
    axis([xmin, xmax, zmin, zmax])
    xlabel(L"$t^*$")
    ylabel(L"$h/c$")
    savefig("forcePlots/plunge-spanwise.png")
    close()

    #Coplot spanwise force coefficients
    cl_all = Float64[]
    cd_all = Float64[]
    cm_all = Float64[]
    for i = 1:nspan
        cl = mat[range,Int((i-1)*8+9)]./mat[range,2]
        cl_all = [cl_all; cl]
        cd = mat[range,Int((i-1)*8+10)]./mat[range,3]
        cd_all = [cd_all; cd]
        cm = mat[range,Int((i-1)*8+11)]./mat[range,4]
        cm_all = [cm_all; cm]
    end

    for i = 1:nspan
        cl = mat[:,Int((i-1)*8+9)]./mat[:,2]
        plot(t, cl)
    end
    xmin = t[1]
    xmax = t[end]
    zmin = minimum(cl_all[:]) - 0.1*abs(minimum(cl_all[:]))
    zmax = maximum(cl_all[:]) + 0.1*abs(maximum(cl_all[:]))
    axis([xmin, xmax, zmin, zmax])
    xlabel(L"$t^*$")
    ylabel(L"$C_l c/C_L c_{ref}$")
    savefig("forcePlots/cl-spanwise.png")
    close()
    
    for i = 1:nspan
        cd = mat[:,Int((i-1)*8+10)]./mat[:,3]
        plot(t, cd)
    end
    xmin = t[1]
    xmax = t[end]
    zmin = minimum(cd_all[:]) - 0.1*abs(minimum(cd_all[:]))
    zmax = maximum(cd_all[:]) + 0.1*abs(maximum(cd_all[:]))
    axis([xmin, xmax, zmin, zmax])
    xlabel(L"$t^*$")
    ylabel(L"$C_d c/C_D c_{ref}$")
    savefig("forcePlots/cd-spanwise.png")
    close()


    for i = 1:nspan
        cm = mat[:,Int((i-1)*8+11)]./mat[:,4]
        plot(t, cm)
    end
    xmin = t[1]
    xmax = t[end]
    zmin = minimum(cm_all[:]) - 0.1*abs(minimum(cm_all[:]))
    zmax = maximum(cm_all[:]) + 0.1*abs(maximum(cm_all[:]))
    axis([xmin, xmax, zmin, zmax])
    xlabel(L"$t^*$")
    ylabel(L"$C_m c^2/C_M c_{ref}^2$")
    savefig("forcePlots/cm-spanwise.png")
    close()
    
    #Coplot spanwise a0 variation
    a0_all = Float64[]
    a03d_all = Float64[]
    for i = 1:nspan
        a0 = mat[:,Int((i-1)*8+8)]
        a0_all = [a0_all; a0]
        a03d = mat[:,Int((i-1)*8+12)]
        a03d_all = [a03d_all; a03d]
    end
    for i = 1:nspan
        a0 = mat[:,Int((i-1)*8+8)]
        plot(t, a0)
    end
    xmin = t[1]
    xmax = t[end]
zmin = minimum(a0_all[:]) - 0.1*abs(minimum(a0_all[:]))
zmax = maximum(a0_all[:]) + 0.1*abs(maximum(a0_all[:]))
axis([xmin, xmax, zmin, zmax])
xlabel(L"$t^*$")
ylabel(L"$A_0$")
savefig("forcePlots/a0-spanwise.png")
close()

for i = 1:nspan
    a03d = mat[:,Int((i-1)*8+12)]
    plot(t, a03d)
end
xmin = t[1]
xmax = t[end]
zmin = minimum(a03d_all[:]) - 0.1*abs(minimum(a03d_all[:]))
zmax = maximum(a03d_all[:]) + 0.1*abs(maximum(a03d_all[:]))
axis([xmin, xmax, zmin, zmax])
xlabel(L"$t^*$")
ylabel(L"$A_0$")
savefig("forcePlots/a03d-spanwise.png")
close()
    
end

function makeVortPlots3D()

    dirvec = readdir()
    dirresults = map(x->(v = tryparse(Float64,x); typeof(v) == Nothing ? 0.0 : v),dirvec)
    #Determine axis limits
    dirmax = maximum(dirresults)
    cd("$(dirmax)")
    
    nspan = 0
    for i = 1:1000
        if isfile("boundv-$i") == true
            nspan += 1
        end
    end
    tip_ind = 1
    root_ind = nspan
    mid_ind= floor(Int, nspan/2)

    tev, _ = DelimitedFiles.readdlm("tev-1", '\t', Float64, header=true)
    if length(tev[:,1]) > 1
        tev = tev[2:end,:]
    end
    lev =   try
        DelimitedFiles.readdlm("lev-1", '\t', Float64, header=true)[1]
    catch
        zeros(0,3)
    end
    if length(lev[:,1]) > 1
        lev = lev[2:end,:]
    end
    bv, _ = DelimitedFiles.readdlm("boundv-1", '\t', Float64, header=true)    
    for i in [mid_ind; root_ind]
        tevt, _ = DelimitedFiles.readdlm("tev-$i", '\t', Float64, header=true)
        if length(tevt[:,1]) > 1
            tevt = tevt[2:end,:]
        end
        tev = [tev; tevt]
        levt = try
            DelimitedFiles.readdlm("lev-$i", '\t', Float64, header=true)[1]
        catch
            zeros(0,3)
        end
        if length(levt[:,1]) > 1
            levt = levt[2:end,:]
        end
        lev = [lev; levt]
        bvt, _ = DelimitedFiles.readdlm("boundv-$i", '\t', Float64, header=true)
        bv = [bv; bvt]
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

            figure()

            tev, _ = DelimitedFiles.readdlm("tev-$root_ind", '\t', Float64, header=true)
            lev =   try
                DelimitedFiles.readdlm("lev-$root_ind", '\t', Float64, header=true)[1]
            catch
                zeros(0,3)
            end
            bv, _ = DelimitedFiles.readdlm("boundv-$root_ind", '\t', Float64, header=true)
            subplot(311)
            viewVort2D(tev, lev, bv)
            axis([xmin-1, xmax+1, zmin-1, zmax+1])
            ylabel("Root")
            
            tev, _ = DelimitedFiles.readdlm("tev-$mid_ind", '\t', Float64, header=true)
            lev =   try
                DelimitedFiles.readdlm("lev-$mid_ind", '\t', Float64, header=true)[1]
            catch
                zeros(0,3)
            end
            bv, _ = DelimitedFiles.readdlm("boundv-$mid_ind", '\t', Float64, header=true)
            subplot(312)
            viewVort2D(tev, lev, bv)
            axis([xmin-1, xmax+1, zmin-1, zmax+1])
            ylabel("Mid")

            tev, _ = DelimitedFiles.readdlm("tev-$tip_ind", '\t', Float64, header=true)
            lev =   try
                DelimitedFiles.readdlm("lev-$tip_ind", '\t', Float64, header=true)[1]
            catch
                zeros(0,3)
            end
            bv, _ = DelimitedFiles.readdlm("boundv-$tip_ind", '\t', Float64, header=true)
            subplot(313)
            viewVort2D(tev, lev, bv)
            axis([xmin-1, xmax+1, zmin-1, zmax+1])
            ylabel("Tip")
            
            savefig("../vortPlots/$(dirresults[i]).png")
            close()
            cd("..")
        end
    end
end

function makeInfoPlots3D()

    mat, _ = DelimitedFiles.readdlm("resultsSummary", '\t', Float64, header=true)
    nspan = Int((length(mat[1,:]) - 4)/8)

    t = mat[:,1]
    len = length(t)
    cl = mat[:,2]
    cdr = mat[:,3]
    cm = mat[:,4]
    
    dirvec = readdir()
    dirresults = map(x->(v = tryparse(Float64,x); typeof(v) == Nothing ? 0.0 : v),dirvec)
    
    #Determine axis limits
    dirmax = maximum(dirresults)

    cd("$(dirmax)")

    tip_ind = 1
    root_ind = nspan
    mid_ind= floor(Int, nspan/2)

    #Axis limits for vorticity plots based on last time step
    tev, _ = DelimitedFiles.readdlm("tev-1", '\t', Float64, header=true)
    if length(tev[:,1]) > 1
        tev = tev[2:end,:]
    end
    lev =   try
        DelimitedFiles.readdlm("lev-1", '\t', Float64, header=true)[1]
    catch
        zeros(0,3)
    end
    lev, _ = DelimitedFiles.readdlm("lev-1", '\t', Float64, header=true)
    if length(lev[:,1]) > 1
        lev = lev[2:end,:]
    end
    bv, _ = DelimitedFiles.readdlm("boundv-1", '\t', Float64, header=true)
    for i in [mid_ind; root_ind]
        tevt, _ = DelimitedFiles.readdlm("tev-$i", '\t', Float64, header=true)
        if length(tevt[:,1]) > 1
            tevt = tevt[2:end,:]
        end
        tev = [tev; tevt]
        levt = try
            DelimitedFiles.readdlm("lev-$i", '\t', Float64, header=true)[1]
        catch
            zeros(0,3)
        end
        if length(levt[:,1]) > 1
            levt = levt[2:end,:]
        end
        lev = [lev; levt]
        bvt, _ = DelimitedFiles.readdlm("boundv-$i", '\t', Float64, header=true)
        bv = [bv; bvt]
    end
    cd("..")
    
    xminv = minimum([tev[:,2];lev[:,2];bv[:,2];])
    zminv = minimum([tev[:,3];lev[:,3];bv[:,3];])
    xmaxv = maximum([tev[:,2];lev[:,2];])
    zmaxv = maximum([bv[:,2];tev[:,3];lev[:,3];bv[:,3];])

    
    #Axis limits for spanwise variations based on all the time steps
    vorstr = Float64[]
    a0str = Float64[]
    gamstr = Float64[]
    for i=1:length(dirresults)
        if dirresults[i] != 0
            dirstr="$(dirresults[i])"
            cd(dirstr)
            
            matind = argmin(abs.(t .- dirresults[i]))
            
            spanvar, _  = DelimitedFiles.readdlm("spanwise-var", '\t', Float64, header=true)
            vorstr = [vorstr; spanvar[:,2]; spanvar[:,3]]
            a0str = [a0str; spanvar[:,5]; spanvar[:,6]; spanvar[:,7]]
            cl_load = mat[matind, ([1:nspan;].-1).*8 .+ 9]
            gamstr = [gamstr; spanvar[:,4]/spanvar[end,4]; cl_load]
            cd("..")
        end
    end
    vormin = minimum(vorstr) - 0.1*abs(minimum(vorstr))
    a0min = minimum(a0str) - 0.1*abs(minimum(a0str))
    gammin = minimum(gamstr) - 0.1*abs(minimum(gamstr))

    vormax = maximum(vorstr) + 0.1*abs(maximum(vorstr))
    a0max = maximum(a0str) + 0.1*abs(maximum(a0str))
    gammax = maximum(gamstr) + 0.1*abs(maximum(gamstr))
    
    if "infoPlots" in dirvec
        rm("infoPlots", recursive=true)
    end
    mkdir("infoPlots")
    
    for i=1:length(dirresults)
        if dirresults[i] != 0
            dirstr="$(dirresults[i])"
            cd(dirstr)
            
            figure(figsize=(12,8))
            
            #Kinematics
            alpha_all = Float64[]
            h_all = Float64[]
            for i = 1:nspan
                alpha = mat[:,Int((i-1)*8+5)]*180/pi
                alpha_all = [alpha_all; alpha]
                h = mat[:,Int((i-1)*8+6)]
                h_all = [h_all; h]
            end
            
            subplot2grid((4, 5), (0, 0))
            for i = 1:nspan
                alpha = mat[:,Int((i-1)*8+5)]*180/pi
                plot(t, alpha)
            end
            xmin = t[1]
            xmax = t[end]
            zmin = minimum(alpha_all) - 0.1*abs(minimum(alpha_all))
            zmax = maximum(alpha_all) + 0.1*abs(maximum(alpha_all))
            axis([xmin, xmax, zmin, zmax])
            xlabel(L"$t^*$")
            ylabel(L"$\alpha$ (deg)")
            
            subplot2grid((4, 5), (0, 1))
            for i = 1:nspan
                h = mat[:,Int((i-1)*8+6)]
                plot(t, h)
            end
            xmin = t[1]
            xmax = t[end]
            zmin = minimum(h_all) - 0.1*abs(minimum(h_all))
            zmax = maximum(h_all) + 0.1*abs(maximum(h_all))
            axis([xmin, xmax, zmin, zmax])
            xlabel(L"$t^*$")
            ylabel(L"$h/c$")

            subplot2grid((4, 5), (1, 0), colspan=4)
            tev, _ = DelimitedFiles.readdlm("tev-$root_ind", '\t', Float64, header=true)
            lev =   try
                DelimitedFiles.readdlm("lev-$root_ind", '\t', Float64, header=true)[1]
            catch
                zeros(0,3)
            end
            bv, _ = DelimitedFiles.readdlm("boundv-$root_ind", '\t', Float64, header=true)
            viewVort2D(tev, lev, bv)
            axis([xminv-1, xmaxv+1, zminv-1, zmaxv+1])
            ylabel("Root")

            subplot2grid((4, 5), (2, 0), colspan=4)
            tev, _ = DelimitedFiles.readdlm("tev-$mid_ind", '\t', Float64, header=true)
            lev =   try
                DelimitedFiles.readdlm("lev-$mid_ind", '\t', Float64, header=true)[1]
            catch
                zeros(0,3)
            end
            bv, _ = DelimitedFiles.readdlm("boundv-$mid_ind", '\t', Float64, header=true)
            viewVort2D(tev, lev, bv)
            axis([xminv-1, xmaxv+1, zminv-1, zmaxv+1])
            ylabel("Mid")

            subplot2grid((4, 5), (3, 0), colspan=4)
            tev, _ = DelimitedFiles.readdlm("tev-$tip_ind", '\t', Float64, header=true)
            lev =   try
                DelimitedFiles.readdlm("lev-$tip_ind", '\t', Float64, header=true)[1]
            catch
                zeros(0,3)
            end
            bv, _ = DelimitedFiles.readdlm("boundv-$tip_ind", '\t', Float64, header=true)
            viewVort2D(tev, lev, bv)
            axis([xminv-1, xmaxv+1, zminv-1, zmaxv+1])
            ylabel("Tip")
            
            #Force plots
            subplot2grid((4, 5), (0, 2))
            plot(t, cl)
            range = Int(round(0.05*len)):Int(round(0.95*len))
            xmin = t[1]
            xmax = t[end]
            zmin = minimum(cl[range]) - 0.1*abs(minimum(cl[range]))
            zmax = maximum(cl[range]) + 0.1*abs(maximum(cl[range]))
            axis([xmin, xmax, zmin, zmax])
            xlabel(L"$t^*$")
            ylabel(L"$C_L$")

            subplot2grid((4, 5), (0, 3))
            plot(t, cdr)
            range = Int(round(0.05*len)):Int(round(0.95*len))
            xmin = t[1]
            xmax = t[end]
            zmin = minimum(cdr[range]) - 0.1*abs(minimum(cdr[range]))
            zmax = maximum(cdr[range]) + 0.1*abs(maximum(cdr[range]))
            axis([xmin, xmax, zmin, zmax])
            xlabel(L"$t^*$")
            ylabel(L"$C_D$")

            subplot2grid((4, 5), (0, 4))
            plot(t, cm)
            range = Int(round(0.05*len)):Int(round(0.95*len))
            xmin = t[1]
            xmax = t[end]
            zmin = minimum(cm[range]) - 0.1*abs(minimum(cm[range]))
            zmax = maximum(cm[range]) + 0.1*abs(maximum(cm[range]))
            axis([xmin, xmax, zmin, zmax])
            xlabel(L"$t^*$")
            ylabel(L"$C_M$")


            #Spanwise variations
            spanvar, _ = DelimitedFiles.readdlm("spanwise-var", '\t', Float64, header=true)
            
            subplot2grid((4, 5), (1, 4)) 
            plot(spanvar[:,1], spanvar[:,2])
            plot(spanvar[:,1], spanvar[:,3])
            xmin = minimum(spanvar[:,1])
            xmax = maximum(spanvar[:,1])
            zmin = vormin
            zmax = vormax
            axis([xmin, xmax, zmin, zmax])
            xlabel(L"$y_{LE}$")
            ylabel(L"tev, lev")
            
            subplot2grid((4, 5), (2, 4)) 
            matind = argmin(abs.(t .- dirresults[i]))
            cl_load = [0.; mat[matind, ([1:nspan;].-1).*8 .+ 9]]
            plot(spanvar[:,1], spanvar[:,4]/spanvar[end,4])
            plot(spanvar[:,1], cl_load)
xmin = minimum(spanvar[:,1])
xmax = maximum(spanvar[:,1])
zmin = gammin
zmax = gammax
axis([xmin, xmax, zmin, zmax])
xlabel(L"$y_{LE}$")
ylabel(L"$\Gamma / \Gamma_r$, $C_l c/C_L c_{ref}$")

subplot2grid((4, 5), (3, 4)) 
plot(spanvar[:,1], spanvar[:,5])
plot(spanvar[:,1], spanvar[:,6])
plot(spanvar[:,1], spanvar[:,7])
xmin = minimum(spanvar[:,1])
xmax = maximum(spanvar[:,1])
zmin = a0min
zmax = a0max
axis([xmin, xmax, zmin, zmax])
xlabel(L"$y_{LE}$")
ylabel(L"$A0_{2D}$, $A0_{3D}$, $A0_{tot}$")

tight_layout()

savefig("../infoPlots/$(dirresults[i]).png")
close()
cd("..")
end
end
end
