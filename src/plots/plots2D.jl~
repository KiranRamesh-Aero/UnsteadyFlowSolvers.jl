function viewVort2D(tev::Array{Float64}, lev::Array{Float64}, bv::Array{Float64})
    scatter(tev[:,2], tev[:,3], s=5, c=tev[:,1], edgecolors="none")
    sc = scatter(lev[:,2], lev[:,3], s=5, c=lev[:,1], edgecolors="none")
    sc2 = scatter(bv[:,2], bv[:,3], s=5, c=bv[:,1], edgecolors="none") 
    plot(bv[:,2], bv[:,3], color = "black", linewidth=1.0)
end

function viewVortConnect2D(tev::Array{Float64}, lev::Array{Float64}, bv::Array{Float64})
    scatter(tev[:,2], tev[:,3], s=5, c=tev[:,1], edgecolors="none")
    sc = scatter(lev[:,2], lev[:,3], s=5, c=lev[:,1], edgecolors="none")
    sc2 = scatter(bv[:,2], bv[:,3], s=5, c=bv[:,1], edgecolors="none") 

    plot(tev[:,2], tev[:,3], color = "red", linewidth=1.0)
    plot(lev[:,2], lev[:,3], color = "blue", linewidth=1.0)
    plot(bv[:,2], bv[:,3], color = "black", linewidth=1.0)
end

function makeVortPlots2D()
    dirvec = readdir()
    dirresults = map(x->(v = tryparse(Float64,x); isnull(v) ? 0.0 : get(v)),dirvec)
    #Determine axis limits
    dirmax = maximum(dirresults)
    cd("$(dirmax)")
    tev = readdlm("tev")
    lev = try
        readdlm("lev")
            catch
        Array{Float64}(0,3)
    end
    bv = readdlm("boundv")
    
    xmin = minimim([tev[:,2];lev[:,2];bv[:,2];])
    zmin = minimum([tev[:,3];lev[:,3];bv[:,3];])
    xmax = maximum([tev[:,2];lev[:,2];])
    zmax = maximum([bv[:,2];tev[:,3];lev[:,3];bv[:,3];])
    
    if "vortPlots" in dirvec
        rm("vortPlots", recursive=true)
    end
    mkdir("vortPlots")
    for i=1:length(dirresults)
        if dirresults[i] != 0
            dirstr="$(dirresults[i])"
            cd(dirstr)
            tev = readdlm("tev")
            lev = try
                readdlm("lev")
            catch
                Array{Float64}(0,3)
            end
            bv = readdlm("boundv")

            viewVort2D(tev, lev, bv)
            axis([xmin-1, xmax+1, zmin-1, zmax+1])
            savefig("../vortPlots/$(dirresults[i]).png")
            close()
            cd("..")
        end
    end
end

