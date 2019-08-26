"""
    compareforces(cfdfile, ldvmfile; cfdskip, ldvmskip, periodic, k, u, c, cfdcol, ldvmcol, figfile)
    
Compares forces between UNSflow simulations and CFD/other data 

cfdfile and ldvm file are paths to the datafiles

skip denotes the number of header lines to be skipped

periodic = 1 for sinusoidal kinematics (k only used in this case)

col denotes the column containing data to be compares

figfile is the output filename of figure

"""
function compareForces(cfdfile, ldvmfile; cfdskip=10, ldvmskip=1, periodic=0, k=1, u=0.1312, c=0.0762, cfdcol=4, ldvmcol=6, figfile="fig.eps")
    
    cfd = readdlm(cfdfile, skipstart=cfdskip)
    ldvm = readdlm(ldvmfile, skipstart=ldvmskip)
    
    if periodic == 1
        kubyc = k*u/c
        cfd = getEndCycle(cfd, kubyc, "cfdEndCycle")
        ldvm = getEndCycle(ldvm, k, "ldvmEndCycle")
    else
        cfd[:,1] *= u/c
    end
    
    cfd_t = cfd[:,1]
    cfd_q = cfd[:,cfdcol]
    ldvm_t = ldvm[:,1]
    ldvm_q = ldvm[:,ldvmcol]
    
    figure(1)
    plot(cfd_t, cfd_q, "k-")
    plot(ldvm_t, ldvm_q, "b--")
    
    savefig(figfile)

end

function compare_cp_edgevel(cfdfile, ldvmfile; cfdskip=1, ldvmskip=1, c=0.0762, u=0.1312, tstar, kinem="undefined", amp=0, cfdx=5, cfdy=6, cfdq1=2, cfdq2=3, cfdp=1, t=zeros(1), alpha=zeros(1), qfile="q.eps", cpfile="cp.eps")

    cfd = readdlm(cfdfile, ',', skipstart=cfdskip)
    ldvm = readdlm(ldvmfile, skipstart=ldvmskip)
    
    if kinem == "plunge"
    
        x = cfd[:,cfdx]
        y = cfd[:,cfdy] .- amp*sin(2*pi*tstar)*c

        cfd_us = cfd[y.>=0,:]
        cfd_ls = cfd[y.<0,:]

        cfd_xu = cfd_us[:,cfdx]/c
        p = sortperm(cfd_xu)
        cfd_xu = cfd_xu[p]
        cfd_qu = sqrt.(cfd_us[p,cfdq1].^2 .+ cfd_us[p,cfdq2].^2)/u
        cfd_cpu = cfd_us[p,cfdp]/(0.5*u^2)
        
        cfd_xl = cfd_ls[:,cfdx]/c
        p = sortperm(cfd_xl)
        cfd_xl = cfd_xl[p]
        cfd_ql = sqrt.(cfd_ls[p,cfdq1].^2 .+ cfd_ls[p,cfdq2].^2)/u
        cfd_cpl = cfd_ls[p,cfdp]/(0.5*u^2)
        
    elseif kinem == "pitch"

        alpha_n = amp*sin(2*pi*tstar)
        var = ([cos(alpha_n) -sin(alpha_n); sin(alpha_n) cos(alpha_n)]*[cfd[:,cfdx]';cfd[:,cfdy]'])'
        x = var[:,1]
        y = var[:,2]

        cfd_us = cfd[y.>=0,:]
        cfd_ls = cfd[y.<0,:]
        
        cfd_xu = cfd_us[:,cfdx]
        cfd_yu = cfd_us[:,cfdy]
        var = ([cos(alpha_n) -sin(alpha_n); sin(alpha_n) cos(alpha_n)]*[cfd_xu';cfd_yu'])'
        cfd_xu = var[:,1]
        cfd_yu = var[:,2]

        p = sortperm(cfd_xu)
        cfd_xu = cfd_xu[p]/c
        cfd_qu = sqrt.(cfd_us[p,cfdq1].^2 .+ cfd_us[p,cfdq2].^2)/u
        cfd_cpu = cfd_us[p,cfdp]/(0.5*u^2)

        cfd_xl = cfd_ls[:,cfdx]
        cfd_yl = cfd_ls[:,cfdy]
        var = ([cos(alpha_n) -sin(alpha_n); sin(alpha_n) cos(alpha_n)]*[cfd_xl';cfd_yl'])'
        cfd_xl = var[:,1]
        cfd_yl = var[:,2]

        p = sortperm(cfd_xl)
        cfd_xl = cfd_xl[p]/c
        cfd_ql = sqrt.(cfd_ls[p,cfdq1].^2 .+ cfd_ls[p,cfdq2].^2)/u
        cfd_cpl = cfd_ls[p,cfdp]/(0.5*u^2)
        
    elseif kinem == "eld"

        alphaspl = Spline1D(t, alpha)
        alpha_n = evaluate(alphaspl, tstar)

        var = ([cos(alpha_n) -sin(alpha_n); sin(alpha_n) cos(alpha_n)]*[cfd[:,cfdx]';cfd[:,cfdy]'])'
        x = var[:,1]
        y = var[:,2]

        cfd_us = cfd[y.>=0,:]
        cfd_ls = cfd[y.<0,:]
        
        cfd_xu = cfd_us[:,cfdx]
        cfd_yu = cfd_us[:,cfdy]
        var = ([cos(alpha_n) -sin(alpha_n); sin(alpha_n) cos(alpha_n)]*[cfd_xu';cfd_yu'])'
        cfd_xu = var[:,1]
        cfd_yu = var[:,2]

        p = sortperm(cfd_xu)
        cfd_xu = cfd_xu[p]/c
        cfd_qu = sqrt.(cfd_us[p,cfdq1].^2 .+ cfd_us[p,cfdq2].^2)/u
        cfd_cpu = cfd_us[p,cfdp]/(0.5*u^2)

        cfd_xl = cfd_ls[:,cfdx]
        cfd_yl = cfd_ls[:,cfdy]
        var = ([cos(alpha_n) -sin(alpha_n); sin(alpha_n) cos(alpha_n)]*[cfd_xl';cfd_yl'])'
        cfd_xl = var[:,1]
        cfd_yl = var[:,2]

        p = sortperm(cfd_xl)
        cfd_xl = cfd_xl[p]/c
        cfd_ql = sqrt.(cfd_ls[p,cfdq1].^2 .+ cfd_ls[p,cfdq2].^2)/u
        cfd_cpl = cfd_ls[p,cfdp]/(0.5*u^2)
    else
        error("undefined kinematic type")
    end

    ldvm_x = ldvm[:,1]
    ldvm_cpu = ldvm[:,2]  
    ldvm_cpl = ldvm[:,3]
    ldvm_qu = ldvm[:,4]
    ldvm_ql = ldvm[:,5]    

    figure(1)
    clf()
    plot(cfd_xu, cfd_qu, "k*")
    plot(cfd_xl, cfd_ql, "b*")
    plot(ldvm_x, ldvm_qu, "k-")
    plot(ldvm_x, ldvm_ql, "b-")

    savefig(qfile)

    figure(2)
    clf()
    plot(cfd_xu, cfd_cpu, "k*")
    plot(cfd_xl, cfd_cpl, "b*")
    plot(ldvm_x, ldvm_cpu, "k-")
    plot(ldvm_x, ldvm_cpl, "b-")

    savefig(cpfile)    
    
end 
