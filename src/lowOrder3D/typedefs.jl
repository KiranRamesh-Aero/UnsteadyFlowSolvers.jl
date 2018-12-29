# mutable struct KinemDef3D
#     alpha :: MotionDef
#     h :: MotionDef
#     u :: MotionDef
#     vartype :: String
#     vary ::Int8
#     add :: Array{Float64,2}

#     function KinemDef3D(alpha :: MotionDef, h::MotionDef, u::MotionDef, vartype = "Constant", vary = 0, add = zeros(1,1))
#         new(alpha, h, u, vartype, vary, add)
#     end
# end

struct ThreeDFieldSimple
    f2d :: Vector{TwoDFlowField}
    function ThreeDFieldSimple()
        f2d = TwoDFlowField[]
    new(f2d)
    end
end

struct ThreeDSurfSimple
    cref :: Float64
    AR :: Float64
    uref :: Float64
    pvt :: Float64
    lespcrit :: Vector{Float64}
    coord_file :: String
    ndiv :: Int8
    nspan :: Int8
    naterm :: Int8
    kindef :: KinemDef
    psi :: Vector{Float64}
    yle :: Vector{Float64}
    s2d :: Vector{TwoDSurf}
    a03d :: Vector{Float64}
    bc :: Vector{Float64}
    nshed :: Vector{Float64}
    bcoeff :: Vector{Float64}
    aterm3d :: Array{Float64}
    levstr :: Vector{Float64}

    function ThreeDSurfSimple(AR, kindef, coord_file, pvt, lespcrit = [10.;]; nspan = 12, cref = 1., uref=1., ndiv=70, naterm=35)

        bref = AR*cref

        psi = zeros(nspan)
        yle = zeros(nspan)
        
        s2d = TwoDSurf[]
        
	for i = 1:nspan
            psi[i] = real(i)*(pi/2)/nspan
            yle[i] = -bref*cos(psi[i])/2.
        end
        
        for i = 1:nspan
            kinem2d = KinemDef(kindef.alpha, kindef.h, kindef.u)
            lespc = lespcrit[1]
            push!(s2d, TwoDSurf(coord_file, pvt,  kinem2d, [lespc;], c=cref, uref=uref, ndiv=ndiv, naterm=naterm))
        end
    
        a03d = zeros(nspan)
        aterm3d = zeros(naterm, nspan)
        
        bc = zeros(nspan)
        nshed = [0.;]
        bcoeff = zeros(nspan)
        levstr = zeros(nspan)
        
        new(cref, AR, uref, pvt, lespcrit, coord_file,  ndiv, nspan, naterm, kindef, psi, yle, s2d, a03d, bc, nshed, bcoeff, aterm3d, levstr)
    end
end

struct KelvinConditionQSLLT
    surf :: ThreeDSurfSimple
    field :: ThreeDFieldSimple
end

function (kelv::KelvinConditionQSLLT)(tev_iter::Array{Float64})

    val = zeros(kelv.surf.nspan)
    
    for i = 1:kelv.surf.nspan
        uprev, wprev = ind_vel([kelv.field.f2d[i].tev[end]], kelv.surf.s2d[i].bnd_x, kelv.surf.s2d[i].bnd_z)
        
        kelv.field.f2d[i].tev[end].s = tev_iter[i]

        unow, wnow = ind_vel([kelv.field.f2d[i].tev[end]], kelv.surf.s2d[i].bnd_x, kelv.surf.s2d[i].bnd_z)
        
        kelv.surf.s2d[i].uind[:] = kelv.surf.s2d[i].uind[:] .- uprev .+ unow
        kelv.surf.s2d[i].wind[:] = kelv.surf.s2d[i].wind[:] .- wprev .+ wnow
        
        #Calculate downwash
        update_downwash(kelv.surf.s2d[i], [kelv.field.f2d[i].u[1], kelv.field.f2d[i].w[1]])
        
        #Calculate first two fourier coefficients
        update_a0anda1(kelv.surf.s2d[i])
        
        kelv.surf.bc[i] = kelv.surf.s2d[i].a0[1] + 0.5*kelv.surf.s2d[i].aterm[1]
    end

    
    calc_a0a13d(kelv.surf)
    
    for i = 1:kelv.surf.nspan
        val[i] = kelv.surf.s2d[i].uref*kelv.surf.s2d[i].c*pi*(kelv.surf.bc[i] + kelv.surf.a03d[i] + 0.5*kelv.surf.aterm3d[1,i] - kelv.surf.s2d[i].a0prev[1] - 0.5*kelv.surf.s2d[i].aprev[1]) + kelv.field.f2d[i].tev[end].s
    end

    return val
end

struct KelvinKuttaQSLLT
    surf :: ThreeDSurfSimple
    field :: ThreeDFieldSimple
    shed_ind :: Vector{Int}
end

function (kelv::KelvinKuttaQSLLT)(v_iter::Array{Float64})

    nshed = length(kelv.shed_ind)
    val = zeros(kelv.surf.nspan+nshed)
    
    levcount = 0
    for i = 1:kelv.surf.nspan
        if i in kelv.shed_ind
            uprev, wprev = ind_vel([kelv.field.f2d[i].tev[end]; kelv.field.f2d[i].lev[end]], kelv.surf.s2d[i].bnd_x, kelv.surf.s2d[i].bnd_z)

            kelv.surf.levstr[i] -= kelv.field.f2d[i].lev[end].s 
            
            kelv.field.f2d[i].tev[end].s = v_iter[i]
            levcount += 1
            kelv.field.f2d[i].lev[end].s = v_iter[kelv.surf.nspan+levcount]

            kelv.surf.levstr[i] += kelv.field.f2d[i].lev[end].s 
            
            unow, wnow = ind_vel([kelv.field.f2d[i].tev[end]; kelv.field.f2d[i].lev[end]], kelv.surf.s2d[i].bnd_x, kelv.surf.s2d[i].bnd_z)
            
        else
            uprev, wprev = ind_vel([kelv.field.f2d[i].tev[end]], kelv.surf.s2d[i].bnd_x, kelv.surf.s2d[i].bnd_z)

            kelv.field.f2d[i].tev[end].s = v_iter[i]
            
            unow, wnow = ind_vel([kelv.field.f2d[i].tev[end]], kelv.surf.s2d[i].bnd_x, kelv.surf.s2d[i].bnd_z)
        end

        kelv.surf.s2d[i].uind[:] = kelv.surf.s2d[i].uind[:] .- uprev .+ unow
        kelv.surf.s2d[i].wind[:] = kelv.surf.s2d[i].wind[:] .- wprev .+ wnow
        
        #Calculate downwash
        update_downwash(kelv.surf.s2d[i], [kelv.field.f2d[i].u[1], kelv.field.f2d[i].w[1]])
        
        #Calculate first two fourier coefficients
        update_a0anda1(kelv.surf.s2d[i])
        
        kelv.surf.bc[i] = kelv.surf.s2d[i].a0[1] + 0.5*kelv.surf.s2d[i].aterm[1]
    end
    
    #calc_a0a13d_wlev(kelv.surf, kelv.shed_ind, v_iter[kelv.surf.nspan+1:end])
    calc_a0a13d_wlev(kelv.surf)

    levcount = 0
    for i = 1:kelv.surf.nspan
        if i in kelv.shed_ind
            val[i] = kelv.surf.s2d[i].uref*kelv.surf.s2d[i].c*pi*(kelv.surf.bc[i] + kelv.surf.a03d[i] + 0.5*kelv.surf.aterm3d[1,i] - kelv.surf.s2d[i].a0prev[1] - 0.5*kelv.surf.s2d[i].aprev[1]) + kelv.field.f2d[i].tev[end].s + kelv.field.f2d[i].lev[end].s
            levcount += 1
            if kelv.surf.s2d[i].a0[1] > 0
                lesp_cond = kelv.surf.s2d[i].lespcrit[1]
            else
                lesp_cond = -kelv.surf.s2d[i].lespcrit[1]
            end
            val[kelv.surf.nspan+levcount] = kelv.surf.s2d[i].a0[1] + kelv.surf.a03d[i] - lesp_cond
            
        else
            val[i] = kelv.surf.s2d[i].uref*kelv.surf.s2d[i].c*pi*(kelv.surf.bc[i] + kelv.surf.a03d[i] + 0.5*kelv.surf.aterm3d[1,i] - kelv.surf.s2d[i].a0prev[1] - 0.5*kelv.surf.s2d[i].aprev[1]) + kelv.field.f2d[i].tev[end].s
        end
    end

    
    return val
end

struct KelvinConditionQSLLT_sep
    surf :: ThreeDSurfSimple
    field :: ThreeDFieldSimple
    ispan :: Int
end

function (kelv::KelvinConditionQSLLT_sep)(tev_iter::Array{Float64})

    i = kelv.ispan
    
    uprev, wprev = ind_vel([kelv.field.f2d[i].tev[end]], kelv.surf.s2d[i].bnd_x, kelv.surf.s2d[i].bnd_z)
    
    kelv.field.f2d[i].tev[end].s = tev_iter[1]

    unow, wnow = ind_vel([kelv.field.f2d[i].tev[end]], kelv.surf.s2d[i].bnd_x, kelv.surf.s2d[i].bnd_z)
    
    kelv.surf.s2d[i].uind[:] = kelv.surf.s2d[i].uind[:] .- uprev .+ unow
    kelv.surf.s2d[i].wind[:] = kelv.surf.s2d[i].wind[:] .- wprev .+ wnow
    
    #Calculate downwash
    update_downwash(kelv.surf.s2d[i], [kelv.field.f2d[i].u[1], kelv.field.f2d[i].w[1]])
    
    #Calculate first two fourier coefficients
    update_a0anda1(kelv.surf.s2d[i])
    
    kelv.surf.bc[i] = kelv.surf.s2d[i].a0[1] + 0.5*kelv.surf.s2d[i].aterm[1]
    
    val = kelv.surf.s2d[i].uref*kelv.surf.s2d[i].c*pi*(kelv.surf.bc[i] + kelv.surf.a03d[i] + 0.5*kelv.surf.aterm3d[1,i] - kelv.surf.s2d[i].a0prev[1] - 0.5*kelv.surf.s2d[i].aprev[1]) + kelv.field.f2d[i].tev[end].s

    return val
end

struct KelvinKuttaQSLLT_sep
    surf :: ThreeDSurfSimple
    field :: ThreeDFieldSimple
    ispan :: Int
end

function (kelv::KelvinKuttaQSLLT_sep)(v_iter::Array{Float64})

    i = kelv.ispan
    val = zeros(2)
    
    uprev, wprev = ind_vel([kelv.field.f2d[i].tev[end]; kelv.field.f2d[i].lev[end]], kelv.surf.s2d[i].bnd_x, kelv.surf.s2d[i].bnd_z)
    
    kelv.field.f2d[i].tev[end].s = v_iter[1]
    kelv.field.f2d[i].lev[end].s = v_iter[2]
    
    unow, wnow = ind_vel([kelv.field.f2d[i].tev[end]; kelv.field.f2d[i].lev[end]], kelv.surf.s2d[i].bnd_x, kelv.surf.s2d[i].bnd_z)

    kelv.surf.s2d[i].uind[:] = kelv.surf.s2d[i].uind[:] .- uprev .+ unow
    kelv.surf.s2d[i].wind[:] = kelv.surf.s2d[i].wind[:] .- wprev .+ wnow
    
    #Calculate downwash
    update_downwash(kelv.surf.s2d[i], [kelv.field.f2d[i].u[1], kelv.field.f2d[i].w[1]])
    
    #Calculate first two fourier coefficients
    update_a0anda1(kelv.surf.s2d[i])
    
    kelv.surf.bc[i] = kelv.surf.s2d[i].a0[1] + 0.5*kelv.surf.s2d[i].aterm[1]
    
    val[1] = kelv.surf.s2d[i].uref*kelv.surf.s2d[i].c*pi*(kelv.surf.bc[i] + kelv.surf.a03d[i] + 0.5*kelv.surf.aterm3d[1,i] - kelv.surf.s2d[i].a0prev[1] - 0.5*kelv.surf.s2d[i].aprev[1]) + kelv.field.f2d[i].tev[end].s + kelv.field.f2d[i].lev[end].s

    if kelv.surf.s2d[i].a0[1] > 0
        lesp_cond = kelv.surf.s2d[i].lespcrit[1]
    else
        lesp_cond = -kelv.surf.s2d[i].lespcrit[1]
    end

    val[2] = kelv.surf.s2d[i].a0[1] + kelv.surf.a03d[i] - lesp_cond
            
    return val
end

mutable struct KinemCantilever
    eta :: Array{Float64}
    etad :: Array{Float64}
    etadd :: Array{Float64}
    eta_pr :: Array{Float64}
    etad_pr :: Array{Float64}
    etadd_pr :: Array{Float64}
    eta_pr2 :: Array{Float64}
    etad_pr2 :: Array{Float64}
    etadd_pr2 :: Array{Float64}
    eta_pr3 :: Array{Float64}
    etad_pr3 :: Array{Float64}
    etadd_pr3 :: Array{Float64}

    function KinemCantilever(eta, etad)
        z = zeros(length(eta))
        new(eta, etad, z, eta, etad, z, eta, etad, z, eta, etad, z)
    end
end
