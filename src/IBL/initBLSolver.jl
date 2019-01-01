function initiliaseViscousPara(ncell::Int64, cfl::Float64, re:: Float64)

    invisicidPara = InvisicidTransport(ncell)
    fluxSplittingPara = FluxSplittingParameters(ncell)
    soln =Solutions(ncell, fluxSplittingPara)
    opCond= OperationalConditions(cfl, re, ncell)

    return invisicidPara, fluxSplittingPara, soln, opCond

end


function coupleInviscidViscous(surf::TwoDSurf, invisidPara::InvisicidTransport, curfield::TwoDFlowField, opCond::OperationalConditions)

    q_u,q_l = UNSflow.calc_edgeVel(surf, [curfield.u[1], curfield.w[1]])

    ncell= length(q_u)
    x = zeros(ncell+2)
    dx= 1.0/(convert(Float64,ncell+2))

    for i=2:ncell+1
        x[i] = (convert(Float64,i)-1.5)*dx
        invisidPara.ue_us[i] = q_u[i-1]
        invisidPara.ue_ls[i] = q_l[i-1]
    end

    x[1] =2*x[2]-x[3]
    x[ncell+2] = 2*x[ncell+1]-x[ncell]

    #x[1] =0.
    #x[ncell+2] = 1.
    #invisidPara.ue_us[1]= 2.0*invisidPara.ue_us[2]-invisidPara.ue_us[3]
    #invisidPara.ue_us[ncell+2]= 2.0*invisidPara.ue_us[ncell+1]-invisidPara.ue_us[ncell+1]
    #invisidPara.ue_us[1]= 0.09512534545251
    #invisidPara.ue_us[2]= 0.009512534545251
    #invisidPara.ue_us[ncell+2]= 0.9947793650842671

    invisidPara.ue_ls[1]= 2.0*invisidPara.ue_ls[2]-invisidPara.ue_ls[3]
    invisidPara.ue_ls[ncell+2]= 2.0*invisidPara.ue_ls[ncell+1]-invisidPara.ue_ls[ncell+1]

    opCond.dx=dx
    opCond.x[:]= x[:]

 end
