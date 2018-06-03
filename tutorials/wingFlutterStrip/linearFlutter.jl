function dyn_linear(dX, X, q, nm, ns, t)

    eta = X(1:nm,1)
    etad = X(nm+1:2*nm,1)
    lbda = X(2*nm+1:2*nm+2*ns,1)
    if t < 0.1
        q2 = 0.5*pi/180
    else
        q2 = q
    end
    etadd = (mmod.m_eta_s_etadd - mmod.m_eta_a_etadd)\((mmod.m_eta_a_etad -
    mmod.m_eta_s_etad)*etad + (mmod.m_eta_a_eta - mmod.m_eta_s_eta)*eta +
    mmod.m_eta_a_lbda*lbda + mmod.m_eta_q*q2)

    lbdad = mmod.m_lbda_etadd*etadd +
    mmod.m_lbda_etad*etad +
    mmod.m_lbda_eta*eta +
    mmod.m_lbda_lbda*lbda +
    mmod.m_lbda_q*q2

    qvec = mmod.m_eta_a_etad*etad + mmod.m_eta_a_eta*eta +
    mmod.m_eta_a_etadd*etadd + mmod.m_eta_a_lbda*lbda +
    mmod.m_eta_q*q2

    dX[1:nm,1] = etad
    dX[nm+1:2*nm,1] = etadd
    dX[2*nm+1:2*nm+2*ns,1] = lbdad

    return dXdt, qvec
end


mat = matread("flatplate_350_40_8_ballast_1_offset_0_nelem_70_aeroelast_ns_10_nU_40.mat")

ns = mat["ns"]
nm = mat["nm"]
mmod = mat["mmod"]

q = 5*pi/180
X0 = zeros(2*nm+2*ns,1)
tspan = (0.0, 3.0)
prob = ODEProblem(dyn_linear, X0, tspan, q, nm, ns, mm)
sol = solve(prob,Tsit5())

eta = X[]:,1:nm]
etad = X[]:,nm+1:2*nm]
lbda = X[]:,2*nm+1:2*nm+2*ns]
