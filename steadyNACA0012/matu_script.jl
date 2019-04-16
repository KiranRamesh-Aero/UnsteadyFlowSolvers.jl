xf = UnsteadyFlowSolvers.createUniformFVM(200)

xf =xf/pi

qf = UnsteadyFlowSolvers.mappingAerofoilToFVGrid(qu, surf, xf)

dqdx = UnsteadyFlowSolvers.spatialDeriv(qf)

matuCyl(1.0, qf, dqdx)
