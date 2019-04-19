using JFVM

Nx = 64
Lx = 3.1413
m = createMesh1D(Nx, Lx)



# BC = createBC(m)
# BC.left.a[:]=BC.right.a[:]=0.0
# BC.left.b[:]=BC.right.b[:]=1.0
# BC.left.c[:]=1.0
# BC.right.c[:]=0.0
# c_init = 0.0 # initial value of the variable
# c_old = createCellVariable(m, 0.0, BC)
# D_val = 1.0 # value of the diffusion coefficient
# D_cell = createCellVariable(m, D_val) # assigned to cells
# # Harmonic average
# D_face = harmonicMean(D_cell)
# N_steps = 20 # number of time steps
# dt= sqrt(Lx^2/D_val)/N_steps # time step
# M_diff = diffusionTerm(D_face) # matrix of coefficient for diffusion term
# (M_bc, RHS_bc)=boundaryConditionTerm(BC) # matrix of coefficient and RHS for the BC
# for i =1:5
#     (M_t, RHS_t)=transientTerm(c_old, dt, 1.0)
#     M=M_t-M_diff+M_bc # add all the [sparse] matrices of coefficient
#     RHS=RHS_bc+RHS_t # add all the RHS's together
#     c_old = solveLinearPDE(m, M, RHS) # solve the PDE
# end
