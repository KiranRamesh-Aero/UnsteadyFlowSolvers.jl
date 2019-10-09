using Roots
using SymPy

@vars E u

F = 4.8274*E.^4 - 5.9816*E.^3 + 4.0274*E.^2 + 0.23247.*E .+ 0.15174

dfde = 4*4.8274*E.^3 - 3*5.9816*E.^2 + 2*4.0274*E .+ 0.2324


det = -dfde - (F - (E + 1)*dfde)

solve(det)


