abstract type delVortDef end

mutable struct delNone <:delVortDef end

mutable struct delSpalart <: delVortDef
    limit :: Int64
    dist :: Float64
    tol :: Float64
end
