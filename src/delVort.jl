abstract type delVortDef end

type delNone <:delVortDef end

type delSpalart <: delVortDef
    limit :: Int64
    dist :: Float64
    tol :: Float64
end

