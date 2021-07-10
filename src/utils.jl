(..)(x::AbstractArray,i...) = getindex.(x,i...)
(..)(x::AbstractArray,i::Symbol) = getfield.(x,i)
