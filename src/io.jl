# JSON

function writejson(path, obj)
    open(path, "w") do io
        JSON.print(io, obj, 4)
    end
end
function readjson(path)
    f = read(path, String)
    return JSON.parse(f)
end


# SHOW

function beautifyshow(s::String; space = Iterators.repeat(" ",4))
    output = ""
    index = 1
    tab = 0
    type = 0
    for si in s
        output *= si
        (si == '}') && (type -= 1)
        (si == '{') && (type += 1)
        type != 0 && continue
        if si == '('
            tab += 1
            output *= "\n"*Iterators.repeat(space,tab)
        end
        (si == ')') && (tab -= 1)
        (si == ',') && (output *= "\n"*Iterators.repeat(space,tab))
    end
    return output
end
function writeshow(path, obj)
    open(path, "w") do io
        s = sprint(show, obj)
        println(io, beautifyshow(s))
    end
end
function readshow(path)
    f = read(path, String)
    expr = Meta.parse(f)
    return eval(expr)
end


# Measurements in Dictionary

function transformdictrecursively!(d::Dict, apply)
    for (k,v) in d
        if v isa Dict
            d[k] = transformdictrecursively!(v, apply)
            continue
        end
        if v isa NamedTuple
            d[k] = NamedTuple{keys(v)}(apply.(Tuple(v)))
            continue
        end
        d[k] = apply.(v)
    end
    return d
end

function ifstringgivemeasurement(x)
    if (x isa String) && (findfirst("±", x) !== nothing)
        return eval(Meta.parse(x))
    end
    return x
end

function ifmeasurementgivestring(x)
    if x isa Measurement
        return string(x)
    end
    return x
end


# dictionary to named tuple

d2nt(x; process=identity) = process(x)
d2nt(d::Dict; process=identity) = NamedTuple{Tuple(Symbol.(keys(d)))}(d2nt.(values(d); process))
d2nt(x::Vector{Any}; process=identity) = vcat(d2nt.(x; process)...)
