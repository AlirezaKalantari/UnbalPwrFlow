using Pkg

#Pkg.add("JuMP")
using JuMP

#Pkg.add("PowerModels")
using PowerModels

#Pkg.add("Ipopt")
using Ipopt

#Pkg.add("Cbc")
using Cbc




function sphere_vol(r)

    return 4/3*pi*r^3
end


quadratic(a, sqr_t, b) = (-b + sqr_t) / 2a


function quadratic1(a::Float64, b::Float64, c::Float64)

    sqr_t = sqrt(b^2-4a*c)
    r1 = quadratic(a, sqr_t, b)
    r2 = quadratic(a, -sqr_t, b)

    r1, r2
end

vol = sphere_vol(3)


quad1, quad2 = quadratic1(1.0, -3.0, -28.0)
println("result 1: ", quad1)

println("result 2: ", quad2)
println("result 3: ", vol)
