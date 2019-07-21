using Pkg

#Pkg.add("JuMP")
using JuMP

#Pkg.add("PowerModels")
using PowerModels

#Pkg.add("Ipopt")
using Ipopt

#Pkg.add("Cbc")
using Cbc


G=collect(1:2)
D=collect(1:2)

Demand=zeros(size(D))
for d in D
    Demand[d]=100*d
end

Cost=zeros(size(G))
Pmax=zeros(size((G)))
for g in G
    Cost[g]=sqrt(g)
    Pmax[g]=300/g
end

myupf=Model(with_optimizer(Cbc.Optimizer))

@variable(myupf,p[g in G]>=0)

@objective(myupf, Min, sum(Cost[g]*p[g] for g in G))


@constraint(myupf,c1[g in G],p[g]<=Pmax[g])

@constraint(myupf,c2,sum(p[g] for g in G)==sum(Demand[d] for d in D))

JuMP.optimize!(myupf)

for g in G
    println("This is dispatch of unit $(g): $(JuMP.value(p[g]))")
end
