using Pkg

#Pkg.add("JuMP")
using JuMP

#Pkg.add("PowerModels")
using PowerModels

#Pkg.add("Ipopt")
using Ipopt

#Pkg.add("Cbc")
using Cbc




G=collect(1:4)    #number of generator
L=collect(1:4)    #number of load

Pg=zeros(size(G))
Pl=zeros(size(L))
p_sp=zeros(size(G))

Qg=zeros(size(G))
Ql=zeros(size(L))
Q_sp=zeros(size(G))



for l in L
    Pl[l]=100*l
end

for l in L
    Ql[l]=50*l

end

for l in L
    Pg[l]=200*l
end

for l in L
    Qg[l]=120*l

end

P_sp=broadcast(-,broadcast(+,Pg,Qg),broadcast(+,Pl,Ql))'





I=zeros(size(G))



    println("This is dispatch of unit $(g): $(JuMP.value(p[g]))")
end


println("salam saeed jan, How can Import sum data such as Pg, Pl, Ybus in julia?")
