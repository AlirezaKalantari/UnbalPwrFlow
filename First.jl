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
Y=zeros(4,4)
J=zeros(4,4)
E_k=[1+0im ,1+0im ,1+0im ,1+0im]
I_cal=[1+0im ,1+0im ,1+0im ,1+0im]

Pg=zeros(size(G))
Pl=zeros(size(L))
P_sp=zeros(size(G))
P_cal=zeros(size(G))
delP=zeros(size(G))

Qg=zeros(size(G))
Ql=zeros(size(L))
Q_sp=zeros(size(G))
Q_cal=zeros(size(G))
delQ=zeros(size(G))


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

P_sp=broadcast(-,Pg,Pl)'
Q_sp=broadcast(-,Qg,Ql)'

Ybus=Y

V_rk=real(E_k)
V_mk=imag(E_k)

G_y=real(Ybus)
B_y=imag(Ybus)

for g in G
    P_cal[g]=V_rk[g]*real(I_cal)[g]'
    Q_cal[g]=V_mk[g]*imag(I_cal)[g]'

<<<<<<< HEAD
    delP[g]=P_sp[g]-P_cal[g]
    delQ[g]=Q_sp[g]-Q_cal[g]

=======
>>>>>>> b18e318ce10efb970ac3bafbf6b76a453c2733c9
end







<<<<<<< HEAD

=======
>>>>>>> b18e318ce10efb970ac3bafbf6b76a453c2733c9




real(I_cal)'




I=zeros(size(G))



    #println("This is dispatch of unit $(g): $(JuMP.value(P[g]))")
#end


println("salam saeed jan, How can Import sum data such as Pg, Pl, Ybus in julia?")
println("Salaam, Check out PowerModels.jl")
