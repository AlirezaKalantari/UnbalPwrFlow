using Pkg

#Pkg.add("JuMP")
using JuMP

#Pkg.add("PowerModels")
using PowerModels

#Pkg.add("Ipopt")
using Ipopt

#Pkg.add("Cbc")
using Cbc


PV_bus=zeros(1,1)
Iter=collect(1:1000)
G=collect(1:4)    #number of generator
j=collect(1:8)
L=collect(1:4)    #number of load
Y=ones(4,4)
J=zeros(4,4)
E_k=[1+0im ,1+0im ,1+0im ,1+0im]
I_cal=[1+0im ,1+0im ,1+0im ,1+0im]
α_k=zeros(4,4)
β_k=zeros(4,4)
γ_k=zeros(4,4)
δ_k=zeros(4,4)
delI_r=zeros(4,1)
delI_m=zeros(4,1)
Pg=zeros(size(G))
Pl=zeros(size(L))
P_sp=zeros(size(G))
P_cal=zeros(size(G))
delP=zeros(size(G))
del_I=zeros(8,1)
J=ones(8,8)
del_V=zeros(8,1)
newdel_V=zeros(8,1)
del_v=zeros(8,1)

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
    α_k[g,g]=(Q_sp[g]*(V_rk[g]^2-V_mk[g]^2)-2*V_mk[g]*V_rk[g]*P_sp[g])/((V_mk[g]^2+V_rk[g]^2)^2)
    β_k[g,g]=(Q_sp[g]*(V_rk[g]^2-V_mk[g]^2)+2*V_mk[g]*V_rk[g]*P_sp[g])/((V_mk[g]^2+V_rk[g]^2)^2)
    δ_k=α_k
    γ_k=-β_k

end
for g in G
    P_cal[g]=V_rk[g]*real(I_cal)[g]'+V_mk[g]*imag(I_cal)[g]'
    Q_cal[g]=V_mk[g]*real(I_cal)[g]'-V_rk[g]*imag(I_cal)[g]'

    delP[g]=P_sp[g]-P_cal[g]
    delQ[g]=Q_sp[g]-Q_cal[g]


    delI_r[g]=(delP[g]*V_rk[g]+V_mk[g]*delQ[g])/(V_mk[g]^2+V_rk[g]^2)
    delI_m[g]=(delP[g]*V_mk[g]+V_rk[g]*delQ[g])/(V_mk[g]^2+V_rk[g]^2)

    del_I[(2*g-1),1]=delI_r[g]
    del_I[(2*g),1]=delI_m[g]


end

for a in G
    for b in G
        if a==b

            J[2*a-1,b]=imag(Ybus[a,a])-α_k[a,a]
            J[2*a-1,b+1]=real(Ybus[a,a])-β_k[a,a]
            a=a+1
            J[a,b]=real(Ybus[a-1,a-1])-α_k[a-1,a-1]
            J[a,b+1]=-real(Ybus[a-1,a-1])+β_k[a-1,a-1]
            a=a-1
        elseif a!==b
            c=2*b-1
            J[a,c]=imag(Ybus[a,b])
            J[a,c+1]=real(Ybus[a,b])
            a=a+1
            J[a,c]=real(Ybus[a-1,b])
            J[a,c+1]=-imag(Ybus[a-1,b])
            a=a-1



        end

    end
end

#PV_bus
#for a in G
#for b in G
#if a!==b && Ql[a]!==0
#    c=2*b-1
#    J[a,c]=imag(Ybus[a,b])-α_k[a,a]-(real(Ybus[a,b])(V_mk[a]/V_rk[a]))
#    J[a,c+1]=real(Ybus[a,b])-β_k[a,a]
#    a=a+1
#    J[a,c]=imag(Ybus[a-1,b])+β_k[a-1,a-1]-(real(Ybus[a-1,b])(V_mk[a-1]/V_rk[a-1]))
#    J[a,c+1]=real(Ybus[a-1,b])-α_k[a-1,a-1]
#    a=a-1
#if a==b
#    J[2*a-1,b]=real(Ybus[a,a])-β_k[a,a]-((imag(Ybus[a,b])-α_k[a,a])(V_mk[a]/V_rk[a])))
#    J[2*a-1,b+1]=V_rk[a]/(V_rk[a]^2+V_mk[a]^2)
#    a=a+1
#    J[a,b]=J[2*a-1,b]=real(Ybus[a-1,a-1])-α_k[a-1,a-1]+((imag(Ybus[a-1,b])+β_k[a-1,a-1])(V_mk[a-1]/V_rk[a-1])))
#    J[a,b+1]=-V_mk[a]/(V_rk[a]^2+V_mk[a]^2)
#    a=a-1
#    delI_m[a]=V_mk[a]/(V_rk[a]^2+V_mk[a]^2)
#    delI_r[a]=V_rk[a]/(V_rk[a]^2+V_mk[a]^2)
#    del_I[(2*a-1),1]=delI_r[g]
#    del_I[(2*a),1]=delI_m[g]


#end
#end


del_v=J\del_I

for g in G
    del_V[2*g-1,1]=V_rk[g]
    del_V[2*g,1]=V_mk[g]

end
newdel_V=del_v+del_V
newdel_V=del_V

println("salam saeed jan, I tell you what has happened to me, I promise you to work harder?")
