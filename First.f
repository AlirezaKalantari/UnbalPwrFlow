using Pkg

#Pkg.add("JuMP")
using JuMP

#Pkg.add("LinearAlgebra")
using LinearAlgebra

#Pkg.add("PowerModels")
using PowerModels

#Pkg.add("Ipopt")
using Ipopt

#Pkg.add("Cbc")
using Cbc

data = PowerModels.parse_file("C:/Users/kalan/Desktop/matpower7.0/data/case5.m")
data["bus"]
data["branch"]
data["gen"]
data["load"]

display(data)
PowerModels.print_summary(data)
ct1 =PowerModels.component_table(data, "bus", "bus_type")'
ct2=PowerModels.component_table(data, "bus", ["vmin", "vmax"])
ct3=PowerModels.component_table(data, "gen", ["pmin", "pmax", "qmin", "qmax"])
ct4=PowerModels.component_table(data, "branch", ["f_bus","t_bus","br_r","br_x"])
ct5=PowerModels.component_table(data, "load", ["load_bus","pd","qd"])
ct6=PowerModels.component_table(data, "gen", ["gen_bus","pg","qg"])


a=size(PowerModels.component_table(data, "bus", "bus_type"),1)
zb_re=zeros(a,a)
yb_re=zeros(a,a)

zb_im=zeros(a,a)
zb_im=complex(zb_im);

yb_im=zeros(a,a)
yb_im=complex(zb_im);

z_RE=ones(a,a)


q=collect(1:(a+1))

for g in q
        yb_re[ct4[g,2],ct4[g,3]]=-real(1/(ct4[g,4]+1im*ct4[g,5]))

end

for l=1:a
    for m=1:a
        if yb_re[l,m]!==0
            yb_re[m,l]=yb_re[l,m]
        end
    end
end

for l=1:a
    yb_re[l,l]=-(yb_re[l,1]+yb_re[l,2]+yb_re[l,3]+yb_re[l,4]+yb_re[l,5])
end

for g in q
    yb_im[ct4[g,2],ct4[g,3]]=-imag(1/(ct4[g,4]+1im*ct4[g,5]))
end

for l=1:a
    for m=1:a
        if yb_im[l,m]!==0
            yb_im[m,l]=yb_im[l,m]
        end
    end
end

for l=1:a
    yb_im[l,l]=-(yb_im[l,1]+yb_im[l,2]+yb_im[l,3]+yb_im[l,4]+yb_im[l,5])
end

YBUS=yb_re+(1im*yb_im)


iter=0

PV_bus=zeros(1,1)
Iter=collect(1:1000)
G=collect(1:(4))    #number of generator #G=collect(1:(24))
j=collect(1:(8))
L=collect(1:4)    #number of load
Y=ones(4,4)
J=zeros(4,4)
E_k=[1+0im ,1+0im ,1+0im ,1+0im] #E_k=ones(a,1); #E_k=E_k(complex)
I_cal=[1+0im ,1+0im ,1+0im ,1+0im]  #I_cal=ones(a,1);   ##I_cal=I_cal(complex)
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
