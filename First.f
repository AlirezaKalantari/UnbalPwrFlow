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
ct3=PowerModels.component_table(data, "gen", ["pg","pmin", "pmax", "qmin", "qmax"])
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
d=size(ct3,1)*ones(1,1)

for i=(1:size(ct3,1))
  b=size(ct3,1)
    if ct3[i,2]!==0.00
      d[1,1]=b-1
          end
end
real_pg=d[1,1]

iter=0


PV_bus=zeros(1,1)
Iter=collect(1:1000)
G=collect(1:real_pg)    #number of generator #G=collect(1:(24))
j=collect(1:(2*a))
L=collect(1:real_pg)    #number of load
Y=ones(a,a)
J=zeros(a,a)
E_k=complex(ones(a,1))
I_cal=complex(ones(a,1))
β_k=zeros(a,a)
γ_k=zeros(a,a)
δ_k=zeros(a,a)
delI_r=zeros(a,1)
delI_m=zeros(a,1)
Pg=zeros(a,1)
Pl=zeros(a,1)
P_sp=zeros(a,1)
P_cal=zeros(a,1)
delP=zeros(size(G))
del_I=zeros((2*a),1)
del_V=zeros((2*a),1)
newdel_V=zeros((2*a),1)
del_v=zeros((2*a),1)

Qg=zeros(a,1)
Ql=zeros(a,1)
Q_sp=zeros(a,1)
Q_cal=zeros(a,1)
delQ=zeros(size(G))
tete=zeros(a,1)

for i=1:(size(ct5,1))
    Pl[ct5[i,2]]=ct5[i,3]
    Ql[ct5[i,2]]=ct5[i,4]
end

for i=1:a
    Pg[ct6[i,2]]=ct6[i,3]
    Qg[ct6[i,2]]=ct6[i,4]
end

P_sp=broadcast(-,Pg,Pl)'
Q_sp=broadcast(-,Qg,Ql)'

for i=1:a
    if ct1[2,i]==1
        println("Bus=s")
    end
end
