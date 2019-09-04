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
ct1 =PowerModels.component_table(data, "bus", "va")'
ct2=PowerModels.component_table(data, "bus", ["vmin", "vmax"])
ct3=PowerModels.component_table(data, "gen", ["pmin", "pmax", "qmin", "qmax"])
ct4=PowerModels.component_table(data, "branch", ["f_bus","t_bus","br_r","br_x"])
ct5=PowerModels.component_table(data, "load", ["load_bus","pd","qd"])
ct6=PowerModels.component_table(data, "gen", ["gen_bus","pg","qg"])


a=size(PowerModels.component_table(data, "bus", "bus_type"),1)
zb_re=zeros(a,a)

zb_im=zeros(a,a)
zb_im=complex(zb_im);

z_RE=ones(a,a)

q=collect(1:(a+1))

for g in q
        zb_re[ct4[g,2],ct4[g,3]]=ct4[g,4]

end

for l=1:a
    for m=1:a
        if zb_re[l,m]!==0
            zb_re[m,l]=zb_re[l,m]
        end
    end
end

for l=1:a
    zb_re[l,l]=-(zb_re[l,1]+zb_re[l,2]+zb_re[l,3]+zb_re[l,4]+zb_re[l,5])
end

for g in q
    zb_im[ct4[g,2],ct4[g,3]]=ct4[g,5]
end

for l=1:a
    for m=1:a
        if zb_im[l,m]!==0
            zb_im[m,l]=zb_im[l,m]

        end
    end
end

for l=1:a
    zb_im[l,l]=-(zb_im[l,1]+zb_im[l,2]+zb_im[l,3]+zb_im[l,4]+zb_im[l,5])
end

ZBUS=zb_re+(1im*zb_im)
YBUS=(ZBUS)^-1
