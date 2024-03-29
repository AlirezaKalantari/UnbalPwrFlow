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

zb_re=zeros(5,5)

zb_im=zeros(5,5)
zb_im=complex(zb_im);

z_RE=ones(5,5)

q=collect(1:6)
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

Qg=zeros(size(G))
Ql=zeros(size(L))
Q_sp=zeros(size(G))
Q_cal=zeros(size(G))
delQ=zeros(size(G))
tete=zeros(4,1)  #tete=zeros(a,1)

#=for l in L
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

end=#
for i=1:3
    Pl[ct5[i,2]]=ct5[i,3]
    Ql[ct5[i,2]]=ct5[i,4]
end
for i=1:4
    Pg[ct6[i,2]]=ct6[i,3]
    Qg[ct6[i,2]]=ct6[i,4]
end

P_sp=broadcast(-,Pg,Pl)'
Q_sp=broadcast(-,Qg,Ql)'

#Ybus=Y
Ybus=YBUS

V_rk=real(E_k)
V_mk=imag(E_k)

G_y=real(Ybus)
B_y=imag(Ybus)

V_rm=[real(E_k)[1,1] real(E_k)[2,1] real(E_k)[3,1] real(E_k)[4,1] imag(E_k)[1,1] imag(E_k)[2,1] imag(E_k)[3,1] imag(E_k)[4,1]]'

#V_rm=V_new
#G=size(PowerModels.component_table(data, "bus", "bus_type"),1))
##V_rm[a+3]=V[a]
#for v_rm[..]=v_rm[a+3] and for v_rm[...]=v_new
# all V_rm[a+3] must be change to [a+4]
for a in G       #all [a+3] must change to [a+4] #Y_bus must change
    for b in G
        if a!==b && Ql[a]!==0
            c=2*b-1
            J[a,c]=imag(Ybus[a,b])-α_k[a,a]-((real(Ybus[a,b]))*(V_rm[a+3,1]/V_rm[a,1]))
            J[a,c+1]=real(Ybus[a,b])-β_k[a,a]
            a=a+1
            J[a,c]=imag(Ybus[a-1,b])+β_k[a-1,a-1]-((real(Ybus[a-1,b]))*(V_rm[a+3,1]/V_rm[a,1]))
            J[a,c+1]=real(Ybus[a-1,b])-α_k[a-1,a-1]
            a=a-1


        end
        if a==b
            J[2*a-1,b]=real(Ybus[a,a])-β_k[a,a]-((imag(Ybus[a,b])-α_k[a,a])*(V_rm[a+3,1]/V_rm[a,1]))
            J[2*a-1,b+1]=V_rm[a,1]/(V_rm[a,1]^2+V_rm[a+3,1]^2)
            a=a+1
            J[a,b]=real(Ybus[a-1,a-1])-α_k[a-1,a-1]+((imag(Ybus[a-1,b])+β_k[a-1,a-1])*(V_rm[a+3,1]/V_rm[a,1]))
            J[a,b+1]=-V_rm[a+3,1]/(V_rm[a,1]^2+V_rm[a+3,1]^2)
            a=a-1
            delI_m[a]=V_rm[a+3,1]/(V_rm[a,1]^2+V_rm[a+3,1]^2)
            delI_r[a]=V_rm[a,1]/(V_rm[a,1]^2+V_rm[a+3,1]^2)

        end
    end

end


if iter==0     #in this part too
    for g in G
        α_k[g,g]=(Q_sp[g]*(V_rm[g,1]^2-V_rm[g+3,1]^2)-2*V_rm[g+3,1]*V_rm[g,1]*P_sp[g])/((V_rm[g+3,1]^2+V_rm[g,1]^2)^2)
        β_k[g,g]=(P_sp[g]*(V_rm[g,1]^2-V_rm[g+3,1]^2)+2*V_rm[g+3,1]*V_rm[g,1]*P_sp[g])/((V_rm[g+3,1]^2+V_rm[g,1]^2)^2)
        δ_k=α_k
        γ_k=-β_k
    end
    for g in G
        P_cal[g]=V_rm[g,1]*real(I_cal)[g]'+V_rm[g+3,1]*imag(I_cal)[g]'
        Q_cal[g]=V_rm[g+3,1]*real(I_cal)[g]'-V_rm[g,1]*imag(I_cal)[g]'

        delP[g]=P_sp[g]-P_cal[g]
        delQ[g]=Q_sp[g]-Q_cal[g]


        delI_r[g]=(delP[g]*V_rm[g,1]+V_rm[g+3,1]*delQ[g])/(V_rm[g+3,1]^2+V_rm[g,1]^2)
        delI_m[g]=(delP[g]*V_rm[g+3,1]+V_rm[g,1]*delQ[g])/(V_rm[g+3,1]^2+V_rm[g,1]^2)

        del_I[(2*g-1),1]=delI_m[g]
        del_I[(2*g),1]=delI_r[g]



    end
    for a in G      #J must be 9*a*a
        for b in G
            if a==b

                J[2*a-1,b]=imag(Ybus[a,a])-α_k[a,a]
                J[2*a-1,b+1]=real(Ybus[a,a])-β_k[a,a]
                a=a+1
                J[a,b]=real(Ybus[a-1,a-1])-α_k[a-1,a-1]
                J[a,b+1]=-imag(Ybus[a-1,a-1])+β_k[a-1,a-1]
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
    #    J[a,c]=imag(Ybus[a,b])-α_k[a,a]-((real(Ybus[a,b])*(V_rm[a+3,1]/V_rm[a,1])))
    #    J[a,c+3]=real(Ybus[a,b])-β_k[a,a]
    #    a=a+1
    #    J[a,c]=imag(Ybus[a-1,b])+β_k[a-1,a-1]-((real(Ybus[a-1,b])*(V_rm[a+3,1]/V_rm[a+3])))
    #    J[a,c+3]=real(Ybus[a-1,b])-α_k[a-1,a-1]
    #    a=a-1
    #if a==b
    #    J[2*a-1,b]=real(Ybus[a,a])-β_k[a,a]-((imag(Ybus[a,b])-α_k[a,a])*(V_rm[a+3,1]/V_rm[a]))
    #    J[2*a-1,b+3]=V_rm[a,1]/(V_rm[a,1]^2+V_rm[a+3,1]^2)
    #    a=a+1
    #    J[a,b]=J[2*a-1,b]=real(Ybus[a-1,a-1])-α_k[a-1,a-1]+((imag(Ybus[a-1,b])+β_k[a-1,a-1])*(V_rm[a+1,1]/V_rm[a+1,1]))
    #    J[a,b+3]=-V_rm[a+3,1]/(V_rm[a,1]^2+V_rm[a+3,1]^2)
    #    a=a-1
    #    delI_m[a]=V_rm[a+3,1]/(V_rm[a,1]^2+V_rm[a+3,1]^2)
    #    delI_r[a]=V_rm[a,1]/(V_rm[a,1]^2+V_rm[a+3,1]^2)
    #    del_I[(2*a-1),1]=delI_m[g]
    #    del_I[(2*a),1]=delI_r[g]


    #end
    #end

    #=  for g in G
           α_k[g,g]=(Q_sp[g]*(V_rm[g,1]^2-V_rm[g+3,1]^2)-2*V_rm[g+3,1]*V_rm[g,1]*P_sp[g])/((V_rm[g+3,1]^2+V_rm[g,1]^2)^2)
           β_k[g,g]=(Q_sp[g]*(V_rm[g,1]^2-V_rm[g+3,1]^2)+2*V_rm[g+3,1]*V_rm[g,1]*P_sp[g])/((V_rm[g+3,1]^2+V_rm[g,1]^2)^2)
           δ_k=α_k
           γ_k=-β_k


       end
        for g in G
           P_cal[g]=V_rm[g,1]*real(I_cal)[g]'+V_rm[g+3,1]*imag(I_cal)[g]'
           Q_cal[g]=V_rm[g+3,1]*real(I_cal)[g]'-V_rm[g,1]*imag(I_cal)[g]'

           delP[g]=P_sp[g]-P_cal[g]
           delQ[g]=Q_sp[g]-Q_cal[g]
           delI_r[g]=(delP[g]*V_rm[g,1]+V_rm[g+3,1]*delQ[g])/(V_rm[g+3,1]^2+V_rm[g,1]^2)
           delI_m[g]=(delP[g]*V_rm[g+3,1]+V_rm[g,1]*delQ[g])/(V_rm[g+3,1]^2+V_rm[g,1]^2)
           del_I[(2*g-1),1]=delI_m[g]
           del_I[(2*g),1]=delI_r[g]
       end
        for a in G
            for b in G
                if a==b
                   J[2*a-1,b]=imag(Ybus[a,a])-α_k[a,a]
                   J[2*a-1,b+1]=real(Ybus[a,a])-β_k[a,a]
                   a=a+1
                   J[a,b]=real(Ybus[a-1,a-1])-α_k[a-1,a-1]
                   J[a,b+1]=-imag(Ybus[a-1,a-1])+β_k[a-1,a-1]
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

    end=#
    #=for a in G
        for b in G
            if a!==b && Ql[a]!==0
                c=2*b-1
                J[a,c]=imag(Ybus[a,b])-α_k[a,a]-((real(Ybus[a,b]))*(V_rm[a+3,1]/V_rm[a,1]))
                J[a,c+1]=real(Ybus[a,b])-β_k[a,a]
                a=a+1
                J[a,c]=imag(Ybus[a-1,b])+β_k[a-1,a-1]-((real(Ybus[a-1,b]))*(V_rm[a+1,1]/V_rm[a-1,1]))
                J[a,c+1]=real(Ybus[a-1,b])-α_k[a-1,a-1]
                a=a-1

            end
            if a==b
                J[2*a-1,b]=real(Ybus[a,a])-β_k[a,a]-((imag(Ybus[a,b])-α_k[a,a])*(V_rm[a+3,1]/V_rm[a,1]))
                J[2*a-1,b+1]=V_rm[a,1]/(V_rm[a,1]^2+V_rm[a+3,1]^2)
                a=a+1
                J[a,b]=real(Ybus[a-1,a-1])-α_k[a-1,a-1]+((imag(Ybus[a-1,b])+β_k[a-1,a-1])*(V_rm[a+1,1]/V_rm[a-1,1]))
                J[a,b+1]=-V_rm[a+3,1]/(V_rm[a,1]^2+V_rm[a+3,1]^2)
                a=a-1
                delI_m[a]=V_rm[a+3,1]/(V_rm[a,1]^2+V_rm[a+3,1]^2)
                delI_r[a]=V_rm[a,1]/(V_rm[a,1]^2+V_rm[a+3,1]^2)

            end
        end

    end=#
 del_v\del_I
    #=for g in G
        del_V[2*g-1,1]=V_rm[g,1]
        del_V[2*g+2,1]=V_rm_k[g+3,1]
        #=for i 1:3
            del_V[i,1]=V_rm[i,1]
            del_V[i+3,1]=V_rm_k[i+3,1]
        end=#
    end=#

    #=
    (PowerModels.component_table(data, "bus", "bus_type"),1)
    del_V=zeros(3*a,1)
    for i 1:a
         for j 1:6
             del_V[6*(i-1)+j,1]=V_new[j,1]
         end
    end=#

    newdel_V=del_v+del_V
    newdel_V=del_V
    for g in G
        if del_v[g]-del_V[g]<10^(-(10)^100)
         iter=1
        end
    end
    for g=1:4
       E_k[g]=V_rm[g,1]+((V_rm[g+4,1])*im)
       tete[g]=atand(V_rm[g+4,1]/V_rm[g,1])
    end
end
#a=size(PowerModels.component_table(data, "bus", "bus_type"),1)
# V_new=complex(zeros(24,1))

#V_rm_k[i]=[V_rka V_rkb V_rkc V_mka V_mkb V_mkc]
#V_I=zeros(6,1)
#VRM=zeros(6,1)
#=for u 1:a
    VRM[u,1]=[V_rka V_rkb V_rkc V_mka V_mkb V_mkc]
end=#

#=for i 1:a
    #=for i 1:6
        V_I[i,1]=V_rm_k[i,1]
    end
    for j 1:3
        V_new[6*(i-1)+j,1]=V_I[j,1]
    end
    for k 1:3
        V_new[6*i+k-3,1]=V_I[k+3,1]
    end
    for l 1:6
        V_rm_k[l,1]=VRM[l,1]
    end

end=#
#V_ra=[V_rka V_rkb V_rkc]
#V_ma=[V_mka V_mkb V_mkc]
#V_rm=[V_rka V_rkb V_rkc V_mka V_mkb V_mkc]'
#o[l,l]=zeros(3,3)
#p[l,l]=zeros(3,3)
#for i 1:3
#    o[l,l]=v_rm[l,1]/(V_rm[l,1]^2+V_rm[l+3,1]^2)
#    p[l,l]=-v_rm[l,l]/(V_rm[l,1]^2+V_rm[l+3,1]^2)
#end
###for i 1:4
###    o[l,l]=v_rm[l,1]/(V_rm[l,1]^2+V_rm[l+4,1]^2)
###   p[l,l]=-v_rm[l,l]/(V_rm[l,1]^2+V_rm[l+4,1]^2)
###end


#ABC TYPE
#V_rm=V_new=#

#=a=size(PowerModels.component_table(data, "bus", "bus_type"),1)
del_V=zeros(3*a,1)
for i 1:a
     for j 1:6
         del_V[6*(i-1)+j,1]=V_new[j,1]
     end
end=#
#V_rm=V_new


println("salam saeed jan, I travel to Gorgan this Afternoon and start it the day ahead")
