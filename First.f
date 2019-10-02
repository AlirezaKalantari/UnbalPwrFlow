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

#PV_bus=zeros(1,1)
Iter=collect(1:1000)
G=collect(1:a)    #number of generator
G_1=G=collect(1:6*a)
G_2=collect(1:3*a)
G_3=collect(1:2*a)
j=collect(1:(6*a))
L=collect(1:a)    #number of load
Y=ones(6*a,6*a)
J=zeros(6*a,6*a)
E_k=complex(ones(3*a,1))
I_cal=complex(ones(6*a,1))
β_k=diagm(0=>ones(3*a))
α_k=diagm(0=>ones(3*a))
γ_k=diagm(0=>ones(3*a))
δ_k=diagm(0=>ones(3*a))
delI_r=zeros(3*a,1)
delI_m=zeros(3*a,1)
Pg=zeros(3*a,1)
Pl=zeros(3*a,1)
P_sp=zeros(3*a,1)
P_cal=zeros(3*a,1)
delP=zeros(3*a,1)
del_I=zeros(6*a,1)
del_V=zeros(6*a,1)
newdel_V=zeros(6*a,1)
del_v=zeros(6*a,1)
V_rm=complex(ones(6*a,1))

Qg=zeros(3*a,1)
Ql=zeros(3*a,1)
Ql_1=zeros(6*a,1)
Q_sp=zeros(3*a,1)
Q_cal=zeros(3*a,1)
delQ=zeros(3*a,1)
tete=zeros(6*a,1)

iter=0
iter=collect(1:10^1000)

for i=1:(size(ct5,1))
    Pl[ct5[i,2]]=ct5[i,3]
    Ql[ct5[i,2]]=ct5[i,4]
end

ct6[2,2]=2

for i=1:a
    Pg[ct6[i,2]]=ct6[i,3]
    Qg[ct6[i,2]]=ct6[i,4]
end

P_sp=broadcast(-,Pg,Pl)
Q_sp=broadcast(-,Qg,Ql)

for j=1:1
    h=0
    for i=1:a
        if ct1[2,i]==1
            h+=1
        end
    end
    println("Numer of S_Bus is =$h")
end

for j=1:1
    h=0
    for i=1:a
        if ct1[2,i]==2
            h+=1
        end
    end
    println("Numer of PQ_Bus is =$h")
end

for j=1:1
    h=0
    for i=1:a
        if ct1[2,i]==3
            h+=1
        end
    end
    println("Numer of PV_Bus is =$h")
end






if iter==0

    for i=1:(3*a)
        V_rm[i,1]=real(E_k)[i,1]
        V_rm[i+3*a,1]=imag(E_k)[i,1]
    end

    for g in G_2
        P_cal[g]=(V_rm[g,1]'*real(I_cal)[g])+(V_rm[g+3*a,1]*imag(I_cal)[g])
        Q_cal[g]=(V_rm[g+3*a,1]'*real(I_cal)[g])-(V_rm[g,1]*imag(I_cal)[g])

        delP[g]=P_sp[g]-P_cal[g]
        delQ[g]=Q_sp[g]-Q_cal[g]
    end

    for g in G_2
        α_k[g,g]=(Q_sp[g]'*(V_rm[g,1]^2-V_rm[g+3*a,1]^2)-2*V_rm[g+3*a,1]*V_rm[g,1]*P_sp[g])/((V_rm[g+3*a,1]^2+V_rm[g,1]^2)^2)
        β_k[g,g]=(P_sp[g]'*(V_rm[g,1]^2-V_rm[g+3*a,1]^2)+2*V_rm[g+3*a,1]*V_rm[g,1]*P_sp[g])/((V_rm[g+3*a,1]^2+V_rm[g,1]^2)^2)
        δ_k=α_k
        γ_k=-β_k
    end

    for g in G_2
        delI_r[g]=(delP[g]'*V_rm[g,1]+delQ[g]'*V_rm[g+3*a,1])/(V_rm[g+3*a,1]^2+V_rm[g,1]^2)
        delI_m[g]=(delP[g]'*V_rm[g+3*a,1]+delQ[g]'*V_rm[g,1])/(V_rm[g+3*a,1]^2+V_rm[g,1]^2)
    end

    for g in G_2
        del_I[(2*g-1),1]=delI_m[g]
        del_I[(2*g),1]=delI_r[g]
    end

    for b=1:a
        for w=1:3
            if  Ql_1[6*(b-1)+w]!==0 || Ql_1[6*(b-1)+w+3]!==0
                J[6*(b-1)+w,6*(b-1)+w]=imag(YBUS[b,b])-(α_k[3*(b-1)+w,3*(b-1)+w])
                J[6*(b-1)+w+3,6*(b-1)+w+3]=-imag(YBUS[b,b])-(α_k[3*(b-1)+w,3*(b-1)+w])
                J[6*(b-1)+w,6*(b-1)+w+3]=real(YBUS[b,b])-(β_k[3*(b-1)+w,3*(b-1)+w])
                J[6*(b-1)+w+3,6*(b-1)+w]=real(YBUS[b,b])+(β_k[3*(b-1)+w,3*(b-1)+w]
                #J[6*(b-1)+w+6,6*(b-1)+w+3]=-V_rm[b+3*a,1]/(V_rm[b,1]^2+V_rm[b+3*a,1]^2)
                #J[6*(b-1)+w+3,6*(b-1)+w+6]=-V_rm[b,1]/(V_rm[b,1]^2+V_rm[b+3*a,1]^2)
                #J[s,b]=real(YBUS[w,w])-(β_k[w,w])
                #J[s,b+1]=-imag(YBUS[w,w])+(α_k[w,w])
                #=elseif s!==b  && Ql[s]!==0
                c=2*b-1
                J[s,c]=imag(YBUS[w,b])
                J[s,c+1]=real(YBUS[w,b])
                #s=s+1
                J[s,c]=real(YBUS[w,b])
                J[s,c+1]=-imag(YBUS[w,b])=#
                #s=s-1
            end
        end
    end

    for b=1:a
        for c=1:a
            for d=1:3
                if b!==a &&  (Ql_1[6*(b-1)+b]!==0 || Ql_1[6*(b-1)+b+3]!==0)
                    J[6*(b-1)+d,6*(c-1)+d]=imag(YBUS[b,c])
                    J[6*(b-1)+d,6*(c-1)+d+3]=real(YBUS[b,c])
                    J[6*(b-1)+d+3,6*(c-1)+d]=real(YBUS[b,c])
                    J[6*(b-1)+d+3,6*(c-1)+d+3]=-imag(YBUS[b,c])
                end
            end
        end
    end
 #check the Ql


    for b=1:a
        if ct1[2,b]==3
            for c=1:a
                if YBUS[c,b]!==0.0
                    if c==b
                        for d=1:3
                            J[6*(c-1)+d,6*(b-1)+d]=real(YBUS[c,b])-(β_k[3*(b-1)+d,3*(b-1)+d])-(imag(YBUS[c,b])-(α_k[3*(b-1)+d,3*(b-1)+d]))*(V_rm[b+3*a,1]/V_rm[c,1])
                            J[6*(c-1)+d,6*(b-1)+d+3]=V_rm[b,1]/(V_rmbs,1]^2+V_rm[b+3*a,1]^2)
                            #J[6*(b-1)+d,6*(c-1)+d+3]=real(YBUS[b,c])
                            #J[6*(b-1)+d+3,6*(c-1)+d]=imag(YBUS[b,c])
                            #J[6*(b-1)+d+3,6*(c-1)+d+3]=imag(YBUS[b,c])
                        end
                    end
                end
            end
        end
    end




    #=for s in G_1
        for b in G_2
            if s!==b && Ql[s]!==0
                c=2*b-1
                J[s,c]=imag(YBUS[s,b])-α_k[s,s]-((real(YBUS[s,b]))*(V_rm[s+3*a,1]/V_rm[s,1]))
                J[s,c+1]=real(YBUS[s,b])-β_k[s,s]
                s=s+1
                J[s,c]=imag(YBUS[s-1,b])+β_k[s-1,s-1]-((real(YBUS[s-1,b]))*(V_rm[s+3*a,1]/V_rm[s,1]))
                J[s,c+1]=real(YBUS[s-1,b])-α_k[s-1,s-1]
                s=s-1
            end
            if  s==b
                J[2*s-1,b]=real(YBUS[s,s])-β_k[s,s]-((imag(YBUS[s,b])-α_k[s,s])*(V_rm[s+3*a,1]/V_rm[s,1]))
                J[2*s-1,b+1]=V_rm[s,1]/(V_rm[s,1]^2+V_rm[s+3*a,1]^2)
                s=s+1
            end
        end
    end=#



    #=for s in G_1
       for b in G_2
           if s!==b && Ql[s]!==0
               c=2*b-1
               J[s,c]=imag(YBUS[s,b])-α_k[s,s]-((real(YBUS[s,b])*(V_rm[s+3*a,1]/V_rm[s,1])))
               J[s,c+3]=real(YBUS[s,b])-β_k[s,s]
               s=s+1
               J[s,c]=imag(YBUS[s-1,b])+β_k[s-1,s-1]-((real(YBUS[s-1,b])*(V_rm[s+3*a,1]/V_rm[s+3*a])))
               J[s,c+3]=real(YBUS[s-1,b])-α_k[s-1,s-1]
               s=s-1
               #delI_m[s]=V_rm[s+3*a,1]/(V_rm[s,1]^2+V_rm[s+3*a,1]^2)
               #delI_r[s]=V_rm[s,1]/(V_rm[s,1]^2+V_rm[s+3*a,1]^2)

               #del_I[(2*b-1),1]=delI_m[s]
               #del_I[(2*b),1]=delI_r[s]
           end
       end
   end=#

    for g=1:(3*a)
       E_k[g]=V_rm[g,1]+((V_rm[g+3*a,1])*im)
       tete[g]=atand(V_rm[g+3*a,1]/V_rm[g,1])
    end

    del_v\del_I

    newdel_V=del_v+del_V
    newdel_V=del_V

    for g in G_2
        if del_v[g]-del_V[g]<10^(-(10)^100)
         iter=1
        end
    end

end

#PV_bus
#=for s in G_1
   for b in G_2
       if s!==b && Ql[s]!==0
           c=2*b-1
           J[s,c]=imag(YBUS[s,b])-α_k[s,s]-((real(YBUS[s,b])*(V_rm[s+3*a,1]/V_rm[s,1])))
           J[s,c+3]=real(YBUS[s,b])-β_k[s,s]
           s=s+1
           J[s,c]=imag(YBUS[s-1,b])+β_k[s-1,s-1]-((real(YBUS[s-1,b])*(V_rm[s+3*a,1]/V_rm[s+3*a])))
           J[s,c+3]=real(YBUS[s-1,b])-α_k[s-1,s-1]
           s=s-1
       end
        #=if s==b
           J[2*s-1,b]=real(YBUS[s,s])-β_k[s,s]-((imag(YBUS[s,b])-α_k[s,s])*(V_rm[s+3,1]/V_rm[s]))
           J[2*s-1,b+3]=V_rm[s,1]/(V_rm[s,1]^2+V_rm[s+3*a,1]^2)
           s=s+1
           J[s,b]=J[2*s-1,b]=real(YBUS[s-1,s-1])-α_k[s-1,s-1]+((imag(YBUS[s-1,b])+β_k[s-1,s-1])*(V_rm[s+1,1]/V_rm[s+1,1]))
           J[s,b+3]=-V_rm[s+3*a,1]/(V_rm[s,1]^2+V_rm[s+3*a,1]^2)
           s=s-1
       end=#
   end
end=#

println("salam saeed jan, because my mother is in the hospital, I have to take care of her, when I return to Tehran, we set a meeting ")




#=for s in G_1
    for b in G_2
        for w=1:a
            if s==b
                J[2*s-1,b]=imag(YBUS[w,w])-(α_k[w,w])
                J[2*s-1,b+1]=real(YBUS[w,w])-(β_k[w,w])
                #s=s+1
                J[s,b]=real(YBUS[w,w])-α_k[w,w]
                J[s,b+1]=-imag(YBUS[w,w])+β_k[w,w]
                #s=s-1
                elseif s!==b
                c=2*b-1
                J[s,c]=imag(YBUS[w,w])
                J[s,c+1]=real(YBUS[w,w])
                #s=s+1
                J[s,c]=real(YBUS[w,w])
                J[s,c+1]=-imag(YBUS[w,w])
                #s=s-1
            end
        end
    end
end=#
