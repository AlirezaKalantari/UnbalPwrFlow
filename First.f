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


data = PowerModels.parse_file("C:/Users/kalan/Desktop/matpower7.0/data/case39.m")        #add case

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


q=collect(1:(size(ct4,1)))

#create Y bus
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

YBUS=yb_re+(1im*yb_im)      #Y buse
d=size(ct3,1)*ones(1,1)

#make unbalance ybus
#add time to any of P


for i=(1:size(ct3,1))
  b=size(ct3,1)
    if ct3[i,2]!==0.00
      d[1,1]=b-1
  end
end

real_pg=d[1,1]
ct6[2,2]=2

T=collect(1:24)             #duration
Iter=collect(1:1000)        #iteraion
G=collect(1:a)              #number of generator
G_1=G=collect(1:6*a)        #collect 1 to 6*size of buse
G_2=collect(1:3*a)          #collect 1 to 3*size of buse
G_3=collect(1:2*a)          #collect 1 to 2*size of buse
j=collect(1:(6*a))          #iteration of Jucobian
L=collect(1:a)              #number of load
Y=ones(6*a,6*a)             #Y buse
Ybus_un=ones(6*a,6*a)       #Y buse unbalance
J=zeros(6*a,6*a)            #Jucobian
E_k=complex(ones(3*a,1))    #voltage
I_cal=complex(ones(6*a,1))  #current
β_k=diagm(0=>ones(3*a))     #beta
α_k=diagm(0=>ones(3*a))     #alpha
γ_k=diagm(0=>ones(3*a))     #ghama
δ_k=(diagm(0=>ones(3*a)),24)#delta
delI_r=zeros(3*a,1)         # real of delta current
delI_m=zeros(3*a,1)         #imag of delta current
Pg=zeros(3*a,1)             #generation active power
Pl=zeros(3*a,1)             #active power of load
Pload_daily=Pl=zeros(3*a,24) #daily active power of load
P_sp=zeros(3*a,1)           #specified active power
P_cal=zeros(3*a,1)          #calculated active power
delP=zeros(3*a,1)           #active power mismatch
del_I=zeros(6*a,1)          #delta current
del_V=zeros(6*a,1)          #delta Voltage
newdel_V=zeros(6*a,1)       #voltage
del_v=zeros(6*a,1)          #delta voltage
V_rm=complex(ones(6*a,1))   #Voltage
solar_power=zeros(5,24)     #number of solar
data_demand=complex(zeros(5,24))   #demand of any time

Qg=zeros(3*a,1)             #generation of reactive power
Ql=zeros(3*a,1)             #reactive power of generation
Ql_daily=zeros(3*a,24)       # daily reactive power of generation
Ql_1=zeros(6*a,1)           #reserve of reactive power generated
Q_sp=zeros(3*a,1)           #specified reactive power
Q_cal=zeros(3*a,1)          #calculated reactive power
delQ=zeros(3*a,1)           #reactive power mismatch
tete=zeros(6*a,1)           #angle

iter=0                      #iteration cycle
iter=collect(1:10^1000)

for i=1:(size(ct5,1))       #active & reactive load of power
    Pl[ct5[i,2]]=ct5[i,3]
    Ql[ct5[i,2]]=ct5[i,4]
end

for i=1:a                   #active & reactive generation of power
    Pg[ct6[i,2]]=ct6[i,3]+solar_power[i,1]  #add solar power as generating power
    Qg[ct6[i,2]]=ct6[i,4]
end

P_sp=broadcast(-,Pg,Pl)     #equation 3
Q_sp=broadcast(-,Qg,Ql)     #equation 4

for j=1:1                   #number of slag bus
    h=0
    for i=1:a
        if ct1[2,i]==1
            h+=1
        end
    end
    println("Numer of S_Bus is =$h")
end

for j=1:1                   #number of PQ bus
    h=0
    for i=1:a
        if ct1[2,i]==2
            h+=1
        end
    end
    println("Numer of PQ_Bus is =$h")
end

for j=1:1                    #number of PV bus
    h=0
    for i=1:a
        if ct1[2,i]==3
            h+=1
        end
    end
    println("Numer of PV_Bus is =$h")
end


for t in t

    while iter==0                           #algorithm solution

        for i=1:(3*a)
            V_rm[i,1]=real(E_k)[i,1]        #equation2
            V_rm[i+3*a,1]=imag(E_k)[i,1]
        end

        for g in G_2
            P_cal[g]=(V_rm[g,1]'*real(I_cal[g])+V_rm[g+3*a,1]*imag(I_cal[g]))     #equation 22
            Q_cal[g]=(V_rm[g+3*a,1]'*real(I_cal[g])-V_rm[g,1]*imag(I_cal[g]))     #equation 23

            delP[g]=P_sp[g]-P_cal[g]                                              #equation 20
            delQ[g]=Q_sp[g]-Q_cal[g]                                              #equation 21
        end

        for g in G_2
            α_k[g,g]=(Q_sp[g]'*(V_rm[g,1]^2-V_rm[g+3*a,1]^2)-2*V_rm[g+3*a,1]*V_rm[g,1]*P_sp[g])/((V_rm[g+3*a,1]^2+V_rm[g,1]^2)^2)       #calculation of alpha
            β_k[g,g]=(P_sp[g]'*(V_rm[g,1]^2-V_rm[g+3*a,1]^2)+2*V_rm[g+3*a,1]*V_rm[g,1]*P_sp[g])/((V_rm[g+3*a,1]^2+V_rm[g,1]^2)^2)       #calculation of beta
            δ_k=α_k
            γ_k=-β_k
        end

        for g in G_2
            delI_r[g]=(delP[g]'*V_rm[g,1]+delQ[g]'*V_rm[g+3*a,1])/(V_rm[g+3*a,1]^2+V_rm[g,1]^2)             #equation 18
            delI_m[g]=(delP[g]'*V_rm[g+3*a,1]+delQ[g]'*V_rm[g,1])/(V_rm[g+3*a,1]^2+V_rm[g,1]^2)             #equation 19
        end

        for g in G_2
            del_I[(2*g-1),1]=delI_m[g]      #equation 18
            del_I[(2*g),1]=delI_r[g]        #equation 18
        end

        for b=1:a
            for w=1:3                                                   #equation 13
                if  Ql_1[6*(b-1)+w]!==0 || Ql_1[6*(b-1)+w+3]!==0
                    J[6*(b-1)+w,6*(b-1)+w]=imag(YBUS[b,b])-(α_k[3*(b-1)+w,3*(b-1)+w])       #equation 14
                    J[6*(b-1)+w+3,6*(b-1)+w+3]=-imag(YBUS[b,b])-(α_k[3*(b-1)+w,3*(b-1)+w])  #equation 17
                    J[6*(b-1)+w,6*(b-1)+w+3]=real(YBUS[b,b])-(β_k[3*(b-1)+w,3*(b-1)+w])     #equation 15
                    J[6*(b-1)+w+3,6*(b-1)+w]=real(YBUS[b,b])+(β_k[3*(b-1)+w,3*(b-1)+w])     #equation 16
                end
            end
        end

        for b=1:a
            for c=1:a               #equation 12
                for d=1:3
                    if b!==a &&  (Ql_1[6*(b-1)+d]!==0 || Ql_1[6*(b-1)+d+3]!==0)
                        J[6*(b-1)+d,6*(c-1)+d]=imag(YBUS[b,c])
                        J[6*(b-1)+d,6*(c-1)+d+3]=real(YBUS[b,c])
                        J[6*(b-1)+d+3,6*(c-1)+d]=real(YBUS[b,c])
                        J[6*(b-1)+d+3,6*(c-1)+d+3]=-imag(YBUS[b,c])
                    end
                end
            end
        end
        #check the Ql

        for b=1:a                            #calculation of Jucobian for pv buse
            if ct1[2,b]==3
                for c=1:a
                    if YBUS[c,b]!==0.0       #equation 28-29
                        if c==b
                            for d=1:3
                                J[6*(c-1)+d,6*(b-1)+d]=real(YBUS[c,b])-(β_k[3*(b-1)+d,3*(b-1)+d])-(imag(YBUS[c,b])-(α_k[3*(b-1)+d,3*(b-1)+d]))*(V_rm[b+3*a,1]/V_rm[c,1])
                                J[6*(c-1)+d,6*(b-1)+d+3]=V_rm[b,1]/(V_rm[b,1]^2+V_rm[b+3*a,1]^2)                                                                             #equation 30
                                J[6*(b-1)+d+3,6*(c-1)+d]=real(YBUS[c,b])+(β_k[3*(b-1)+d,3*(b-1)+d])+(imag(YBUS[b,b])+(α_k[3*(b-1)+d,3*(b-1)+d]))*(V_rm[b+3*a,1]/V_rm[c,1])
                                J[6*(b-1)+d+3,6*(c-1)+d+3]=-V_rm[b+3*a,1]/(V_rm[b,1]^2+V_rm[b+3*a,1]^2)                                                                      #equation 31
                            end
                        end
                    end
                end
            end
        end

        for b=1:a
            if ct1[2,b]==3
                for c=1:a
                    if YBUS[c,b]!==0.0
                        if c!==b
                            for d=1:3                    #equation 33-34
                                J[6*(c-1)+d,6*(b-1)+d]=real(YBUS[c,b])-(β_k[3*(b-1)+d,3*(b-1)+d])-(imag(YBUS[c,b])-(α_k[3*(b-1)+d,3*(b-1)+d]))*(V_rm[b+3*a,1]/V_rm[b,1])
                                J[6*(c-1)+d,6*(b-1)+d+3]=0
                                J[6*(b-1)+d+3,6*(c-1)+d]=real(YBUS[c,b])+(β_k[3*(b-1)+d,3*(b-1)+d])+(imag(YBUS[b,b])+(α_k[3*(b-1)+d,3*(b-1)+d]))*(V_rm[b+3*a,1]/V_rm[b,1])
                                J[6*(b-1)+d+3,6*(c-1)+d+3]=0
                                del_V[6*b+d,1]=del_V[d+3*b,1]
                                del_V[d+3*b,1]=delQ[3*(b-1)+d,1]
                            end
                        end
                    end
                end
            end
        end

        for g=1:(3*a)
           E_k[g]=V_rm[g,1]+((V_rm[g+3*a,1])*im)
           tete[g]=atand(V_rm[g+3*a,1]/V_rm[g,1])
        end

        J\del_I=del_v       #calculation of Jucobian for obtaining of Unknowns parameters


        newdel_V=del_v+del_V
        newdel_V=del_V

        for g in G_2        #constrain of voltage
            if -del_v[g]+del_V[g]<10^(-(10)^100)
             iter=1
            else
             iter=0
            end
        end
    end
end



#println("سلام سعید جان، وقت بخیر، من توضیحات کد رو نوشتم جلوی هر خط، لب تابمو رو دادم که درست کنن بلندگوشو تا باهمدیگه صحبت کنیم، امروز بعدازظهر دارم میرم کربلا، خوبی بدی دیدی حلال کن، اگه قابل باشم حتما برای خودت و خونواده ت دعا میکنم، انشااله برگشتم با قرار میذاریم، فقط سوالی که من در مورد کد داشتم این بود که چه طوری ولتاژهای نامتعادل رو وارد میکنیم توی کد. من کدهای تک تک معادلات یمقاله نوشتم، میخواستم خواهش کنم که اگه ممکنه منطق کدهارو هم یه نگاه کنی که درست کد زدم. خیلی مخلصیم. ارادتمند")


## Testing ##

#P_sp=zeros(3*a,24)           #specified active power
#P_cal=zeros(3*a,24)
#Q_cal=zeros(3*a,24)
#Q_sp=zeros(3*a,24)
#delP=zeros(3*a,24)
#delQ=zeros(3*a,24)
#Qg=zeros(3*a,24)
#Pg=zeros(3*a,24)
#Pl=zeros(3*a,24)
#Ql=zeros(3*a,24)
#delI_r=zeros(3*a,24)
#delI_m=zeros(3*a,24)
#I_cal=complex(ones(6*a,1),24)
#β_k=diagm(0=>ones(3*a))
#α_k=diagm(0=>ones(3*a))
#γ_k=diagm(0=>ones(3*a))
#δ_k=diagm(0=>ones(3*a))




#=for t=1:24
    for i=1:(size(ct5,1))       #active & reactive load of power of any time
        Pl[ct5[i,2],t]=real(data_demand[i,t])
        Ql[ct5[i,2],t]=imag(data_demand[i,t])
    end
end=#

#=for t=1:24
    for i=1:a                                     #active & reactive generation of power of any time
        Pg[ct6[i,2],t]=ct6[i,3]+solar_power[i,t]  #add solar power as generating power
        Qg[ct6[i,2],t]=ct6[i,4]
    end
end=#

#=for t=1:24
    for g in G_2
        P_cal[g,t]=(V_rm[g,1]'*real(I_cal[g])+V_rm[g+3*a,1]*imag(I_cal[g]))     #equation 22
        Q_cal[g,t]=(V_rm[g+3*a,1]'*real(I_cal[g])-V_rm[g,1]*imag(I_cal[g]))     #equation 23

        delP[g,t]=P_sp[g,t]-P_cal[g,t]                                          #equation 20
        delQ[g,t]=Q_sp[g,t]-Q_cal[g,t]                                          #equation 21
    end
end=#

#=for t=1:24
    for g in G_2
        delI_r[g,t]=(delP[g,t]'*V_rm[g,1]+delQ[g,t]'*V_rm[g+3*a,1])/(V_rm[g+3*a,1]^2+V_rm[g,1]^2)             #equation 18
        delI_m[g,t]=(delP[g,t]'*V_rm[g+3*a,1]+delQ[g,t]'*V_rm[g,1])/(V_rm[g+3*a,1]^2+V_rm[g,1]^2)             #equation 19
    end
end=#

#=for t=1:24
    for g in G_2
        del_I[(2*g-1),t]=delI_m[g,t]
        del_I[(2*g),t]=delI_r[g,t]
    end
end=#

#=for g in G_2
    α_k[g,g]=(Q_sp[g,t]'*(V_rm[g,1]^2-V_rm[g+3*a,1]^2)-2*V_rm[g+3*a,1]*V_rm[g,1]*P_sp[g,t])/((V_rm[g+3*a,1]^2+V_rm[g,1]^2)^2)       #calculation of alpha
    β_k[g,g]=(P_sp[g,t]'*(V_rm[g,1]^2-V_rm[g+3*a,1]^2)+2*V_rm[g+3*a,1]*V_rm[g,1]*P_sp[g,t])/((V_rm[g+3*a,1]^2+V_rm[g,1]^2)^2)       #calculation of beta
    δ_k=α_k
    γ_k=-β_k
end=#

#=for t=1:24
    for b=1:a
        for w=1:3                                                   #equation 13
            if  Ql_1[6*(b-1)+w,t]!==0 || Ql_1[6*(b-1)+w+3,t]!==0
                J[6*(b-1)+w,6*(b-1)+w]=imag(YBUS[b,b])-(α_k[3*(b-1)+w,3*(b-1)+w])       #equation 14
                J[6*(b-1)+w+3,6*(b-1)+w+3]=-imag(YBUS[b,b])-(α_k[3*(b-1)+w,3*(b-1)+w])  #equation 17
                J[6*(b-1)+w,6*(b-1)+w+3]=real(YBUS[b,b])-(β_k[3*(b-1)+w,3*(b-1)+w])     #equation 15
                J[6*(b-1)+w+3,6*(b-1)+w]=real(YBUS[b,b])+(β_k[3*(b-1)+w,3*(b-1)+w])     #equation 16
            end
        end
    end
end=#

#=for t=1:24
    for d=1:3
        if b!==a &&  (Ql_1[6*(b-1)+d]!==0 || Ql_1[6*(b-1)+d+3]!==0)
            J[6*(b-1)+d,6*(c-1)+d]=imag(YBUS[b,c])
            J[6*(b-1)+d,6*(c-1)+d+3]=real(YBUS[b,c])
            J[6*(b-1)+d+3,6*(c-1)+d]=real(YBUS[b,c])
            J[6*(b-1)+d+3,6*(c-1)+d+3]=-imag(YBUS[b,c])
        end
    end
end=#

#=P_sp=broadcast(-,Pg,Pl)
Q_sp=broadcast(-,Qg,Ql)=#
