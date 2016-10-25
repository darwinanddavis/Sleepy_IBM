# note this is all for 35 degrees and needs to be corrected to whatever the 2min body temp is
Wdry=251.8 #g
Wwet=Wdry/0.3 #g
Lmax=7.929 #cm
Vmax=Lmax^3 #cm3
p_M.spec=31.8 #J/cm3/d
p_M=p_M.spec*Vmax #J/d
p_M.h=p_M/24 #J/h
mrate.watts=p_M/(3600*24) #J/s
mrate.O2.h=mrate.watts/20.1*3600 # convert to ml O2, make it per hour
mrate.O2.g.h=mrate.O2.h/Wwet
VO2max=0.722 #ml O2/g/h
J.Loc.h=VO2max*Wwet*20.1-p_M.h #J/h
J.Loc.2min=J.Loc.h/30 # J/min activity, extra activity cost to add when moving
# last step - adjust J.Loc.2min to what it would be fore 20 degrees, using the arrhenius correction factor

-------------------------------------------
  # body temperature dependent direct movement cost (cost of locomotion)  
  #initial par values
  
 step<-2/1440 # step size (2 mins). For hourly: 1/24
p_M<-32 * step #J

T_REF = 20
TA = 8085
TAL = 18721
TAH = 9e+4
TL = 288 
TH = 315
Tb=seq(0,45)
Tcorr = exp(TA * (1/(273 + T_REF) - 1/(273 + Tb)))/(1 + exp(TAL * (1/(273 + Tb) - 1/TL)) + exp(TAH * (1/TH - 1/(273 + Tb))))
plot(Tb,Tcorr)

# correct [p_M] from 20 to 35 degrees
Tb=35
Tcorr = exp(TA * (1/(273 + T_REF) - 1/(273 + Tb)))/(1 + exp(TAL * (1/(273 + Tb) - 1/TL)) + exp(TAH * (1/TH - 1/(273 + Tb))))
p_M=p_M*Tcorr # now corrected to 25 deg C

L_max<-7.997 #cm
V_max<-L_max^3 #cm3
p_M<-p_M*V_max #J/d
p_M.h<-p_M/24 #J/h
mrate.watts<-p_M/(3600*24) #J/s
mrate.O2.h<-mrate.watts/20.1*3600 # convert to ml O2, make it per hour
mrate.O2.g.h<-mrate.O2.h/mass # ml 02 per h by total wetmass
VO2max<-0.722 #ml O2/g/h

T_REF = 35
Tb=20
Tcorr = exp(TA * (1/(273 + T_REF) - 1/(273 + Tb)))/(1 + exp(TAL * (1/(273 + Tb) - 1/TL)) + exp(TAH * (1/TH - 1/(273 + Tb))))

VO2max<-VO2max*Tcorr
J.Loc.h<-VO2max*mass*20.1-p_M.h #J/h
J.Loc.2min<-J.Loc.h/30 # J/2min activity, extra activity cost to add when moving

#28-7-16
# direct movement cost. convert o2 to grams then multiple by mass
#loco_km<-0.921/20.1*800 # J/km
#loco<-loco_km/1000 #J/m

#-------------------------------------------



# ------- New loco cost 16-10-16 -------
mass<-800 # mass
V_pres = 3.9752^3 #structure 
step<-1/24 #hourly
#step<-2/1440 #2min

p_M<-32*step #J/h
p_M<-p_M*V_pres # loco cost * structure

# movement cost for time period
VO2<-0.45 # O2/g/h JohnAdler etal 1986

# multiple p_M by structure = movement cost (diff between p_M with loco cost and structure for movement period)
# p_M with loco cost 
loco<-VO2*mass*20.1 # convert ml O2 to J = J/h 
loco<-loco+p_M # add to p_M = J/h
loco<-loco/30 ; loco #J/2min


