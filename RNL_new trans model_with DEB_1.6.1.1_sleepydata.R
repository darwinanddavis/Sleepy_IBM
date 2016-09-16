

#*******************************************************************
# ****************** TRANSIENT MODEL SETUP *************************

# Transient model, interpolated to particular time of day from microclimate model prediction
# RNL_new trans model_with DEB_1.4


# Transient model function set up for a constant environment and a starting condition
# Michael Kearney & Warren Porter developed this R function on 28 July 2014.
# Michael modified it on 1st Aug to work as an analytical model, without evaporation

transient<-function(t,y,thresh,input){
  unlist(input)
  sigma<-0.0000000567 #Stefan-Boltzman, W/(m.K)
  Zenith<-Zen*pi/180 # zenith angle in radians
  Tc<-y # core temperature, deg C
  Tskin<-y+0.1 # make skin temperature very close to core temperature
  if(vel<0.01){
    vel<-0.01 # don't let wind speed go too low - always some free convection
  }
  S2<-0.0001 # initializing
  DENSTY<-press/(287.04*(Tair+273)) # air density, kg/m3
  THCOND<-0.02425+(7.038*10^-5*Tair) # air thermal conductivity, W/(m.K)
  VISDYN<-(1.8325*10^-5*((296.16+120)/((Tair+273)+120)))*(((Tair+273)/296.16)^1.5) # dynamic viscosity of air, kg/(m.s)

  m<-mass/1000 # convert mass to kg
  C<-m*cp # thermal capacitance, J/K
  V<-m/rho # volume, m3
  Qgen<-q*V # total metabolic heat, J
  L<-V^(1./3.) # characteristic dimension, m
  
  # geometry section ############################################################
  # FLAT PLATE geometry
  if(lometry==0){
    ALENTH<-(V/shape_b*shape_c)^(1./3.) # length, m
    AWIDTH<-ALENTH*shape_b # width, m
    AHEIT<-ALENTH*shape_c # height, m    
    ATOT<-ALENTH*AWIDTH*2.+ALENTH*AHEIT*2.+AWIDTH*AHEIT*2. # total area, m2
    ASILN<-ALENTH*AWIDTH # max silhouette area, m2
    ASILP<-AWIDTH*AHEIT # min silhouette area, m2
    L<-AHEIT # characteristic dimension, m
    if(AWIDTH<=ALENTH){
      L<-AWIDTH
    }else{
      L<-ALENTH
    }
    R<-ALENTH/2. # 'radius', m
  } 
  
  # CYLINDER geometry      
  if(lometry==1){
    R1<-(V/(pi*shape_b*2))^(1./3.) # radius, m
    ALENTH<-2*R1*shape_b # length, m
    ATOT<-2*pi*R1^2+2*pi*R1*ALENTH # total surface area, m2
    AWIDTH<-2.*R1 # width, m
    ASILN<-AWIDTH*ALENTH # max silhouette area, m2
    ASILP<-pi*R1^2 # min silhouette area, m2
    L<-ALENTH # characteristic dimension, m
    R2<-L/2
    if(R1>R2){ # choose shortest dimension as R
      R<-R2
    }else{
      R<-R1
    }
  }
  
  # Ellipsoid geometry
  if(lometry==2){
    A1<-((3./4.)*V/(pi*shape_b*shape_c))^0.333 # axis A, m  
    B1<-A1*shape_b # axis B, m
    C1<-A1*shape_c # axis C, m
    P1<-1.6075 # a constant
    ATOT<-(4*pi*(((A1^P1*B1^P1+A1^P1*C1^P1+B1^P1*C1^P1))/3)^(1/P1)) # total surface area, m2
    ASILN<-max(pi*A1*C1,pi*B1*C1) # max silhouette area, m2
    ASILP<-min(pi*A1*C1,pi*B1*C1) # min silhouette area, m2
    S2<-(A1^2*B1^2*C1^2)/(A1^2*B1^2+A1^2*C1^2+B1^2*C1^2) # fraction of semi-major and minor axes, see Porter and Kearney 2009 supp1
    kflesh<-0.5 + 6.14*B1 + 0.439 # thermal conductivity of flesh as a function of radius, see Porter and Kearney 2009
  }              
  
  # Lizard geometry - DESERT IGUANA (PORTER ET AL. 1973 OECOLOGIA)
  if(lometry==3){
    ATOT<-(10.4713*mass^.688)/10000. # total surface area, m2
    AV<-(0.425*mass^.85)/10000. # ventral surface area, m2   
    # NORMAL AND POINTING @ SUN SILHOUETTE AREA: PORTER & TRACY 1984   
    ASILN<-(3.798*mass^.683)/10000. # Max. silhouette area (normal to the sun), m2
    ASILP<-(0.694*mass^.743)/10000. # Min. silhouette area (pointing toward the sun), m2
    R<-L
  }
  
  # Frog geometry - LEOPARD FROG (C.R. TRACY 1976 ECOL. MONOG.)
  if(lometry==4){
    ATOT = (12.79*mass^.606)/10000. # total surface area, m2
    AV = (0.425*mass^.85)/10000. # ventral surface area, m2  
    # NORMAL AND POINTING @ SUN SILHOUETTE AREA: EQ'N 11 TRACY 1976
    ZEN<-0.
    PCTN<-1.38171E-06*ZEN^4-1.93335E-04*ZEN^3+4.75761E-03*ZEN^2-0.167912*ZEN+45.8228  
    ASILN<-PCTN*ATOT/100. # Max. silhouette area (normal to the sun), m2
    ZEN<-90. 
    PCTP<-1.38171E-06*ZEN^4-1.93335E-04*ZEN^3+4.75761E-03*ZEN^2-0.167912*ZEN+45.8228  
    ASILP<-PCTP*ATOT/100. # Min. silhouette area (pointing toward the sun), m2
    R<-L
  }
  
  # user defined geometry
  if(lometry==5){
    ATOT = (customallom[1]*mass^customallom[2])/10000. # total surface area, m2 
    AV = (customallom[3]*mass^customallom[4])/10000. # ventral surface area, m2   
    # NORMAL AND POINTING @ SUN SILHOUETTE AREA: PORTER & TRACY 1984   
    # User must define Max. silhouette area (normal to the sun)
    ASILN = (customallom[5]*mass^customallom[6])/10000. # Max. silhouette area (normal to the sun), m2
    # User must define Min. silhouette area (pointing toward the sun)         
    ASILP = (customallom[7]*mass^customallom[8])/10000. # Min. silhouette area (pointing toward the sun), m2
    R<-L
  }
  # end geometry section ############################################################
   
  if(Zen>=90){
    Qnorm<-0
  }else{
    Qnorm <- (Qsol / cos(Zenith)) 
  }
  if(Qnorm>1367){
    Qnorm<-1367 #making sure that low sun angles don't lead to solar values greater than the solar constant
  }
  if(posture=='p'){
    Qabs<-(Qnorm*(1-pctdif)*ASILP+Qsol*pctdif*FATOSK*ATOT+Qsol*sub_reflect*FATOSB*ATOT)*abs
  }
  if(posture=='n'){
    Qabs<-(Qnorm*(1-pctdif)*ASILN+Qsol*pctdif*FATOSK*ATOT+Qsol*sub_reflect*FATOSB*ATOT)*abs
  }
  if(posture=='b'){
    Qabs<-(Qnorm*(1-pctdif)*(ASILN+ASILP)/2+Qsol*pctdif*FATOSK*ATOT+Qsol*sub_reflect*FATOSB*ATOT)*abs
  }
  
  Rrad<-((Tskin+273)-(Trad+273))/(emis*sigma*Fo_e*ATOT*((Tskin+273)^4-(Trad+273)^4)) # radiation resistance
  Re<-DENSTY*vel*L/VISDYN # Reynolds number
  PR<-1005.7*VISDYN/THCOND # Prandlt number
  
  if(lometry==0){
    NUfor<-0.102*Re^0.675*PR^(1./3.)
  }
  if(lometry==3|lometry ==5){
    NUfor<-0.35*Re^0.6
  }
  if(lometry==1){ 
    #       FORCED CONVECTION OF A CYLINDER    
    #       ADJUSTING NU - RE CORRELATION FOR RE NUMBER (P. 260 MCADAMS,1954) 
    if(Re<4.){ 
      NUfor=.891*Re**.33  
    }else{
      if(Re<40.){
        NUfor=.821*Re**.385 
      }else{
        if(Re<4000.){
          NUfor=.615*Re**.466 
        }else{ 
          if(Re<40000.){  
            NUfor=.174*Re**.618 
          }else{
            if(Re<400000.){ 
              NUfor=.0239*Re**.805
            }else{
              NUfor=.0239*Re**.805
            }}}}} 
  }
  if(lometry==2|lometry==4){
    NUfor<-0.35*Re^(0.6) # Nusselt number, forced convection
  }
  hc_forced<-NUfor*THCOND/L # convection coefficent, forced
  
  GR<-abs(DENSTY^2*(1/(Tair+273.15))*9.80665*L^3*(Tskin-Tair)/VISDYN^2) # Grashof number
  Raylei<-GR*PR # Rayleigh number
  
  # get Nusselt for Free Convect
  if(lometry==0){
    NUfre=0.55*Raylei^0.25
  }
  if(lometry==1|lometry==3|lometry==5){
    if(Raylei<1.0e-05){  
      NUfre=0.4
    }else{
      if(Raylei<0.1){   
        NUfre=0.976*Raylei^0.0784
      }else{
        if(Raylei<100){  
          NUfre=1.1173*Raylei^0.1344   
        }else{  
          if(Raylei<10000.){  
            NUfre=0.7455*Raylei^0.2167   
          }else{  
            if(Raylei<1.0e+09){  
              NUfre=0.5168*Raylei^0.2501
            }else{ 
              if(Raylei<1.0e+12){  
                NUfre=0.5168*Raylei^0.2501
              }}}}}}
  }
  
  if(lometry==2|lometry==4){
    Raylei=(GR^.25)*(PR^.333)
    NUfre=2.+0.60*Raylei
  }
  hc_free<-NUfre*THCOND/L # convection coefficent, forced
  hc_comb<-hc_free+hc_forced
  Rconv<-1/(hc_comb*ATOT)
  Nu<-hc_comb*L/THCOND # Nu combined
  hr<-4*emis*sigma*((Tc+Trad)/2+273)^3 # radiation resistance
  hc<-hc_comb
  
  if(lometry==2){
    j<-(Qabs+Qgen+hc*ATOT*((q*S2)/(2*kflesh)+Tair)+hr*ATOT*((q*S2)/(2*kflesh)+Trad))/C
  }else{
    j<-(Qabs+Qgen+hc*ATOT*((q*R^2)/(2*kflesh)+Tair)+hr*ATOT*((q*S2)/(2*kflesh)+Trad))/C
  }
  kTc<-ATOT*(Tc*hc+Tc*hr)/C
  k<-ATOT*(hc+hr)/C
  Tcf<-j/k # final Tc
  Tci<-Tc
  Tc<-(Tci-Tcf)*exp(-1*k*t)+Tcf # Tc at time t
  timethresh<-log((thresh-Tcf)/(Tci-Tcf))/(-1*k)
  tau<-(rho*V*cp)/(ATOT*(hc+hr)) # time constant
  dTc<-j-kTc
  
  return(list(Tc=Tc,Tcf=Tcf,tau=tau,dTc=dTc,timethresh=timethresh))
}
   




# constants 

#----------------------------------------- 2-12-14 new params for one lump transient
kflesh<-0.5 # thermal conductivity of flesh W/mK
q<-0

cp<-4185 #specific heat of flesh, J/kg-C
emis<-0.95 #emissivity of skin, -
sigma<-0.0000000567 #Stefan-Boltzman, W/mK
Fo_e<-0.8 #config factor, object to IR environment, -
rho<-1000 #animal density, kg/m3
abs<-0.95 #animal solar absorptivity
# 'lometry' determines whether standard or custom shapes/surface area/volume relationships are used.
# 0=plate,1=cyl,2=ellips,3=lizard (desert iguana),4=frog (leopard frog),
# 5=custom (cylinder geometry is automatically invoked when container model operates)
lometry<-3 # organism shape (see above)
# 'custallom' below operates if lometry=5, and consists of 4 pairs of values representing 
# the parameters a and b of a relationship AREA=a*mass^b, where AREA is in cm2 and mass is in g.
# The first pair are a and b for total surface area, then a and b for ventral area, then for  
# sillhouette area normal to the sun, then sillhouette area perpendicular to the sun
customallom<-c(10.4713,.688,0.425,0.85,3.798,.683,0.694,.743) # custom allometry coefficients (see above)
shape_a<-1. 
shape_b<-3.16666666667
shape_c<-0.6666666667
posture<-'n' # pointing normal 'n' or parallel 'p' to the sun's rays?
FATOSK<-0.4 # configuration factor to sky
FATOSB<-0.4 # configuration factor to substrate

press<-101325 #atmospheric pressure, pa
sub_reflect<-0.2 # solar reflectance of substrate
pctdif<-0.1 # proportion of solar energy that is diffuse (rather than direct beam)
# RNL_new trans model_with DEB_1.6
# RNL_new trans model_with DEB_1.6

# RNL_new trans model_with DEB_1.6

tzone<-paste("Etc/GMT-",10,sep="")
 
 
 


metout<-read.csv('/Applications/Programs/NetLogo 5.0.5/NicheMapR/Transient heat budget model/metout.csv')
soil<-read.csv('/Applications/Programs/NetLogo 5.0.5/NicheMapR/Transient heat budget model/soil.csv')
shadmet<-read.csv('/Applications/Programs/NetLogo 5.0.5/NicheMapR/Transient heat budget model/shadmet.csv')
shadsoil<-read.csv('/Applications/Programs/NetLogo 5.0.5/NicheMapR/Transient heat budget model/shadsoil.csv')
micro_sun_all<-cbind(metout[,2:5],metout[,9],metout[,11],metout[,14:16])
colnames(micro_sun_all)<-c('dates','JULDAY','TIME','TALOC','VLOC','TS','ZEN','SOLR','TSKYC')
micro_shd_all<-cbind(metout[,2],shadmet[,2:4],shadmet[,8],shadmet[,10],shadmet[,13:15])
colnames(micro_shd_all)<-c('dates','JULDAY','TIME','TALOC','VLOC','TS','ZEN','SOLR','TSKYC')


# choose a day(s) to simulate

daystart<-paste('09/11/10',sep="") # yy/mm/dd
dayfin<-paste('10/12/31',sep="") # yy/mm/dd

micro_sun<-subset(micro_sun_all, format(as.POSIXlt(micro_sun_all$dates), "%y/%m/%d")>=daystart & format(as.POSIXlt(micro_sun_all$dates), "%y/%m/%d")<=dayfin)
micro_shd<-subset(micro_shd_all, format(as.POSIXlt(micro_shd_all$dates), "%y/%m/%d")>=daystart & format(as.POSIXlt(micro_shd_all$dates), "%y/%m/%d")<=dayfin)
days<-as.numeric(as.POSIXlt(dayfin)-as.POSIXlt(daystart))

time<-seq(0,(days+1)*60*24,60) #60 minute intervals from microclimate output
time<-time[-1]
times2<-seq(0,(days+1)*60*24,2) #two minute intervals for prediction
time<-time*60 # minutes to seconds
times2<-times2*60 # minutes to seconds


  Qsolf_sun<- approxfun(time, micro_sun[,8], rule = 2)
  Tradf_sun<- approxfun(time, rowMeans(cbind(micro_sun[,6],micro_sun[,9])), rule = 2)
  velf_sun<- approxfun(time, micro_sun[,5], rule = 2)
  Tairf_sun<- approxfun(time, micro_sun[,4], rule = 2)
  Zenf_sun<- approxfun(time, micro_sun[,7], rule = 2)

  Qsolf_shd<- approxfun(time, micro_shd[,8]*.1, rule = 2)
  Tradf_shd<- approxfun(time, rowMeans(cbind(micro_shd[,6],micro_shd[,9])), rule = 2)
  velf_shd<- approxfun(time, micro_shd[,5], rule = 2)
  Tairf_shd<- approxfun(time, micro_shd[,4], rule = 2)
  Zenf_shd<- approxfun(time, micro_shd[,7], rule = 2)

Tc_min<- 3.7 # from Bundey field site (2-9-14)
Tc_max<- 51 # from max(results$Tb)


# ***************** end TRANSIENT MODEL SETUP ***************
#************************************************************







#************************************************************
# *********************** DEB MODEL *************************


# physiological traits
TMAXPR<-36 #34 ^ degrees C, voluntary thermal maximum (upper body temperature for foraging) Pamula 1997 - where frequency dropped substantially, rather than extreme (Fig. 3.42)
TMINPR<-26.0 #26.0 # ^ degrees C, voluntary thermal minimum (lower body temperature for foraging) Pamula 1997 (Fig. 3.42)
TBASK<-19#26.#23.1 # degrees C, minimum basking temperature Pamula Table 3.14
TEMERGE<-8.5#8.5 # degrees C, temperature at which animal will move to a basking site *based on Kerr and Bull 2004
ctmax<-43.0 # ^ degrees C, critical thermal maximum (used by program to determine depth selected when inactive and burrowing) (43.0, Bennett, A.F. & John-Alder, H. (1986) Thermal Relations of Some Australian Skinks (Sauria: Scincidae). Copeia, 1986, 57-64.)
ctmin<-3.5 # ^ degrees C, critical thermal minimum (used by program to determine depth selected when inactive and burrowing) (3.5, Bennett, A.F. & John-Alder, H. (1986) Thermal Relations of Some Australian Skinks (Sauria: Scincidae). Copeia, 1986, 57-64.)
ctminthresh<-12 #number of consecutive hours below CTmin that leads to death
ctkill<-1 #if 1, animal dies when it hits critical thermal limits
TPREF<-33.5 # ^ preferred body temperature (animal will attempt to regulate as close to this value as possible) (mean 31.9, range 29.4-34.3, Bennett, A.F. & John-Alder, H. (1986) Thermal Relations of Some Australian Skinks (Sauria: Scincidae). Copeia, 1986, 57-64.), mode in Pamula Fig. 3.14 around 33.5
DELTAR<-0.1 # degrees C, temperature difference between expired and inspired air
skinwet<-0.5 # estimated from data in Bently 1959 at 23 degrees and 34.5 degrees #0.2#0.35 # %, of surface area acting like a free water surface (e.g. most frogs are 100% wet, many lizards less than 5% wet)
extref<-20. # %, oxygen extraction efficiency (need to check, but based on 35 deg C for a number of reptiles, from Perry, S.F., 1992. Gas exchange strategies in reptiles and the origin of the avian lung. In: Wood, S.C., Weber, R.E., Hargens, A.R., Millard, R.W. (Eds.), Physiological Adaptations in Vertebrates: Respiration, Circulation, andMetabo -  lism. Marcel Dekker, Inc., New York, pp. 149-167.)
PFEWAT<-73. # %, fecal water (from Shine's thesis, mixed diet 75% clover, 25% mealworms)
PTUREA<-0. # %, water in excreted nitrogenous waste
FoodWater<-82#82 # 82%, water content of food (from Shine's thesis, clover)
minwater<-9.5 # %, minimum tolerated dehydration (% of wet mass) - prohibits foraging if greater than this
raindrink<-5. # daily rainfall (mm) required for animal to rehydrate from drinking (zero means standing water always available)
gutfill<-75. # % gut fill at which satiation occurs - if greater than 100%, animal always tries to forage

# behavioural traits
dayact<-1 # diurnal activity allowed (1) or not (0)?
nocturn<-0 # nocturnal activity allowed (1) or not (0)?
crepus<-0 # crepuscular activity allowed (1) or not (0)?
burrow<-1 # shelter in burrow allowed (1) or not (0)?
shdburrow<-1 #
mindepth<-2 # minimum depth (soil node) to which animal can retreat if burrowing
maxdepth<-10 # maximum depth (soil node) to which animal can retreat if burrowing
CkGrShad<-1 # shade seeking allowed (1) or not (0)?
climb<-0 # climbing to seek cooler habitats allowed (1) or not (0)?
fosorial<-0 # fossorial activity (1) or not (0)
rainact<-0 # activity is limited by rainfall (1) or not (0)?
actrainthresh<-0.1 # threshold mm of rain causing activity if rainact=1
breedactthresh<-1 # threshold numbers of hours active after start of breeding season before eggs can be laid (simulating movement to the breeding site)


fract<-1
f<-1.
MsM<-186.03*6. # produces a stomach volume of 5.3 cm3/100 g, as measured for Disosaurus dorsalis
z<-7.9286*fract
delta<-0.2401
kappa_X<-0.85#0.85
v_dotref<-0.063289/24.
kappa<-0.84416 
p_Mref<-31.804/24.
E_G<-7772
k_R<-0.95
k_J<-0.002/24.
E_Hb<-102980
E_Hj<-E_Hb*fract^3
E_Hp<-249580
h_aref<-2.8503e-10/(24.^2) 
s_G<-0.01

E_Egg<-1068400*fract^4# J, initial energy of one egg # this includes the residual yolk, which is eaten upon hatching
svl_met<-11 # mm, snout vent length at metamorphosis
E_m<-(p_Mref*z/kappa)/v_dotref
p_Xm<-12420/24 # J/h.cm2, maximum intake rate when feeding
K<-10 # half-saturation constant
X<-3265 # food density J/cm2

# for insect model
metab_mode<-0 # 0 = off, 1 = holometabolous with Dyar's rule scaling, 2 = holometabolous linear scaling, 3 = hemimetabolous with Dyar's rule scaling, 4 = hemimetabolous linear scaling
stages<-8 # number of stages (max = 8) = number of instars plus 1 for egg + 1 for pupa + 1 for imago
p_Am1<-0.9296852/24*100
p_AmIm<-2.068836/24*100
disc<-0.0307
gam<-1.6

# these next five parameters control the thermal response, effectively generating a thermal response curve
T_REF<-20 # degrees C, reference temperature - correction factor is 1 for this temperature
TA<-8819.8
TAL<-18265
TAH<-20844
TL<-288.
TH<-315.

arrhenius<-matrix(data = 0, nrow = 8, ncol = 5)
arrhenius[,1]<-TA # critical thermal minimum
arrhenius[,2]<-TAL # critical thermal maximum
arrhenius[,3]<-TAH # voluntary thermal minimum
arrhenius[,4]<-TL # voluntary thermal maximum
arrhenius[,5]<-TH # basking threshold 

thermal_stages<-matrix(data = 0, nrow = 8, ncol = 6)
thermal_stages[,1]<-ctmin # critical thermal minimum
thermal_stages[,2]<-ctmax # critical thermal maximum
thermal_stages[,3]<-TMINPR # voluntary thermal minimum
thermal_stages[,4]<-TMAXPR # voluntary thermal maximum
thermal_stages[,5]<-TBASK # basking threshold
thermal_stages[,6]<-TPREF # preferred body temperature

behav_stages<-matrix(data = 0, nrow = 8, ncol = 14)

behav_stages[,1]<-dayact
behav_stages[,2]<-nocturn
behav_stages[,3]<-crepus
behav_stages[,4]<-burrow
behav_stages[,5]<-shdburrow
behav_stages[,6]<-mindepth
behav_stages[,7]<-maxdepth
behav_stages[,8]<-CkGrShad
behav_stages[,9]<-climb
behav_stages[,10]<-fosorial
behav_stages[,11]<-rainact
behav_stages[,12]<-actrainthresh
behav_stages[,13]<-breedactthresh
behav_stages[,14]<-flyer

water_stages<-matrix(data = 0, nrow = 8, ncol = 8)

water_stages[,1]<-skinwet
water_stages[,2]<-extref
water_stages[,3]<-PFEWAT
water_stages[,4]<-PTUREA
water_stages[,5]<-FoodWater
water_stages[,6]<-minwater
water_stages[,7]<-raindrink
water_stages[,8]<-gutfill

# composition related parameters
andens_deb<-1. # g/cm3, density of structure 
d_V<-0.3 # density of structure (reflects fraction of mass that is dry)
d_E<-0.3 # density of reserve (reflects fraction of mass that is dry)
eggdryfrac<-0.3 # decimal percent, dry mass of eggs
mu_X<-525000 # J/cmol, chemical potential of food
mu_E<-585000 # J/cmol, chemical potential of reserve
mu_V<-500000 # J/cmol, chemical potential of structure 
mu_P<-480000 # J/cmol, chemical potential of product (faeces)
kappa_X_P<-0.1 # fraction of food energy into faeces

          # elemental maxtrix of organics  
nX<-c(1,1.8,0.5,.15) # composition of food (atoms per carbon atoms for CHON)
nE<-c(1,1.8,0.5,.15) # composition of reserve (atoms per carbon atoms for CHON)
nV<-c(1,1.8,0.5,.15) # composition of structure (atoms per carbon atoms for CHON)
nP<-c(1,1.8,0.5,.15) # composition of product/faeces (atoms per carbon atoms for CHON)
N_waste<-c(5,4,3,4) # chemical formula for nitrogenous waste product, CHON, e.g. Urea c(0,3,0,1), Uric acid c(5,4,3,4)


# breeding life history
clutchsize<-2. # clutch size
eggmass<-3.787 # initial dry mass of an egg (g)
viviparous<-1 # 1=yes, 0=no
batch<-1 # invoke Pequerie et al.'s batch laying model?

# the following four parameters apply if batch = 1, i.e. animal mobilizes
breedrainthresh<-0 # rain dependent breeder? 0 means no, otherwise enter rainfall threshold in mm
# photoperiod response triggering ovulation, none (0), summer solstice (1), autumnal equinox (2),  
# winter solstice (3), vernal equinox (4), specified daylength thresholds (5)
photostart<- 5 # photoperiod initiating breeding
photofinish<- 5 # photoperiod terminating breeding
daylengthstart<- 12.5 # threshold daylength for initiating breeding
daylengthfinish<- 13.8 # threshold daylength for terminating breeding
photodirs <- 1 # is the start daylength trigger during a decrease (0) or increase (1) in day length?
photodirf <- 1 # is the finish daylength trigger during a decrease (0) or increase (1) in day length?
startday<-1 # make it 90 for T. rugosa loop day of year at which DEB model starts
breedtempthresh<-200 # body temperature threshold below which breeding will occur
breedtempcum<-24*7 # cumulative time below temperature threshold for breeding that will trigger breeding

reset<-0 # reset options, 0=quit simulation upon death, 1=restart at emergence, 2=restart at first egg laid, 3=restart at end of breeding season, 4=reset at death

# frog breeding mode 0 is off, 
# 1 is exotrophic aquatic (eggs start when water present in container and within breeding season)
# 2 is exotrophic terrestrial/aquatic (eggs start at specified soil node within breeding season, 
# diapause at birth threshold, start larval phase if water present in container)
# 3 endotrophic terrestrial (eggs start at specified soil node within breeding season and continue
# to metamorphosis on land)
# 4 turtle mode (eggs start at specified soil node within breeding season, hatch and animals enter
# water and stay there for the rest of their life, but leave the water if no water is present)
frogbreed<-0 # frog breeding mode
frogstage<-0 # 0 is whole life cycle, 1 is just to metamorphosis ({ reset and start again)

# metabolic depression
aestivate<-0
depress<-0.2



#*********************************** DEB model initial conditions **************************************


v_init<-3e-9
E_init<-E_Egg/v_init
E_H_init<-0
stage<-0
v_init<-(3.9358^3)*fract^3 #hatchling
E_init<-E_m
E_H_init<-E_Hb+5
stage<-1
v_init<-(498.4132^3)*fract^3*0.85 # adult 
E_init<-E_m
E_H_init<-E_Hp+1
stage<-3
ma<-1e-4  # hourly active mortality rate (probability of mortality per hour)
mi<-0  # hourly inactive mortality rate (probability of mortality per hour)
mh<-0.5   # survivorship of hatchling in first year
# DEB model initial conditions
V_init_baby<-3e-9
E_init_baby<-E_Egg/V_init_baby
E_baby_init<-E_init_baby
V_baby_init<-V_init_baby
ms_init<-0.
cumrepro_init<-0.
q_init<-0.
hs_init<-0.
cumbatch_init<-0.
pregnant<-0
E_m<-(p_Mref*z/kappa)/v_dotref

# conversions from percent to proportion
PTUREA1<-PTUREA/100
PFEWAT1<-PFEWAT/100
FoodWater1<-FoodWater/100
water_stages[,3]<-water_stages[,3]/100
water_stages[,4]<-water_stages[,4]/100
water_stages[,5]<-water_stages[,5]/100
eggmass<-0 # initial dry mass of an egg (g) - no longer used so delete



#******************** DEB mass balance calculations ************************

nO<-cbind(nX,nV,nE,nP) # matrix of composition of organics, i.e. food, structure, reserve and faeces
CHON<-c(12,1,16,14)
wO<-CHON%*%nO
w_V=wO[3]
M_V<-d_V/w_V
yEX<-kappa_X*mu_X/mu_E # yield of reserve on food
yXE<-1/yEX # yield of food on reserve
yVE<-mu_E*M_V/E_G  # yield of structure on reserve
yPX<-kappa_X_P*mu_X/mu_P # yield of faeces on food
yXP<-1/yPX # yield of food on faeces
yPE<-yPX/yEX # yield of faeces on reserve  0.143382353
nM<-matrix(c(1,0,2,0,0,2,1,0,0,0,2,0,N_waste),nrow=4)
N_waste_inv<-c(-1*N_waste[1]/N_waste[4],(-1*N_waste[2])/(2*N_waste[4]),(4*N_waste[1]+N_waste[2]-2*N_waste[3])/(4*N_waste[4]),1/N_waste[4])
nM_inv<-matrix(c(1,0,-1,0,0,1/2,-1/4,0,0,0,1/2,0,N_waste_inv),nrow=4)
JM_JO<--1*nM_inv%*%nO
etaO<-matrix(c(yXE/mu_E*-1,0,1/mu_E,yPE/mu_E,0,0,-1/mu_E,0,0,yVE/mu_E,-1/mu_E,0),nrow=4)
w_N<-CHON%*%N_waste


############################# end input #########################################

# fixes Mike made to get it to work in R, including case-sensitive issues and other symbol changes

w_X=wO[1]
w_E=wO[3]
w_V=wO[2]
w_P=wO[4]

T_A<-TA
T_ref<-T_REF
E_egg<-E_Egg
ANDENS_deb<-andens_deb
k_Jref<-k_J
zfact<-z
vdotref<-v_dotref
p_Xmref<-p_Xm
waiting<-0
hour<-1
daycount<-1
lambda=6./12.
breeding<-1
surviv_init<-1
delta_deb<-delta
halfsat<-K


funct<-f
fecundity<-0
clutches<-0

cumrepro_prev<-0
cumbatch_prev<-0
cumbatch_init<-0
cumrepro_init<-0


stage<-3
dead<-0

E_pres<-E_init
V_pres<-v_init
E_H_pres<-E_H_init
q_pres<-q_init
hs_pres<-hs_init
surviv_pres<-surviv_init
ms_pres<-ms_init


month<-"nov1"
Tairf<-Tairf_shd
Tc_init<-Tairf(1)+0.1 # Initial core temperature
mass<-800             # Weight in grams
NL_shade<-100          # Shade patches 
NL_food<-10          # Food patches 
NL_T_b<-Tc_init       # Initial T_b 
NL_T_b_min<-Tc_min         # Min critical T_b
NL_T_b_max<-Tc_max        # Max critical T_b
NL_T_opt_l<-28
NL_T_opt_u<-35
NL_ctminthresh<-ctminthresh # No. of consecutive hours below CTmin that leads to death
NL_reserve<-E_m        # Initial reserve density
NL_max_reserve<-E_m    # Maximum reserve level
NL_maint<-0               # Maintenance cost
NL_move<-0                # Movement cost
NL_zen<-Zenf(1*60*60)     # Zenith angle 


# -------------- extracting waddle data ---------------------
install.packages(c("raster","zoo","stringr","rgdal"))
library(raster)
library(zoo)
library(stringr)
library(rgdal)
# pull waddle data 
lizard<-read.csv('/Users/camel/Desktop/Matt2016/Data/waddleometer/11885_2009_ALL.csv',stringsAsFactors=FALSE)
lizard$Month<-str_trim(lizard$Month)
lizard$Day<-str_trim(lizard$Day)
lizard<-as.data.frame(lizard)
#change days and months to two chars
lizard$Day[nchar(lizard$Day)==1]<-paste(0,lizard$Day[nchar(lizard$Day)==1],sep="")
lizard$Month[nchar(lizard$Month)==1]<-paste(0,lizard$Month[nchar(lizard$Month)==1],sep="")

surf<-raster('/Users/camel/Desktop/Matt2016/Data/Church_DSM.tif')
# load dtm
terr<-raster('/Users/camel/Desktop/Matt2016/Data/Church_DTM_corrected_projection.tif')
# calc the difference between the two layers to get vegetation heights
difr<-surf-terr 
# write this to a raster file
setwd("/Users/camel/Desktop/Matt2016/Data/")
writeRaster(difr,filename='veg.tif')

X1<-343000
X2<-344500
Y1<-6248400
Y2<-6249715 # from max(difr)

x<-cbind(lizard$Easting,lizard$Northing); x<-as.data.frame(x)
colnames(x)<-c("X","Y")
# extract from the raster the above lizard's easting and northing points. DTM and DSM files must contain these GPS points
shade<-extract(difr,x) 
# add this shade column to the waddle data
lizard<-cbind(lizard,shade) 
lizard[,3:10]<-na.locf(lizard[,3:10],na.rm=FALSE)
lizard[,17]<-na.locf(lizard[,17],na.rm=FALSE)
lizard<-na.omit(lizard)
#lizard<-subset(lizard,Easting<600000)
#lizard<-subset(lizard,Northing<6249000)
#lizard<-subset(lizard,Northing<6250000)
#lizard<-subset(lizard,Steps>5)
#lizard<-subset(lizard,Hours>12 | Hours<10)

# remove GPS outliers   ...... ???
lizard$Easting<-as.numeric(lizard$Easting); lizard$Northing<-as.numeric(lizard$Northing)
lizard<-subset(lizard,subset=lizard$Easting <  344765 & lizard$Easting > 342747) #from max(difr)
lizard<-subset(lizard,subset=lizard$Northing < 6249715 & lizard$Northing > 6247878) #from max(difr)
max(lizard$Northing); max(lizard$Easting)


#plot(difr,xlim=c(min(x[,1],na.rm=TRUE),max(x[,1],na.rm=TRUE)),ylim=c(min(x[,2],na.rm=TRUE),max(x[,2],na.rm=TRUE)),zlim=c(0.5,5)) # plot entire site
#plot(difr,xlim=c(343450,343500),ylim=c(6248800,6248850),zlim=c(-30,0.05)) # plot a subregion
plot(difr,xlim=c(X1,X2),ylim=c(Y1,Y2),zlim=c(0.5,5)) # plot a subregion
colvec = adjustcolor(lizard$Month, alpha = 0.1)
with(lizard,points(lizard$Northing~lizard$Easting,cex=0.5,pch=20,col=colvec)) # plot sleepy GPS


library(rgdal) 
library(adehabitatHR)
# home range. Use percent = 95 to eliminate outliers
lizard$Easting<-as.numeric(lizard$Easting); lizard$Northing<-as.numeric(lizard$Northing)
lizard$Month<-as.numeric(lizard$Month)
liz<-cbind(lizard$Easting,lizard$Northing)
lizm<-lizard$Month

spdf<-SpatialPointsDataFrame(c(lizard$Easting,lizard$Northing), lizard$Month)# creates a spatial points data frame (adehabitatHR package)
homerange<-mcp(spdf,percent=100)

# plot home range polygon of real lizard
liz<-cbind(lizard["Easting"],lizard["Northing"])
df<-data.frame(rep.int(1,length(liz[[1]])))
length(df)
spdf<-SpatialPointsDataFrame(liz[1:2],data=df)
homerange<-mcp(spdf,100)
colvec = adjustcolor(c("light blue"), alpha = 0.5)
plot(homerange,col=colvec,border=colvec,
     xlim=c(X1,X2),
     ylim=c(Y1,Y2),
     add=T, 
     #axes=T,
     )



kud <- kernelUD(spdf, grid=100, extent=0.2, h="href", same4all=T) # class 'estUD'
image(kud)
ver <- getverticeshr(kud,95); class(ver)
plot(ver, col=adjustcolor(ver$id,alpha=0.6),border=ver$id,xlab="Easting",ylab="Northing",
     main="getverticeshr() function of kernelUD()", 
     #col=rainbow(4), 
) # plot kernelUD 


# change this with each sim
hr7<-homerange
hrpath7<-spdf 

# original code. NB: using min/maxpcors doesn't work for exporting plot to pdf or jpeg
minpxcor<-NLReport("min-pxcor");maxpxcor<-NLReport("max-pxcor")
minpycor<-NLReport("min-pycor");maxpycor<-NLReport("max-pycor")
# plot 75% and 100% gutfull HRs
colvec = adjustcolor(c("black"), alpha = 0.5)
plot(hr7,lty=1,bty="o",pch=21,col=colvec, xlim=c(minpxcor,maxpxcor),ylim=c(minpycor,maxpycor),axes=F)
plot(hrpath7,pch=3,col=turtles$days,add=T)

colvec = adjustcolor(c("red"), alpha = 0.5)
plot(hr7,lty=1,bty="o",col=colvec,add=T)
plot(hrpath7,pch=3,col=turtles$days,add=T) 

# **********************************************************
# ******************** start sleepy sim  ***********

# change parenthesised number to number of days to simulate  
#dayf<-tail(seq(0,(2)*60*24,60),1); dayf
#define start and finish dates
#daystart<-paste(substr(lizard[1,11],3,4),lizard[1,12],lizard[1,13],sep="/") # yy/mm/dd
#dayfin<-paste(substr(lizard[2880,11],3,4),lizard[2880,12],lizard[2880,13],sep="/") # yy/mm/dd

# choose a day(s) to simulate
daystart<-paste('09/09/05',sep="") # yy/mm/dd
dayfin<-paste('10/12/31',sep="") # yy/mm/dd

install.packages("NicheMapR")
library(NicheMapR)
source('/Users/camel/Desktop/Matt2016/NicheMapR/onelump/DEB.R')
source('/Users/camel/Desktop/Matt2016/NicheMapR/onelump/onelump_varenv.R')

# read in microclimate data
tzone<-paste("Etc/GMT-",10,sep="")
metout<-read.csv('/Users/camel/Desktop/Matt2016/NicheMapR/microclim/metout.csv')
soil<-read.csv('/Users/camel/Desktop/Matt2016/NicheMapR/microclim/soil.csv')
shadmet<-read.csv('/Users/camel/Desktop/Matt2016/NicheMapR/microclim/shadmet.csv')
shadsoil<-read.csv('/Users/camel/Desktop/Matt2016/NicheMapR/microclim/shadsoil.csv')
micro_sun_all<-cbind(metout[,2:5],metout[,9],metout[,11],metout[,14:16])
colnames(micro_sun_all)<-c('dates','JULDAY','TIME','TALOC','VLOC','TS','ZEN','SOLR','TSKYC')
micro_shd_all<-cbind(metout[,2],shadmet[,2:4],shadmet[,8],shadmet[,10],shadmet[,13:15])
colnames(micro_shd_all)<-c('dates','JULDAY','TIME','TALOC','VLOC','TS','ZEN','SOLR','TSKYC')


micro_sun<-subset(micro_sun_all, format(as.POSIXlt(micro_sun_all$dates), "%y/%m/%d")>=daystart & format(as.POSIXlt(micro_sun_all$dates), "%y/%m/%d")<=dayfin)
micro_shd<-subset(micro_shd_all, format(as.POSIXlt(micro_shd_all$dates), "%y/%m/%d")>=daystart & format(as.POSIXlt(micro_shd_all$dates), "%y/%m/%d")<=dayfin)
days<-as.numeric(as.POSIXlt(dayfin)-as.POSIXlt(daystart));days

time<-seq(0,(days+1)*60*24,60) #60 minute intervals from microclimate output
time<-time[-1]
times2<-seq(0,(days+1)*60*24,2) #two minute intervals for prediction
time<-time*60 # minutes to seconds
times2<-times2*60 # minutes to seconds


  Qsolf_sun<- approxfun(time, micro_sun[,8], rule = 2)
  Tradf_sun<- approxfun(time, rowMeans(cbind(micro_sun[,6],micro_sun[,9])), rule = 2)
  velf_sun<- approxfun(time, micro_sun[,5], rule = 2)
  Tairf_sun<- approxfun(time, micro_sun[,4], rule = 2)
  Zenf_sun<- approxfun(time, micro_sun[,7], rule = 2)

  Qsolf_shd<- approxfun(time, micro_shd[,8]*.1, rule = 2)
  Tradf_shd<- approxfun(time, rowMeans(cbind(micro_shd[,6],micro_shd[,9])), rule = 2)
  velf_shd<- approxfun(time, micro_shd[,5], rule = 2)
  Tairf_shd<- approxfun(time, micro_shd[,4], rule = 2)
  Zenf_shd<- approxfun(time, micro_shd[,7], rule = 2)

#sleepy_ticks<-days / (2 / 60 / 24)  No. of NL ticks (measurement of days)
#for (i in Sleepy_ticks) { }

# choose sun or shade
tick<-i
times3<-c(times2[tick],times2[tick+1])

# ----------------------------------- 2-12-14 new one_lump_trans params 
Qsol<-Qsolf(mean(times3))
vel<-velf(mean(times3)) ;vel
Tair<-Tairf(mean(times3));Tair
Trad<-Tradf(mean(times3)); Trad
Zen<-Zenf(mean(times3)); Zen

#indata<-list(cp=cp,emis=emis,sigma=sigma,Fo_e=Fo_e,rho=rho,abs=abs,lometry=lometry,customallom=customallom,shape_a=shape_a,shape_b=shape_b,shape_c=shape_c,posture=posture,FATOSK=FATOSK,FATOSB=FATOSB,mass=mass,sub_reflect=sub_reflect,pctdif=pctdif)
input<-list(kflesh=kflesh,q=q,cp=cp,emis=emis,Fo_e=Fo_e,rho=rho,abs=abs,lometry=lometry,customallom=customallom,shape_a=shape_a,shape_b=shape_b,shape_c=shape_c,posture=posture,FATOSK=FATOSK,FATOSB=FATOSB,mass=mass,sub_reflect=sub_reflect,pctdif=pctdif,Qsol=Qsol,vel=vel,Tair=Tair,Trad=Trad,Zen=Zen)

#lizTb<-lizard$Temperature

#new trans function to update lizard$Temp
trans2 <- function(temp,input)  {
	
# transient is function(t,y,thresh,input)
# calc Tb params at 2 mins interval 
	Tbs<-transient(120,temp,30,input) # (old ode model) as.data.frame(ode(y=Tc_init,times=times3,func=transient,parms=indata))

}

	Tb<-Tbs$Tc 
	rate<-Tbs$dTc
	Tc_init<-Tb
	thresh<-c(25) # threshold body temperature at which time is required (deg C)


Tempnew <-0
Tempnew[1]<-lizard[1,1]
for (j in 2:length(lizard$shade)){
    if(lizard$shade[j] <= 0){
      Qsolf<-Qsolf_sun
      Tradf<-Tradf_sun
      velf<-velf_sun  
      Tairf<-Tairf_sun
      Zenf<-Zenf_sun
    }else{
      Qsolf<-Qsolf_shd 
      Tradf<-Tradf_shd
      velf<-velf_shd   
      Tairf<-Tairf_shd
      Zenf<-Zenf_shd
    }
    offset<-length(seq(0,(lizard$Hours[1]*60*60+lizard$Minutes[1]*60),120))
      Qsol<-Qsolf(times2[j+offset])
      Trad<-Tradf(times2[j+offset])
      vel<-velf(times2[j+offset]) 
      Tair<-Tairf(times2[j+offset])
      Zen<-Zenf(times2[j+offset])
	input<-list(kflesh=kflesh,q=q,cp=cp,emis=emis,Fo_e=Fo_e,rho=rho,abs=abs,lometry=lometry,customallom=customallom,shape_a=shape_a,shape_b=shape_b,shape_c=shape_c,posture=posture,FATOSK=FATOSK,FATOSB=FATOSB,mass=mass,sub_reflect=sub_reflect,pctdif=pctdif,Qsol=Qsol,vel=vel,Tair=Tair,Trad=Trad,Zen=Zen)
	  Tempnew[j] <- as.numeric(trans2(Tempnew[j - 1]),input)
#		  lizard$Tempnew[j] <- trans2(Tbs$Tc)  
  }

lizard$Tempnew<-Tempnew



### ------------- plotting ------------------

     
# ------------------------ export plots from NL --------------------------

month<-"sep"

#dir = /Applications/Programs/NetLogo 5.0.5/Soft foraging model/Simulations

#spatial plot
sfh<-paste(month,NL_days,round(mass,0),NL_shade,NL_food*10,"_move","",sep="");sfh
#NLCommand("export-plot \"Spatial coordinates of transition between activity states\" \"Simulations/spatialplot.csv\"")
NLCommand(paste("export-plot \"Spatial coordinates of transition between activity states\" \"Simulations/",sfh,".csv\"",sep=""))
# home range
hfh<-paste(month,NL_days,round(mass,0),NL_shade,NL_food*10,"_homerange","",sep="");hfh
# reserve plot
rfh<-paste(month,NL_days,round(mass,0),NL_shade,NL_food*10,"_reserve","",sep="");rfh
#NLCommand("export-plot \"Reserve level and starvation reserve over time\" \"Simulations/reserveplot.csv\"")
NLCommand(paste("export-plot \"Reserve level and starvation reserve over time\" \"Simulations/",rfh,".csv\"",sep=""))
#temp plot
tfh<-paste(month,NL_days,round(mass,0),NL_shade,NL_food*10,"_temp",sep="");tfh
#NLCommand("export-plot \"Body temperature (T_b)\" \"Simulations/tempplot.csv\"")
NLCommand(paste("export-plot \"Body temperature (T_b)\" \"Simulations/",tfh,".csv\"",sep=""))
# activity budget
afh<-paste(month,NL_days,round(mass,0),NL_shade,NL_food*10,"_act","",sep="");afh
#NLCommand("export-plot \"Global time budget\" \"Simulations/activitybudget.csv\"")
NLCommand(paste("export-plot \"Global time budget\" \"Simulations/",afh,".csv\"",sep=""))
# world view
NLCommand("export-view \"Simulations/worldview.png\"")
# text output
xfh<-paste(month,NL_days,round(mass,0),NL_shade,NL_food*10,"_txt",sep="");xfh
NLCommand(paste("export-output \"Simulations/",xfh,".csv\"",sep=""))
# gut level
gfh<-paste(month,NL_days,round(mass,0),NL_shade,NL_food*10,"_gut","",sep="");gfh
#NLCommand("export-plot \"Global time budget\" \"Simulations/activitybudget.csv\"")
NLCommand(paste("export-plot \"Gutfull\" \"Simulations/",gfh,".csv\"",sep=""))




#-----------------------------------------------------------
#------------------------ plot NETLOGO RESULTS -------------
#-----------------------------------------------------------

###
### add export function that saves plot as e.g. "sep15"
###

#fh<-paste(daystart,NL_days,mass,NL_shade,NL_food,sep="")
ttl<-paste("From ",daystart," + ",NL_days,"days ;","Weight =",mass,"g ; Shade =", NL_shade,"; Food =", NL_food,"; VTMIN", NL_T_b_min,"C")
pdf(paste("/Applications/Programs/NetLogo 5.0.5/NL_transient_outputs/",fh,".pdf",sep=""),width=15,height=15,paper="a4r",title=ttl)

par(mfrow=c(1,1)) # new plot window with x rows and y columns, fills by rows
par(mar=c(5,6,5,5))

#------------------ regular results plot ------------------------

#NL_shade<-NLReport("Shade-patches"); NL_food<-NLReport("Food-patches")
results<-as.data.frame(results)
ticktime<-results$tick * 2 / 60 / 24 # convert to real days

# Plot with title
# with(results,plot(Tb~ticktime,type='l',las=1,xlab="Days",ylab = expression(paste("Body temperature (" * degree,"C)")),main=paste("From ",daystart," + ",NL_days,"days ;","Weight =",mass,"g ; Shade =", NL_shade,"; Food =", NL_food,"; VTMIN", NL_T_b_min,"C")))
with(results,plot(Tb~ticktime,type='l',las=1,xlab="Days",ylab = expression(paste("Body temperature (" * degree,"C)"))))
shade.results<-results[results$shade %in% 1,] ; head(shade.results)
sticktime<-shade.results$tick * 2 / 60 / 24 # convert to real days
shade.results_Tb<-results[results$shade %in% 1,'Tb'] ; head(shade.results_Tb)
#with(results,points(shade.results$Tb~shade.results$tick, type='p', col="blue"))
with(results,points(shade.results$Tb~sticktime, type='p', col="blue"))
#with(results,points(shade*20~tick,type='p',col="red"))
abline(h = c(NL_T_opt_l,NL_T_opt_u), col = "red", lty = 3)
text(0,28, "Activity range", col = "red", adj = c(.3, 1))

# dev.off()  # use only with jpeg function above

hist(shade.results$Tb, main="Proportion of Tb when in shade", xlab="Tb (C) when in shade")
plot(shade.results$tick,shade.results$Tb, col="blue")

# ------------------ Tb plot per month ------------------------


if (exists("results")){
  newones<-results
}
# change this to reflect new data period and save as dataframe
dec15<-newones

ticktime<-newones$tick * 2 / 60 / 24 # convert to real days
# define shade data frame
shade.newones<-newones[newones$shade %in% 1,]
shade.newones_Tb<-newones[newones$shade %in% 1,'Tb'] ; head(shade.newones_Tb)
sticktime<-shade.newones$tick * 2 / 60 / 24 # convert to real days

#sep15: lty=3
#nov1: lty=2
#dec15: lty=1

month<-"dec15"
par(pty="m")
fh<-paste(month,NL_days,mass,NL_shade,NL_food*10,"_2",sep="")
ttl<-paste("From ",daystart," + ",NL_days,"days ;","Weight =",mass,"g ; Shade =", NL_shade,"; Food =", NL_food,"; VTMIN", NL_T_b_min,"C")
pdf(paste("/Users/matthewmalishev/Documents/Manuscripts/Malishev and Kearney/Figures/Simulations/Tb plot/",fh,".pdf",sep=""),width=15,height=15,paper="a4r",title=ttl)
#plot.new()

#par(new=T)

# make data points transparent
colvec = adjustcolor(c("red"), alpha = 0.5)
col=colvec[sep15$Tb]

# Plot with title
# with(results,plot(Tb~ticktime,type='l',las=1,xlab="Days",ylab = expression(paste("Body temperature (" * degree,"C)")),main=paste("From ",daystart," + ",NL_days,"days ;","Weight =",mass,"g ; Shade =", NL_shade,"; Food =", NL_food,"; VTMIN", NL_T_b_min,"C")))
with(newones,plot(Tb~ticktime,type='l',las=1,lwd=1,lty=1,col=colvec,xlim=c(0,NL_days),ylim=c(0,45),xlab="Days",ylab = expression(paste("Body temperature (" * degree,"C)"))))
#with(results,points(shade.results$Tb~shade.results$tick, type='p', col="blue"))
with(newones,points(shade.newones$Tb~sticktime, type='p', col=colvec))
#with(results,points(shade*20~tick,type='p',col="red"))
abline(h = c(NL_T_opt_l,NL_T_opt_u), col = "red", lty = 3)
text(0,28, "Activity range", col = "red", adj = c(.3, 1))

dev.off()



# ------------------ simultaneous simulation plots --------------------

par(mfrow=c(1,1),mar=c(5,6,5,5),pty="m")
month<-"all"
fh<-paste(month,NL_days,mass,NL_shade,NL_food*10,"",sep="");fh
ttl<-paste("From ",daystart," + ",NL_days,"days ;","Weight =",mass,"g ; Shade =", NL_shade,"; Food =", NL_food,"; VTMIN", NL_T_b_min,"C")
pdf(paste("/Users/matthewmalishev/Documents/Manuscripts/Malishev and Kearney/Figures/Simulations/Tb plot/",fh,".pdf",sep=""),width=15,height=15,paper="a4r",title=ttl)
#plot.new()

# sep15

# make data points transparent
colvec = adjustcolor(c("black"), alpha = 0.3)
#col=colvec[sep15$Tb]

ticktime<-sep15$tick * (2 / 60 / 24) # convert to real days
sticktime<-shade.sep15$tick * 2 / 60 / 24 # convert to real days
with(sep15,plot(Tb~ticktime,type='l',las=1,lwd=2,lty=1,,col=colvec,xlim=c(0,max(ticktime)),ylim=c(0,45),xlab="",ylab="",axes=F))
shade.sep15<-sep15[sep15$shade %in% 1,] ; head(shade.sep15)
shade.sep15_Tb<-sep15[sep15$shade %in% 1,'Tb'] ; head(shade.sep15_Tb)
#with(results,points(shade.results$Tb~shade.results$tick, type='p', col="blue"))
with(sep15,points(shade.sep15$Tb~sticktime, type='p', col=colvec))


#nov1

colvec = adjustcolor(c("blue"), alpha = 0.3)

ticktime<-nov1$tick * (2 / 60 / 24) # convert to real days
sticktime<-shade.nov1$tick * 2 / 60 / 24 # convert to real days
par(new=T)
with(nov1,plot(Tb~ticktime,type='l',las=1,lwd=2,lty=1,col=colvec,xlim=c(0,max(ticktime)),ylim=c(0,45),xlab="",ylab="",axes=F))
shade.nov1<-nov1[nov1$shade %in% 1,] ; head(shade.nov1)
shade.nov1_Tb<-nov1[nov1$shade %in% 1,'Tb'] ; head(shade.nov1_Tb)
#with(results,points(shade.results$Tb~shade.results$tick, type='p', col="blue"))
with(nov1,points(shade.nov1$Tb~sticktime, type='p', col=colvec))

#dec15

colvec = adjustcolor(c("red"), alpha = 0.3)

ticktime<-dec15$tick * (2 / 60 / 24) # convert to real days
sticktime<-shade.dec15$tick * 2 / 60 / 24 # convert to real days
par(new=T)
with(dec15,plot(Tb~ticktime,type='l',las=1,lwd=2,lty=1,col=colvec,xlim=c(0,max(ticktime)),ylim=c(0,45),xlab="Days",ylab = expression(paste("Body temperature (" * degree,"C)"))))
shade.dec15<-dec15[dec15$shade %in% 1,] ; head(shade.dec15)
shade.dec15_Tb<-nov1[dec15$shade %in% 1,'Tb'] ; head(shade.dec15_Tb)
#with(results,points(shade.results$Tb~shade.results$tick, type='p', col="blue"))
with(dec15,points(shade.dec15$Tb~sticktime, type='p', col=colvec))
abline(h = c(NL_T_opt_l,NL_T_opt_u), col = "red", lty = 3)
text(0,30, "Activity range", col = "red", adj = c(.3, 1))



dev.off()  # use only with jpeg function above

hist(shade.results$Tb, main="Proportion of Tb when in shade", xlab="Tb (C) when in shade")
plot(shade.results$tick,shade.results$Tb, col="blue")




# ----------- home range plot ------------------


# get homerange ------------------
#X<-NLReport("[X] of turtle 0"); head(X)
#Y<-NLReport("[Y] of turtle 0"); head(Y)
#  turtles<-data.frame(X,Y)
#  who1<-rep(who,NL_ticks); who #
#who1<-rep(who,NL_ticks); who
#  turtle<-data.frame(ID = who1)
#  turtles<-cbind(turtles,turtle)
#
#

# need to make new variable for new spdf (movement paths)  ........ under progress
spdf<-SpatialPointsDataFrame(turtles[1:2], turtles[3]) # creates a spatial points data frame (adehabitatHR package)
homerange<-mcp(spdf,percent=100)

# draw-homerange ------------------
#calculate homerange points
temp <- slot(homerange,'polygons')[[which(slot(homerange,'data')$id == who)]]@Polygons[[1]]@coords
tempX<-temp[,1]
tempY<- temp[,2]
#tempXY<-cbind(tempX,tempY); tempXY<-as.list(tempXY)
tempXY<-NLCommand("set tempXY (map [list ?1 ?2]",tempX,tempY, ")" )
# Setting as.list makes no difference
#as.list(tempXY)

# The below command doesn't work yet, so the model calls the NL 'draw-homerange' procedure instead
#NLCommand("set tempX", tempX, "set tempY", tempY, "set tempXY", tempXY)
NLCommand("draw-homerange")

# hatch a turtle to draw the homerange boundary points
NLCommand("ask turtle 0 [hatch-homeranges 1", "hide-turtle","set color red]")
# draw the homerange
NLCommand("foreach tempXY [ask homeranges", "[move-to patch (item 0 ?) (item 1 ?)", "pd]]")
# close the homerange polygon
NLCommand("ask homeranges [let lastpoint first tempXY","move-to patch (item 0 lastpoint) (item 1 lastpoint)","]")

# to save-homerange ---------------

# change this with each sim
sep15_homerange<-homerange

month<-"all"
par(new=T)
par(pty="s")
#save the homerange polygon into a pdf
fh<-paste(month,NL_days,mass,NL_shade,NL_food*10,"_homerange","",sep="");fh
ttl<-paste("From ",daystart," + ",NL_days,"days ;","Weight =",mass,"g ; Shade =", NL_shade,"; Food =", NL_food,"; VTMIN", NL_T_b_min,"C")
#pdf(paste("/Users/matthewmalishev/Documents/Manuscripts/Malishev and Kearney/Figures/Simulations/Home range polygon/",fh,".pdf",sep=""),width=15,height=15,paper="a4r",title=ttl)
pdf(paste("/Users/matthewmalishev/Documents/Manuscripts/Malishev and Kearney/Figures/Simulations/Home range polygon/",fh,".pdf",sep=""))
jpeg(paste("/Users/matthewmalishev/Documents/Manuscripts/Malishev and Kearney/Figures/Simulations/Home range polygon/",fh,".jpeg",sep=""))

# plot polygon
colvec_h = adjustcolor(c("black"), alpha = 0.5)
plot(sep15_homerange,lty=1,bty="o",col=colvec_h,lwd=2,xlim=c(-25,25),ylim=c(-25,25),axes=F)
# plot animal tracks in polygon. NB: col = white. ........ under progress
plot(spdf,col=as.data.frame(spdf)[,1],pch=4,add=T)

# adding 2nd homerange polygon
colvec_h = adjustcolor(c("blue"), alpha = 0.2)
plot(nov1_homerange,lty=1,bty="o",col=colvec_h,lwd=2,xlim=c(-25,25),ylim=c(-25,25),axes=F,add=T)
# plot animal tracks in polygon. NB: col = white......... under progress
plot(spdf,col=as.data.frame(spdf)[,1],pch=4,add=T)

# adding 3rd homerange polygon
colvec_h = adjustcolor(c("red"), alpha = 0.2)
plot(dec15_homerange,lty=1,bty="o",col=colvec_h,lwd=2,xlim=c(-25,25),ylim=c(-25,25),axes=F,add=T)
# plot animal tracks in polygon. NB: col = white......... under progress
plot(spdf,col=as.data.frame(spdf)[,1],pch=4,add=T)
text(-16,5,labels="A",cex=1.5);text(0,5,labels="B",cex=1.5); text(14,5,labels="C",cex=1.5)

# use if axes = F
realaxis<-c("-50","-25","0","25","50")
axis(1,at=c(-25,-12.5,0,12.5,25),labels=realaxis)
axis(2,at=c(-25,-12.5,0,12.5,25),labels=realaxis)
par(new=T); plot(0, type="n",xlab="X",ylab="Y",axes=F)

box(which="plot")

dev.off()

#save a plot of the homerange area level into a pdf
pdf(paste("/Users/matthewmalishev/Documents/Manuscripts/Malishev and Kearney/Figures/Simulations/Home range polygon/",month,"_homearea",".pdf",sep=""))
plot(mcp.area(spdf,unout='m2'),colpol="blue",main=paste("Homerange polygon for ",NL_days," days with ",NL_shade," shade and ",NL_food," food patches"))
dev.off()

# original code. NB: using min/maxpcors doesn't work for exporting plot to pdf or jpeg
minpxcor<-NLReport("min-pxcor");maxpxcor<-NLReport("max-pxcor")
minpycor<-NLReport("min-pycor");maxpycor<-NLReport("max-pycor")
plot(nov1_homerange,lty=3,lwd=2,xlab="X", ylab="Y", xlim=c(minpxcor,maxpxcor), ylim=c(minpycor,maxpycor), axes=T)
par(new=T); plot(0, type="n",xlab="X coordinates",ylab="Y coordinates",axes=F)




#-----------------------------------------------------------
#------------------------ plot activity budget ------------
#-----------------------------------------------------------
library(data.table)

setwd("/Users/matthewmalishev/Documents/Manuscripts/Malishev and Kearney/Figures/Simulations/Activity budget/")
sep15_act<-read.csv("sep15/sep15.csv",header=T,sep=",",row.names=NULL,skip=20)

# sep15
names(sep15_act);head(sep15_act)
sfeed<-sep15_act[,1:2]
#setnames(sfeed,"x","X"); setnames(sfeed,"y","Y")
ssearchi<-sep15_act[,5:6]
# trans overlaps searching in 'diffr' plot (diff bw sep and dec activity)
#strans<-sep15_act[,9:10]
srest<-sep15_act[,13:14]

sep15_act<-as.data.frame(sep15_act)
activity_time<-sfeed$x*(2/60/24)

par(mfrow=c(1,1))

par(new=T)
activity<-sfeed
colvec_a = adjustcolor(c("black"), alpha = 1)
plot(activity_time,activity$y,
     type="s",
     lwd=1,
     ylab="",
     xlab="",
     xlim=c(0,max(activity_time)),
     ylim=c(0,5000),
     col=colvec_a,
     axes=F
     )


# use if axes = F
actaxis_time<-c("0","1","2","3","4","5","6","7")
actaxis_freq<-c("0","1000","2000","3000","4000","5000")
axis(1,at=c(0:7),labels=actaxis_time)
axis(2,at=c(0,1000,2000,3000,4000,5000),labels=actaxis_freq)
par(new=T); plot(0, type="n",xlab="Days",ylab="Frequency of activity state",axes=F)
box(which="plot")

# nov1
nov1_act<-read.csv("nov1/nov1.csv",header=T,sep=",",row.names=NULL,skip=20)

nfeed<-nov1_act[,1:2]
nsearchi<-nov1_act[,5:6]
#ntrans<-nov1_act[,9:10]
nrest<-nov1_act[,13:14]

# dec15
dec15_act<-read.csv("dec15/dec15.csv",header=T,sep=",",row.names=NULL,skip=20)

dfeed<-dec15_act[,1:2]
dsearchi<-dec15_act[,5:6]
#dtrans<-dec15_act[,9:10]
drest<-dec15_act[,13:14]


# plot difference between time spent in activity states in early vs late season
par(mfrow=c(1,1),mar=c(5,6,5,5),pty="m")
month<-"diffr"
fh<-paste(month,NL_days,mass,NL_shade,NL_food*10,"_act","",sep="");fh
pdf(paste("/Users/matthewmalishev/Documents/Manuscripts/Malishev and Kearney/Figures/Simulations/Activity budget/",fh,".pdf",sep=""))
jpeg(paste("/Users/matthewmalishev/Documents/Manuscripts/Malishev and Kearney/Figures/Simulations/Activity budget/",fh,".jpeg",sep=""))

# plotting order
#rest
#searchi
#feed

par(mfrow=c(1,1),mar=c(5,6,5,5),pty="m")
diffra<-"rest"
diffr<-eval(parse(text=paste("s",diffra,"-","d",diffra,sep="")))

par(new=T)
activity<-diffr
colvec_a = adjustcolor(c("dark green"), alpha = 1)
plot(activity_time,activity$y,
     type="s",
     lwd=1,
     ylab="",
     xlab="",
     xlim=c(0,7),
     ylim=c(-400,400),
     col=colvec_a,
     axes=F
     )

# find y==0 value (pivotal point when animal switches activity)  ........ under progress
diffr23<-subset(diffr,subset=srest$x * (2 / 60 / 24)==c(2:3))
diffrz<-subset(diffr23,subset=diffr23$y==0)

# use if axes = F
actaxis_time<-c("0","1","2","3","4","5","6","7")
actaxis_freq<-c("-400","-200","0","200","400")
axis(1,at=c(0:7),labels=actaxis_time,las=1)
axis(2,at=c(-400,-200,0,200,400),labels=actaxis_freq,las=1)
mtext("Days",1,line=3)
mtext("Difference in time spent within activity state",2,line=4)
abline(v=2.41111,col="red",lty=3)
box(which="plot")

dev.off()

#------------- reserve plot

par(mfrow=c(1,1))
#Netlogo outputs
setwd("/Users/matthewmalishev/Documents/Manuscripts/Malishev and Kearney/Figures/Simulations/Reserve plot")

setwd("/Users/matthewmalishev/Documents/Manuscripts/Malishev and Kearney/Figures/Simulations/")


reserve<-read.table("sep157800100100_reserve.csv",header=T,sep=",",row.names=NULL,skip=17)
# change for each sim month
sep15reserve<-reserve

reserve<-reserve[,1:2]
reserve<-as.data.frame(reserve)
reserve$Time<-reserve$x * (2 / 60 / 24)
plot(reserve$Time,reserve$y, main="Reserve level (J)",xlim=c(0,max(reserve$Time)),xlab="Days",ylab="Reserve level (J)", col="red",las=1,type="l")


#-----------------------------------------------------------
#------------------------ plot movement paths --------------
#-----------------------------------------------------------

# get all food and shade patches
foodh <- NLGetAgentSet(c("pxcor","pycor"),
                       "patches with [pcolor = green]",
                       as.data.frame=TRUE,
                       # df.col.names=c("x","y")
                       )

foodl <- NLGetAgentSet(c("pxcor","pycor"),
                       "patches with [pcolor = yellow]",
                       as.data.frame=TRUE,
                       )

shadep1 <- NLGetAgentSet(c("pxcor","pycor"),
                         "patches with [pcolor = black]",
                         as.data.frame=TRUE,
                         )

shadep2 <- NLGetAgentSet(c("pxcor","pycor"),
                         "patches with [pcolor = black + 2]",
                         as.data.frame=TRUE,
                         )

# remove NA's
foodh<-na.omit(foodh)
foodl<-na.omit(foodl)

shadeall<-rbind(shadep1,shadep2)
foodall<-rbind(foodh,foodl)


#x.minmax <- NLReport("(list min-pxcor max-pxcor)")
#y.minmax <- NLReport("(list min-pycor max-pycor)")

# retrieve turtle location
#liz <- NLGetAgentSet(c("xcor","ycor"),
# "turtles","pd",
# as.data.frame=TRUE,
#)
#par(new=T)
# plot(liz, xlim=x.minmax, ylim=y.minmax,col="red", pch=19)

# ------------ formatting data

setwd("/Users/matthewmalishev/Documents/Manuscripts/Malishev and Kearney/Figures/Simulations/Movement trajectories/")

# -----sep15
spatialplot<-read.table("sep15/sep15.csv",header=T,sep=",",skip=19)
spatialplot<-as.data.frame(spatialplot); head(spatialplot)
names(spatialplot)
# format data
spen_move<-spatialplot[,1:2];head(spen_move)
sfeed_move<-spatialplot[,5:6];head(sfeed_move)
ssearch_move<-spatialplot[,9:10];head(ssearch_move)
srest_move<-spatialplot[,13:14]; head(srest_move)
#round off food and shade data points
sfeed_move<-round(sfeed_move,0)
srest_move<-round(srest_move,0)

# ------nov1
spatialplot<-read.table("nov1/nov1.csv",header=T,sep=",",skip=19)
spatialplot<-as.data.frame(spatialplot); head(spatialplot)
names(spatialplot)
# format data
npen_move<-spatialplot[,1:2];head(npen_move)
nfeed_move<-spatialplot[,5:6];head(nfeed_move)
nsearch_move<-spatialplot[,9:10];head(nsearch_move)
nrest_move<-spatialplot[,13:14]; head(nrest_move)
#round off food and shade data points
nfeed_move<-round(nfeed_move,0)
nrest_move<-round(nrest_move,0)

# -----dec15
spatialplot<-read.table("dec15/dec15.csv",header=T,sep=",",skip=19)
spatialplot<-as.data.frame(spatialplot); head(spatialplot)
names(spatialplot)
# format data
dpen_move<-spatialplot[,1:2];head(dpen_move)
dfeed_move<-spatialplot[,5:6];head(dfeed_move)
dsearch_move<-spatialplot[,9:10];head(dsearch_move)
drest_move<-spatialplot[,13:14]; head(drest_move)
#round off food and shade data points
dfeed_move<-round(dfeed_move,0)
drest_move<-round(drest_move,0)


# remove NA's
sfeed_move<-na.omit(sfeed_move)


#subset feeding data by food patch type (high or low)
green<-subset(spatialplot[,5:7],color.1==55,select=c(x.1,y.1))
green<-green[,1:2]
yellow<-subset(spatialplot[,5:7],color.1==45,select=c(x.1,y.1))
yellow<-yellow[,1:2]
green<-round(green,0) ;yellow<-round(yellow,0)


# ------------------ plotting all food and shade patches

# plot food patches
# colvec_m = adjustcolor(c("yellow"), alpha = 0.8)
plot(foodl, xlim=c(-25,25), ylim=c(-25,25),pch=22,cex=2,col="grey",bg="yellow",xlab="",ylab="",axes=F)
par(new=T)
plot(foodh, xlim=c(-25,25), ylim=c(-25,25),pch=22, cex=2,col="grey",bg="dark green",xlab="",ylab="",axes=F)
par(new=T) #if plotting food patches above

# plot all shade patches
par(new=T)
plot(shadeall, xlim=c(-25,25), ylim=c(-25,25),pch=22, cex=2,col="grey",bg="black",xlab="",ylab="",axes=F)
par(new=T)

# ------------------- plotting movement paths

# change input to reflect simulation month
par(pty="s")
month<-"sep15"
mr<-"s"
feed_cex<-eval(parse(text=paste("rbind(",mr,"feed_move$x,",mr,"feed_move$y)",sep="")))
rest_cex<-eval(parse(text=paste("rbind(",mr,"rest_move$x,",mr,"rest_move$y)",sep="")))

# feed_cex<-rbind(sfeed_move$x,sfeed_move$y);feed_cex
# rest_cex<-rbind(srest_move$x,srest_move$y)

# remove duplicated elements of matrix
#rest_cex<-rest_cex[, !duplicated(t(rest_cex))]

#export plot to file
fh<-paste(month,NL_days,mass,NL_shade,NL_food*10,"_move","",sep="");fh
pdf(paste("/Users/matthewmalishev/Documents/Manuscripts/Malishev and Kearney/Figures/Simulations/Movement trajectories/",fh,".pdf",sep=""))
jpeg(paste("/Users/matthewmalishev/Documents/Manuscripts/Malishev and Kearney/Figures/Simulations/Movement trajectories/",fh,".jpeg",sep=""))

# feeding
# plot by food level. Uses opposite plot colours.
colvec_m = adjustcolor(c("yellow"), alpha = 0.4)
plot(green,pch=".",cex=feed_cex,xlim=c(-25,25), ylim=c(-25,25), xlab="",ylab="",col=colvec_m,axes=F)
par(new=T)
colvec_m = adjustcolor(c("green"), alpha = 0.4)
plot(yellow,pch=".",cex=feed_cex,xlim=c(-25,25), ylim=c(-25,25), xlab="",ylab="",col=colvec_m,axes=F)
#or plot without food level
par(new=T)
colvec_m = adjustcolor(c("brown"), alpha = 0.4)
plot(sfeed_move,pch=".",cex=feed_cex,xlim=c(-25,25), ylim=c(-25,25), xlab="",ylab="",col=colvec_m,axes=F)

# plot(sfeed_move,pch=".",cex=c(sfeed_move$x),xlim=c(-25,25), ylim=c(-25,25), xlab="",ylab="",col=colvec_m,axes=F)
# par(new=T)
# plot(sfeed_move,pch=".",cex=c(sfeed_move$y),xlim=c(-25,25), ylim=c(-25,25),xlab="",ylab="", col=colvec_m,axes=F)

# searching
par(new=T)
colvec_m = adjustcolor(c("blue"), alpha = 0.5)
plot(ssearch_move,type="l",xlim=c(-25,25), ylim=c(-25,25), xlab="",ylab="",col=colvec_m,axes=F)

# resting
par(new=T)
colvec_m = adjustcolor(c("black"), alpha = 0.8)
plot(srest_move,pch=".",cex=rest_cex,xlim=c(-25,25), ylim=c(-25,25),xlab="",ylab="",col=colvec_m,axes=F)

#plot(srest_move,pch=".",cex=srest_move$y,xlim=c(-25,25), ylim=c(-25,25),xlab="",ylab="",col=colvec_m,axes=F)
#par(new=T)
#plot(srest_move,pch=".",cex=srest_move$x,xlim=c(-25,25), ylim=c(-25,25),xlab="",ylab="",col=colvec_m,axes=F)

# turtle start position
turtle_start<-ssearch_move[1,]
#turtle_start<-data.frame(-8,19)
par(new=T)
plot(turtle_start,pch=8,cex=3,col="pink",xlim=c(-25,25), ylim=c(-25,25),xlab="",ylab="",axes=F)

# use if axes = F
realaxis<-c("-50","-25","0","25","50")
axis(1,at=c(-25,-12.5,0,12.5,25),labels=realaxis)
axis(2,at=c(-25,-12.5,0,12.5,25),labels=realaxis)
par(new=T); plot(0, type="n",xlab="X",ylab="Y",axes=F)
box(which="plot")


dev.off()



# cool random plot
colvec_m = adjustcolor(c("pink","black"), alpha = 0.8)
plot(dsearch_move$x, dsearch_move$y,pch=19,cex=dsearch_move$x,xlim=c(-25,25), ylim=c(-25,25), col=colvec_m,xlab="X coord",ylab="Y coord")




#-----------------------------------------------------------
#------------------------ plot DEB RESULTS -----------------
#-----------------------------------------------------------


colnames(all_results_deb)<-c("hour","E_pres","V_pres","E_H_pres","q_pres","hs_pres","surviv_pres","ms_pres","cumrepro","cumbatch","O2FLUX","CO2FLUX","MLO2","GH2OMET","DEBQMET","DRYFOOD","FAECES","NWASTE","wetgonad","wetstorage","wetfood","wetmass","gutfreemass")
all_results_deb<-as.data.frame(all_results_deb)
plot(all_results_deb$wetmass~all_results_deb$hour,type='l') # Plot mass over time
plot((all_results_deb$gutfreemass-all_results_deb$wetgonad)~all_results_deb$hour,type='l') # plot mass excluding reproduction buffer
plot(all_results_deb$GH2OMET~all_results_deb$hour,type='l') # plot metabolic water over time
plot(all_results_deb$E_pres~all_results_deb$hour,type='l') # Plot reserve density over time







# Lab

VAR<-"ms_pres";VAR<-as.vector(VAR)
par(new=TRUE)
plot(all_results_deb$hour, paste("all_results_deb$",VAR,sep=""),type="l")
eval(parse(text=paste("plot(all_results_deb$hour,all_results_deb$",VAR,",type='l')",sep="")))
VAR
paste("all_results_deb$","ms_pres",sep="")

# Output variables (from str(all_results))
# $ hour       : num  1 2 3 4 5 6 7 8 9 10 ...
# $ E_pres     : num  6813 6806 6806 6806 6806 ...
# $ V_pres     : num  299 300 300 300 300 ...
# $ E_H_pres   : num  137501 137501 137501 137501 137501 ...
# $ q_pres     : num  0.00 2.26e-16 4.52e-16 6.77e-16 9.03e-16 ...
# $ hs_pres    : num  0.00 0.00 2.26e-16 6.77e-16 1.35e-15 ...
# $ surviv_pres: num  1 1 1 1 1 ...
# $ ms_pres    : num  0 334288 334299 334311 334322 ...
# $ cumrepro   : num  0 -11.2 0 -11 0 ...
# $ cumbatch   : num  0 220 439 659 879 ...
# $ O2FLUX     : num  148 177 177 177 177 ...
# $ CO2FLUX    : num  117 140 140 140 140 ...
# $ MLO2       : num  87.7 105 105 105 105 ...
# $ GH2OMET    : num  0.05 0.0599 0.0599 0.0599 0.0599 ...
# $ DEBQMET    : num  0.168 0.18 0.18 0.18 0.18 ...
# $ DRYFOOD    : num  0 0.124 0.124 0.124 0.124 ...
# $ FAECES     : num  0 0.0136 0.0136 0.0136 0.0136 ...
# $ NWASTE     : num  0.0212 0.0254 0.0254 0.0254 0.0254 ...
# $ wetgonad   : num  0 0.0284 0.0598 0.0883 0.1197 ...
# $ wetstorage : num  278 278 278 278 278 ...
# $ wetfood    : num  18.9 18.9 18.9 18.9 18.9 ...
# $ wetmass    : num  596 596 596 596 596 ...
# $ gutfreemass: num  577 577 577 577 577 ...






#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------

# Experimental R lab

# 27-1-15
x <- rchisq(100, df = 4)
qqplot(x, qchisq(ppoints(x), df = 4)); abline(0, 1, col = 2, lty = 2)


points(shade.results$Tb~shade.results$tick, type='p', col="blue")

# Replace file with new file name if it already exists

current.file<-'NL turtle Tb with shade.jpg'
if (exists(current.file))
{
  for (i in 1:100){
    replace.file<-i
    new.file<-c(current.file[replace.file],current.file[replace.file+1])
  }
  
  jpeg(paste('NL turtle Tb with shade','new.file','.jpg')
       dev.off()
}


# Calling T_opt_lower and T_opt_upper values from NL to create dynamics activity range line in plot

abline(h = c(paste(NLGetAgentSet("T_opt_lower","turtles")),NLGetAgentSet("T_opt_upper","turtles")), col = "red", lty = 3)


# ------ plot function -------

month<-function(sim) {
  sim<-newresults
  par(new=T)
  ticktime<-sim$tick * (2 / 60 / 24) # convert to real days
  sticktime<-shade.sim$tick * 2 / 60 / 24 # convert to real days
  # Plot with title
  # with(results,plot(Tb~ticktime,type='l',las=1,xlab="Days",ylab = expression(paste("Body temperature (" * degree,"C)")),main=paste("From ",daystart," + ",NL_days,"days ;","Weight =",mass,"g ; Shade =", NL_shade,"; Food =", NL_food,"; VTMIN", NL_T_b_min,"C")))
  with(sim,plot(Tb~ticktime,type='l',las=1,lwd=2,lty=2,col="grey",xlab="Days",ylab = expression(paste("Body temperature (" * degree,"C)"))))
  shade.sim<-sim[sim$shade %in% 1,]
  shade.sim_Tb<-sim[sim$shade %in% 1,'Tb'] ; head(shade.sim_Tb)
  #with(results,points(shade.results$Tb~shade.results$tick, type='p', col="blue"))
  with(sim,points(shade.sim$Tb~sticktime, type='p', col="lightblue"))
  #with(results,points(shade*20~tick,type='p',col="red"))
  abline(h = c(NL_T_opt_l,NL_T_opt_u), col = "red", lty = 3)
  text(0,28, "Activity range", col = "red", adj = c(.3, 1))
}

# -------- replace old results -----------
if (exists("results")){
  newones<-sep15plus
}

ticktime<-newones$tick * (2 / 60 / 24) # convert to real days
shade.newones<-newones[newones$shade %in% 1,]
shade.newones_Tb<-newones[newones$shade %in% 1,'Tb'] ; head(shade.newones_Tb)
sticktime<-shade.newones$tick * 2 / 60 / 24 # convert to real days

# Plot with title
# with(results,plot(Tb~ticktime,type='l',las=1,xlab="Days",ylab = expression(paste("Body temperature (" * degree,"C)")),main=paste("From ",daystart," + ",NL_days,"days ;","Weight =",mass,"g ; Shade =", NL_shade,"; Food =", NL_food,"; VTMIN", NL_T_b_min,"C")))
with(newones,plot(Tb~ticktime,type='l',las=1,lwd=2,lty=2,col="grey",xlab="Days",ylab = expression(paste("Body temperature (" * degree,"C)"))))
#with(results,points(shade.results$Tb~shade.results$tick, type='p', col="blue"))
with(newones,points(shade.newones$Tb~sticktime, type='p', col="lightblue"))
#with(results,points(shade*20~tick,type='p',col="red"))
abline(h = c(NL_T_opt_l,NL_T_opt_u), col = "red", lty = 3)
text(0,28, "Activity range", col = "red", adj = c(.3, 1))

# --------------- movement path plots -----------------
# comparing values between sfeed_move and foodh to find high/low food xy coords in sfeed_move and plot on map
# 1.
sfeed_moveh<-data.frame(sfeed_move, foodh[match(sfeed_move, foodh), 2])
sfeed_movel<-data.frame(sfeed_move, foodl[match(sfeed_move, foodl), 2])
names(foodh)
head(sfeed_movel)

#2.
library(compare)
comparisonh <- compare(sfeed_move,foodh,allowAll=TRUE)
sfeed_moveh<-comparisonh$tM
# make df
sfeed_moveh<-as.data.frame(sfeed_moveh); colnames(sfeed_moveh)<-c("x","y")
comparisonh$tM

comparisonl <- compare(sfeed_move,foodl,allowAll=TRUE)
sfeed_movel<-comparisonl$tM
# make df
sfeed_movel<-as.data.frame(sfeed_movel); colnames(sfeed_movel)<-c("x","y")
comparisonl$tM

diff_h <-data.frame(lapply(1:ncol(sfeed_move),function(i)setdiff(sfeed_move[,i],comparisonh$tM[,i])))
colnames(diff_h) <- colnames(sfeed_move)

diff_l <-data.frame(lapply(1:ncol(sfeed_move),function(i)setdiff(sfeed_move[,i],comparisonl$tM[,i])))
colnames(diff_l) <- colnames(sfeed_move)

#2.1.
comparisonh <- compare(foodh,sfeed_move,allowAll=TRUE)
comparisonh$tM; colnames(comparisonh$tM)<-c("x","y")

comparisonl <- compare(foodl,sfeed_move,allowAll=TRUE)
comparisonl$tM; colnames(comparisonl$tM)<-c("x","y")

diff_h <-data.frame(lapply(1:ncol(sfeed_move),function(i)setdiff(sfeed_move[,i],comparisonh$tM[,i])))
colnames(diff_h) <- colnames(sfeed_move)

diff_l <-data.frame(lapply(1:ncol(sfeed_move),function(i)setdiff(sfeed_move[,i],comparisonl$tM[,i])))
colnames(diff_l) <- colnames(sfeed_move)

#3.
library(sqldf)
# select rows of sfeed_move not in foodh
a1NotIna2_h  <- sqldf('SELECT * FROM sfeed_move EXCEPT SELECT * FROM foodh')
a1NotIna2_l  <- sqldf('SELECT * FROM sfeed_move EXCEPT SELECT * FROM foodl')

# select rows from sfeed_move that also appear in foodh
a1Ina2_h  <- sqldf('SELECT * FROM sfeed_move INTERSECT SELECT * FROM foodh')
a1Ina2_l  <- sqldf('SELECT * FROM sfeed_move INTERSECT SELECT * FROM foodl')

# matching plot functions
# high and low food plot
#1.
par(new=T)
plot(sfeed_moveh,pch=".",cex=c(feed_cex),xlim=c(-25,25), ylim=c(-25,25), xlab="",ylab="",col="green",axes=F)
par(new=T)
plot(sfeed_movel,pch=".",cex=c(feed_cex),xlim=c(-25,25), ylim=c(-25,25), xlab="",ylab="",col="yellow",axes=F)

#2.
par(new=T)
plot(diff_h,pch=".",cex=c(feed_cex),xlim=c(-25,25), ylim=c(-25,25), xlab="",ylab="",col="green",axes=F)
par(new=T)
plot(diff_l,pch=".",cex=c(feed_cex),xlim=c(-25,25), ylim=c(-25,25), xlab="",ylab="",col="yellow",axes=F)

#3.
par(new=T)
plot(a1NotIna2_h,pch=".",cex=c(feed_cex),xlim=c(-25,25), ylim=c(-25,25), xlab="",ylab="",col="green",axes=F)
par(new=T)
plot(a1NotIna2_l,pch=".",cex=c(feed_cex),xlim=c(-25,25), ylim=c(-25,25), xlab="",ylab="",col="yellow",axes=F)


# --------------- exporting plots ------------------------
NLCommand("setup")
# need to somehow change working dir
NLCommand("export-plot \"Spatial coordinates of transition between activity states\" \"Users/matthewmalishev/Documents/Manuscripts/Malishev and Kearney/Figures/Simulations/spatialplot.csv\"")

#spatial plot
NLCommand("export-plot \"Spatial coordinates of transition between activity states\" \"Simulations/spatialplot.csv\"")
# reserve plot
NLCommand("export-plot \"Reserve level and starvation reserve over time\" \"Simulations/reserveplot.csv\"")
#temp plot
NLCommand("export-plot \"Body temperature (T_b)\" \"Simulations/tempplot.csv\"")
# activity budget
NLCommand("export-plot \"Global time budget\" \"Simulations/activitybudget.csv\"")
# world view
NLCommand("export-view \"Simulations/worldview.png\"")


# clear all plots
graphics.off()



