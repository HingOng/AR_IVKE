; Numerical analyses of IVT, VKE, IVKE, and tendency terms from MERRA2 data
; Here is a template for executing this script in a PBS job array:
;#PBS -J 20190101-20191231
;### Run program
;ncl 'date='$PBS_ARRAY_INDEX
;    'f_qdt=addfile("MERRA2_400.tavg3_3d_qdt_Np.'$PBS_ARRAY_INDEX'.nc4","r")'
;    'f_udt=addfile("MERRA2_400.tavg3_3d_udt_Np.'$PBS_ARRAY_INDEX'.nc4","r")'
;    'f_asm=addfile("MERRA2_400.inst3_3d_asm_Np.'$PBS_ARRAY_INDEX'.nc4","r")'
;   'f_asm1=addfile("MERRA2_400.inst3_3d_asm_Np.'$((PBS_ARRAY_INDEX+1))'.nc4","r")'
;    f_asma=addfile("MERRA2_400.tavg3_3d_asm_Nv.'$PBS_ARRAY_INDEX'.nc4","r")' ARanalysisMERRA2.ncl
; Hing Ong, 2 Oct 2024

; Parameter setting
g   =    9.81
pB  = 1000.
pT  =  200.
pTp =  150.

; Add input files
if (ismissing(f_qdt)) then
    exit
end if

print("Reading begins")
time_avg = f_qdt->time
n_tim_avg= dimsizes(time_avg)
lev      = f_asm->lev     (          {pB:pTp}    )
n_lev    = dimsizes(lev)
lat      = f_asm->lat
n_lat    = dimsizes(lat)
lon      = f_asm->lon
n_lon    = dimsizes(lon)
time     = f_asm->time(0:n_tim_avg-1)
n_tim    = dimsizes(time)
ps0      = f_asm->PS      (0:n_tim-1,         :,:)
q0       = f_asm->QV      (0:n_tim-1,{pB:pTp},:,:)
u0       = f_asm->U       (0:n_tim-1,{pB:pTp},:,:)
v0       = f_asm->V       (0:n_tim-1,{pB:pTp},:,:)
h0       = f_asm->H       (0:n_tim-1,{pB:pTp},:,:)
omega0   = f_asm->OMEGA   (0:n_tim-1,{pB:pTp},:,:)
ps1      = f_asm1->PS     (0:1      ,         :,:)
q1       = f_asm1->QV     (0:1,      {pB:pTp},:,:)
u1       = f_asm1->U      (0:1,      {pB:pTp},:,:)
v1       = f_asm1->V      (0:1,      {pB:pTp},:,:)
h1       = f_asm1->H      (0:1,      {pB:pTp},:,:)
omega1   = f_asm1->OMEGA  (0:1,      {pB:pTp},:,:)
ps       = array_append_record(   ps0,   ps1,0)
q        = array_append_record(    q0,    q1,0)
u        = array_append_record(    u0,    u1,0)
v        = array_append_record(    v0,    v1,0)
h        = array_append_record(    h0,    h1,0)
omega    = array_append_record(omega0,omega1,0)
PL       = f_asma->PL     (0:n_tim-1,::-1,:,:)
PS       = f_asma->PS     (0:n_tim-1,     :,:)
Q0       = f_asma->QV     (0:n_tim-1,::-1,:,:)
U0       = f_asma->U      (0:n_tim-1,::-1,:,:)
V0       = f_asma->V      (0:n_tim-1,::-1,:,:)
H0       = f_asma->H      (0:n_tim-1,::-1,:,:)
OMEGA0   = f_asma->OMEGA  (0:n_tim-1,::-1,:,:)
DQVDTANA = f_qdt->DQVDTANA(0:n_tim-1,{pB:pTp},:,:)
DQVDTDYN = f_qdt->DQVDTDYN(0:n_tim-1,{pB:pTp},:,:)
DQVDTMST = f_qdt->DQVDTMST(0:n_tim-1,{pB:pTp},:,:)
DQVDTTRB = f_qdt->DQVDTTRB(0:n_tim-1,{pB:pTp},:,:)
DQVDTCHM = f_qdt->DQVDTCHM(0:n_tim-1,{pB:pTp},:,:)
DUDTANA  = f_udt->DUDTANA (0:n_tim-1,{pB:pTp},:,:)
DUDTDYN  = f_udt->DUDTDYN (0:n_tim-1,{pB:pTp},:,:)
DUDTMST  = f_udt->DUDTMST (0:n_tim-1,{pB:pTp},:,:)
DUDTTRB  = f_udt->DUDTTRB (0:n_tim-1,{pB:pTp},:,:)
DUDTGWD  = f_udt->DUDTGWD (0:n_tim-1,{pB:pTp},:,:)
DVDTANA  = f_udt->DVDTANA (0:n_tim-1,{pB:pTp},:,:)
DVDTDYN  = f_udt->DVDTDYN (0:n_tim-1,{pB:pTp},:,:)
DVDTMST  = f_udt->DVDTMST (0:n_tim-1,{pB:pTp},:,:)
DVDTTRB  = f_udt->DVDTTRB (0:n_tim-1,{pB:pTp},:,:)
DVDTGWD  = f_udt->DVDTGWD (0:n_tim-1,{pB:pTp},:,:)
print("Reading ends")

print("Calculation begins")

; Instantaneous steps
p        = dble2flt(lev)*100.
pB       = pB*100.
pT       = pT*100.
LAT      = dble2flt(lat)
LON      = dble2flt(lon)
q2       = q ; copy attributes
k        = q ; copy attributes
qu       = q ; copy attributes
qv       = q ; copy attributes
vke      = q ; copy attributes
VKE      = DQVDTANA ; copy attributes
VKETEND  = DQVDTANA ; copy attributes
DOKE     = DQVDTANA ; copy attributes
HAKE     = DQVDTANA ; copy attributes
VAKE     = DQVDTANA ; copy attributes
PEKE     = DQVDTANA ; copy attributes
MOKE     = DQVDTANA ; copy attributes
TOKE     = DQVDTANA ; copy attributes
GOKE     = DQVDTANA ; copy attributes
DOV      = DQVDTANA ; copy attributes
HAV      = DQVDTANA ; copy attributes
VAV      = DQVDTANA ; copy attributes
POV      = DQVDTANA ; copy attributes
TOV      = DQVDTANA ; copy attributes
COV      = DQVDTANA ; copy attributes
AOVKE    = DQVDTANA ; copy attributes
HCV      = DQVDTANA ; copy attributes
VCV      = DQVDTANA ; copy attributes
q2       = q*q
k        = (u*u+v*v)*0.5
qu       = q*u
qv       = q*v
vke      = q2*k
iqu      = vibeta(p,qu (time|:,lat|:,lon|:,lev|:),1,ps,pB,pT) /g
iqv      = vibeta(p,qv (time|:,lat|:,lon|:,lev|:),1,ps,pB,pT) /g
IVKE     = vibeta(p,vke(time|:,lat|:,lon|:,lev|:),1,ps,pB,pT) /g
IVT      = sqrt(iqu*iqu+iqv*iqv)

; Time average steps
dUdp0    = center_finite_diff_n(U0       ,PL,False,0,1)
WdUdp0   = OMEGA0*dUdp0
dVdp0    = center_finite_diff_n(V0       ,PL,False,0,1)
WdVdp0   = OMEGA0*dVdp0
dQdp0    = center_finite_diff_n(Q0       ,PL,False,0,1)
WdQdp0   = OMEGA0*dQdp0
dQWdp0   = center_finite_diff_n(Q0*OMEGA0,PL,False,0,1)

gradKx0  = Q0
gradKy0  = Q0
gradQx0  = Q0
gradQy0  = Q0
gradsf((U0*U0+V0*V0)*0.5,gradKx0,gradKy0)
gradsf( Q0              ,gradQx0,gradQy0)
UgradQx0 = U0*gradQx0
VgradQy0 = V0*gradQy0

divQ0    = uv2dvF(Q0*U0,Q0*V0)

Q        = linint1_n_Wrap(PL,Q0    ,False,p,0,1)
U        = linint1_n_Wrap(PL,U0    ,False,p,0,1)
V        = linint1_n_Wrap(PL,V0    ,False,p,0,1)
H        = linint1_n_Wrap(PL,H0    ,False,p,0,1)
Q2       = Q*Q
K        = (U*U+V*V)*0.5
VKE      = Q2*K
Q2U      = Q2*U
Q2V      = Q2*V
KQ2      = 2.*K*Q

WdUdp    = linint1_n_Wrap(PL,WdUdp0,False,p,0,1)
WdVdp    = linint1_n_Wrap(PL,WdVdp0,False,p,0,1)
WdQdp    = linint1_n_Wrap(PL,WdQdp0,False,p,0,1)
dQWdp    = linint1_n_Wrap(PL,dQWdp0,False,p,0,1)
VAKE     = -Q2U*WdUdp-Q2V*WdVdp
VAV      = -KQ2*WdQdp
VCV      = -KQ2*dQWdp

PHI      = H*g
gradPHI  = grad_latlon_cfd(PHI,LAT,LON,True,False)
gradKx   = linint1_n_Wrap(PL, gradKx0,False,p,0,1)
gradKy   = linint1_n_Wrap(PL, gradKy0,False,p,0,1)
UgradQx  = linint1_n_Wrap(PL,UgradQx0,False,p,0,1)
VgradQy  = linint1_n_Wrap(PL,VgradQy0,False,p,0,1)
divQ     = linint1_n_Wrap(PL,   divQ0,False,p,0,1)
PEKE     = -  Q2U*gradPHI[1]-  Q2V*gradPHI[0]
HAKE     = -  Q2U*gradKx    -  Q2V*gradKy  
HAV      = - KQ2*(UgradQx + VgradQy)
HCV      = - KQ2*  divQ

vke_missing_filled = linmsg_n(VKE,-1,1)

; Tendency from instantaneous steps
;VKETEND  = ( vke(1:n_tim,:,:,:)- vke(0:n_tim-1,:,:,:))/10800.
DQDT     = (   q(1:n_tim,:,:,:)-   q(0:n_tim-1,:,:,:))/10800.
DUDT     = (   u(1:n_tim,:,:,:)-   u(0:n_tim-1,:,:,:))/10800.
DVDT     = (   v(1:n_tim,:,:,:)-   v(0:n_tim-1,:,:,:))/10800.
VKETEND  = KQ2*DQDT + Q2U*DUDT + Q2V*DVDT
IVKETEND = (IVKE(1:n_tim,  :,:)-IVKE(0:n_tim-1,  :,:))/10800.
pstend   = (  ps(1:n_tim,  :,:)-  ps(0:n_tim-1,  :,:))/10800.
SPTE     = pstend*vke_missing_filled(:,0,:,:)/g

; Tendency steps
MOKE     = Q2U*DUDTMST+Q2V*DVDTMST
TOKE     = Q2U*DUDTTRB+Q2V*DVDTTRB
GOKE     = Q2U*DUDTGWD+Q2V*DVDTGWD
DOKE     = Q2U*DUDTDYN+Q2V*DVDTDYN
AOKE     = Q2U*DUDTANA+Q2V*DVDTANA
POV      = KQ2*DQVDTMST
TOV      = KQ2*DQVDTTRB
COV      = KQ2*DQVDTCHM
DOV      = KQ2*DQVDTDYN
AOV      = KQ2*DQVDTANA
AOVKE    = AOKE+AOV

IVKEAVG  = vibeta(p,VKE (time|:,lat|:,lon|:,lev|:),1,PS,pB,pT) /g
IPEKE    = IVKEAVG

IDOKE    = vibeta(p,DOKE (time|:,lat|:,lon|:,lev|:),1,PS,pB,pT) /g
IHAKE    = vibeta(p,HAKE (time|:,lat|:,lon|:,lev|:),1,PS,pB,pT) /g
IVAKE    = vibeta(p,VAKE (time|:,lat|:,lon|:,lev|:),1,PS,pB,pT) /g
IPEKE(:,1:n_lat-2,:) = vibeta(p,PEKE (time|:,lat|1:n_lat-2,lon|:,lev|:),1,PS(:,1:n_lat-2,:),pB,pT) /g
IMOKE    = vibeta(p,MOKE (time|:,lat|:,lon|:,lev|:),1,PS,pB,pT) /g
ITOKE    = vibeta(p,TOKE (time|:,lat|:,lon|:,lev|:),1,PS,pB,pT) /g
IGOKE    = vibeta(p,GOKE (time|:,lat|:,lon|:,lev|:),1,PS,pB,pT) /g
IDOV     = vibeta(p,DOV  (time|:,lat|:,lon|:,lev|:),1,PS,pB,pT) /g
IHAV     = vibeta(p,HAV  (time|:,lat|:,lon|:,lev|:),1,PS,pB,pT) /g
IVAV     = vibeta(p,VAV  (time|:,lat|:,lon|:,lev|:),1,PS,pB,pT) /g
IHCV     = vibeta(p,HCV  (time|:,lat|:,lon|:,lev|:),1,PS,pB,pT) /g
IVCV     = vibeta(p,VCV  (time|:,lat|:,lon|:,lev|:),1,PS,pB,pT) /g
IPOV     = vibeta(p,POV  (time|:,lat|:,lon|:,lev|:),1,PS,pB,pT) /g
ITOV     = vibeta(p,TOV  (time|:,lat|:,lon|:,lev|:),1,PS,pB,pT) /g
ICOV     = vibeta(p,COV  (time|:,lat|:,lon|:,lev|:),1,PS,pB,pT) /g
IAOVKE   = vibeta(p,AOVKE(time|:,lat|:,lon|:,lev|:),1,PS,pB,pT) /g
I_VKETEND= vibeta(p,VKETEND(time|:,lat|:,lon|:,lev|:),1,PS,pB,pT) /g

print("Calculation ends")

print("Writing begins")

; Instantaneous steps
outN = "./analysis_AR_inst"+date+".nc"
outF = addfile(outN ,"c")

dim_names     = (/"time","lat","lon"/)
dim_sizes     = (/n_tim,n_lat,n_lon/)
dim_unlimited = (/True,False,False/)
filedimdef(outF,dim_names,dim_sizes,dim_unlimited)
filevardef(outF,"time",typeof(time),"time")
outF->time    = time
filevardef(outF, "lat" ,typeof(lat) ,"lat")
outF->lat     = lat
filevardef(outF, "lon" ,typeof(lon) ,"lon")
outF->lon     = lon
VarAtts       = 0.0

filevardef(outF,"IVKE"    ,"float",(/"time","lat","lon"/))
VarAtts@long_name = "Integrated vapor kinetic energy"
VarAtts@units = "kg / s^2"
filevarattdef(outF, "IVKE"    ,VarAtts)
outF->IVKE     = (/IVKE(0:n_tim-1,:,:)/)

filevardef(outF,"IVT"     ,"float",(/"time","lat","lon"/))
VarAtts@long_name = "Integrated vapor transport"
VarAtts@units = "kg / m / s"
filevarattdef(outF, "IVT"     ,VarAtts)
outF->IVT      = (/IVT(0:n_tim-1,:,:)/)

; Tendency steps
outN = "./analysis_AR_tavg"+date+".nc"
outF = addfile(outN ,"c")

dim_names1     = (/"time","lev","lat","lon"/)
dim_sizes1     = (/n_tim_avg,n_lev,n_lat,n_lon/)
dim_unlimited1 = (/True,False,False,False/)
filedimdef(outF,dim_names1,dim_sizes1,dim_unlimited1)
filevardef(outF,"time",typeof(time_avg),"time")
outF->time    = time_avg
filevardef(outF, "lev" ,typeof(lev) ,"lev")
outF->lev     = lev
filevardef(outF, "lat" ,typeof(lat) ,"lat")
outF->lat     = lat
filevardef(outF, "lon" ,typeof(lon) ,"lon")
outF->lon     = lon
VarAtts       = 0.0

filevardef(outF,"VKE"     ,"float",(/"time","lev","lat","lon"/))
VarAtts@long_name = "Vapor kinetic energy"
VarAtts@units = "m^2 / s^2"
filevarattdef(outF, "VKE"     ,VarAtts)
outF->VKE      = (/VKE/)

filevardef(outF,"IVKEAVG","float",(/"time","lat","lon"/))
VarAtts@long_name = "Integrated vapor kinetic energy"
VarAtts@units = "kg / s^2"
filevarattdef(outF, "IVKEAVG",VarAtts)
outF->IVKEAVG = (/IVKEAVG/)

filevardef(outF,"IVKETEND","float",(/"time","lat","lon"/))
VarAtts@long_name = "IVKE tendency"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IVKETEND",VarAtts)
outF->IVKETEND = (/IVKETEND/)

filevardef(outF,"SPTE"    ,"float",(/"time","lat","lon"/))
VarAtts@long_name = "Surface pressure tendency effect"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "SPTE"    ,VarAtts)
outF->SPTE     = (/SPTE/)

filevardef(outF,"I_VKETEND","float",(/"time","lat","lon"/))
VarAtts@long_name = "Integral of VKE tendency"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "I_VKETEND",VarAtts)
outF->I_VKETEND = (/I_VKETEND/)

filevardef(outF,"IDOKE"   ,"float",(/"time","lat","lon"/))
VarAtts@long_name = "IVKE tendency due to dynamics on KE"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IDOKE"   ,VarAtts)
outF->IDOKE    = (/IDOKE/)

filevardef(outF,"IHAKE"   ,"float",(/"time","lat","lon"/))
VarAtts@long_name = "IVKE tendency due to horizontal advection of KE"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IHAKE"   ,VarAtts)
outF->IHAKE    = (/IHAKE/)

filevardef(outF,"IVAKE"   ,"float",(/"time","lat","lon"/))
VarAtts@long_name = "IVKE tendency due to vertical advection of KE"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IVAKE"   ,VarAtts)
outF->IVAKE    = (/IVAKE/)

filevardef(outF,"IPEKE"   ,"float",(/"time","lat","lon"/))
VarAtts@long_name = "IVKE tendency due to potential energy conversion to KE"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IPEKE"   ,VarAtts)
outF->IPEKE    = (/IPEKE/)

filevardef(outF,"IMOKE"   ,"float",(/"time","lat","lon"/))
VarAtts@long_name = "IVKE tendency due to moist convection on KE"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IMOKE"   ,VarAtts)
outF->IMOKE    = (/IMOKE/)

filevardef(outF,"ITOKE"   ,"float",(/"time","lat","lon"/))
VarAtts@long_name = "IVKE tendency due to turbulence on KE"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "ITOKE"   ,VarAtts)
outF->ITOKE    = (/ITOKE/)

filevardef(outF,"IGOKE"   ,"float",(/"time","lat","lon"/))
VarAtts@long_name = "IVKE tendency due to gravity wave drag on KE"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IGOKE"   ,VarAtts)
outF->IGOKE    = (/IGOKE/)

filevardef(outF,"IDOV"    ,"float",(/"time","lat","lon"/))
VarAtts@long_name = "IVKE tendency due to dynamics on vapor"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IDOV"    ,VarAtts)
outF->IDOV     = (/IDOV/)

filevardef(outF,"IHAV"    ,"float",(/"time","lat","lon"/))
VarAtts@long_name = "IVKE tendency due to horizontal advection of vapor"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IHAV"    ,VarAtts)
outF->IHAV     = (/IHAV/)

filevardef(outF,"IVAV"    ,"float",(/"time","lat","lon"/))
VarAtts@long_name = "IVKE tendency due to vertical advection of vapor"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IVAV"    ,VarAtts)
outF->IVAV     = (/IVAV/)

filevardef(outF,"IHCV"    ,"float",(/"time","lat","lon"/))
VarAtts@long_name = "IVKE tendency due to horizontal convergence of vapor transport"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IHCV"    ,VarAtts)
outF->IHCV     = (/IHCV/)

filevardef(outF,"IVCV"    ,"float",(/"time","lat","lon"/))
VarAtts@long_name = "IVKE tendency due to vertical convergence of vapor transport"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IVCV"    ,VarAtts)
outF->IVCV     = (/IVCV/)

filevardef(outF,"IPOV"    ,"float",(/"time","lat","lon"/))
VarAtts@long_name = "IVKE tendency due to precipitation on vapor"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IPOV"    ,VarAtts)
outF->IPOV     = (/IPOV/)

filevardef(outF,"ITOV"    ,"float",(/"time","lat","lon"/))
VarAtts@long_name = "IVKE tendency due to turbulence on vapor"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "ITOV"    ,VarAtts)
outF->ITOV     = (/ITOV/)

filevardef(outF,"ICOV"    ,"float",(/"time","lat","lon"/))
VarAtts@long_name = "IVKE tendency due to chemistry on vapor"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "ICOV"    ,VarAtts)
outF->ICOV     = (/ICOV/)

filevardef(outF,"IAOVKE"   ,"float",(/"time","lat","lon"/))
VarAtts@long_name = "IVKE tendency due to analysis on VKE"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IAOVKE"  ,VarAtts)
outF->IAOVKE   = (/IAOVKE/)

filevardef(outF,"VKETEND" ,"float",(/"time","lev","lat","lon"/))
VarAtts@long_name = "Vapor kinetic energy tendency"
VarAtts@units = "m^2 / s^3"
filevarattdef(outF, "VKETEND" ,VarAtts)
outF->VKETEND  = (/VKETEND/)

filevardef(outF,"DOKE"    ,"float",(/"time","lev","lat","lon"/))
VarAtts@long_name = "VKE tendency due to dynamics on KE"
VarAtts@units = "m^2 / s^3"
filevarattdef(outF, "DOKE"    ,VarAtts)
outF->DOKE     = (/DOKE/)

filevardef(outF,"HAKE"    ,"float",(/"time","lev","lat","lon"/))
VarAtts@long_name = "VKE tendency due to horizontal advection of KE"
VarAtts@units = "m^2 / s^3"
filevarattdef(outF, "HAKE"    ,VarAtts)
outF->HAKE     = (/HAKE/)

filevardef(outF,"VAKE"    ,"float",(/"time","lev","lat","lon"/))
VarAtts@long_name = "VKE tendency due to vertical advection of KE"
VarAtts@units = "m^2 / s^3"
filevarattdef(outF, "VAKE"    ,VarAtts)
outF->VAKE     = (/VAKE/)

filevardef(outF,"PEKE"    ,"float",(/"time","lev","lat","lon"/))
VarAtts@long_name = "VKE tendency due to potential energy conversion to KE"
VarAtts@units = "m^2 / s^3"
filevarattdef(outF, "PEKE"    ,VarAtts)
outF->PEKE     = (/PEKE/)

filevardef(outF,"MOKE"    ,"float",(/"time","lev","lat","lon"/))
VarAtts@long_name = "VKE tendency due to moist convection on KE"
VarAtts@units = "m^2 / s^3"
filevarattdef(outF, "MOKE"    ,VarAtts)
outF->MOKE     = (/MOKE/)

filevardef(outF,"TOKE"    ,"float",(/"time","lev","lat","lon"/))
VarAtts@long_name = "VKE tendency due to turbulence on KE"
VarAtts@units = "m^2 / s^3"
filevarattdef(outF, "TOKE"    ,VarAtts)
outF->TOKE     = (/TOKE/)

filevardef(outF,"GOKE"    ,"float",(/"time","lev","lat","lon"/))
VarAtts@long_name = "VKE tendency due to gravity wave drag on KE"
VarAtts@units = "m^2 / s^3"
filevarattdef(outF, "GOKE"    ,VarAtts)
outF->GOKE     = (/GOKE/)

filevardef(outF,"DOV"     ,"float",(/"time","lev","lat","lon"/))
VarAtts@long_name = "VKE tendency due to dynamics on vapor"
VarAtts@units = "m^2 / s^3"
filevarattdef(outF, "DOV"     ,VarAtts)
outF->DOV      = (/DOV/)

filevardef(outF,"HAV"     ,"float",(/"time","lev","lat","lon"/))
VarAtts@long_name = "VKE tendency due to horizontal advection of vapor"
VarAtts@units = "m^2 / s^3"
filevarattdef(outF, "HAV"     ,VarAtts)
outF->HAV      = (/HAV/)

filevardef(outF,"VAV"     ,"float",(/"time","lev","lat","lon"/))
VarAtts@long_name = "VKE tendency due to vertical advection of vapor"
VarAtts@units = "m^2 / s^3"
filevarattdef(outF, "VAV"     ,VarAtts)
outF->VAV      = (/VAV/)

filevardef(outF,"HCV"     ,"float",(/"time","lev","lat","lon"/))
VarAtts@long_name = "VKE tendency due to horizontal convergence of vapor transport"
VarAtts@units = "m^2 / s^3"
filevarattdef(outF, "HCV"     ,VarAtts)
outF->HCV      = (/HCV/)

filevardef(outF,"VCV"     ,"float",(/"time","lev","lat","lon"/))
VarAtts@long_name = "VKE tendency due to vertical convergence of vapor transport"
VarAtts@units = "m^2 / s^3"
filevarattdef(outF, "VCV"     ,VarAtts)
outF->VCV      = (/VCV/)

filevardef(outF,"POV"     ,"float",(/"time","lev","lat","lon"/))
VarAtts@long_name = "VKE tendency due to precipitation on vapor"
VarAtts@units = "m^2 / s^3"
filevarattdef(outF, "POV"     ,VarAtts)
outF->POV      = (/POV/)

filevardef(outF,"TOV"     ,"float",(/"time","lev","lat","lon"/))
VarAtts@long_name = "VKE tendency due to turbulence on vapor"
VarAtts@units = "m^2 / s^3"
filevarattdef(outF, "TOV"     ,VarAtts)
outF->TOV      = (/TOV/)

filevardef(outF,"COV"     ,"float",(/"time","lev","lat","lon"/))
VarAtts@long_name = "VKE tendency due to chemistry on vapor"
VarAtts@units = "m^2 / s^3"
filevarattdef(outF, "COV"     ,VarAtts)
outF->COV      = (/COV/)

filevardef(outF,"AOVKE"   ,"float",(/"time","lev","lat","lon"/))
VarAtts@long_name = "VKE tendency due to analysis"
VarAtts@units = "m^2 / s^3"
filevarattdef(outF, "AOVKE"   ,VarAtts)
outF->AOVKE    = (/AOVKE/)

print("Writing ends")

