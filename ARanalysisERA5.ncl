; Numerical analyses of IVT, VKE, IVKE, and tendency terms from ERA5 data
; Here is a template for executing this script in a PBS job array:
;
;#PBS -J 20190101-20191231
;MONTH=$((PBS_ARRAY_INDEX/100))
;PLDIR=/glade/collections/rda/data/ds633.0/e5.oper.an.pl/$MONTH/
;SFCFN=`ls /glade/collections/rda/data/ds633.0/e5.oper.an.sfc/$MONTH/e5.oper.an.sfc.128_134_sp.ll025sc.*.nc`
;### Run program
;ncl 'date='$PBS_ARRAY_INDEX
;    'f_q=addfile("'$PLDIR'e5.oper.an.pl.128_133_q.ll025sc.'$PBS_ARRAY_INDEX'00_'$PBS_ARRAY_INDEX'23.nc","r")'
;    'f_u=addfile("'$PLDIR'e5.oper.an.pl.128_131_u.ll025uv.'$PBS_ARRAY_INDEX'00_'$PBS_ARRAY_INDEX'23.nc","r")'
;    'f_v=addfile("'$PLDIR'e5.oper.an.pl.128_132_v.ll025uv.'$PBS_ARRAY_INDEX'00_'$PBS_ARRAY_INDEX'23.nc","r")'
;    'f_w=addfile("'$PLDIR'e5.oper.an.pl.128_135_w.ll025sc.'$PBS_ARRAY_INDEX'00_'$PBS_ARRAY_INDEX'23.nc","r")'
;    'f_z=addfile("'$PLDIR'e5.oper.an.pl.128_129_z.ll025sc.'$PBS_ARRAY_INDEX'00_'$PBS_ARRAY_INDEX'23.nc","r")'
;    'f_sp=addfile("'$SFCFN'","r")' ARanalysisERA5.ncl
;
; Hing Ong, 2 Oct 2024

; Parameter setting
g   =    9.81
pB  = 1000.
pT  =  200.
pTp =  150.

; Add input files
if (ismissing(f_q)) then
    exit
end if

print("Reading begins")
time     = f_q->time
n_tim    = dimsizes(time)
lev      = f_q->level   (          {pB:pTp}    )
n_lev    = dimsizes(lev)
lat      = f_q->latitude
n_lat    = dimsizes(lat)
lon      = f_q->longitude
n_lon    = dimsizes(lon)
ps       = f_sp->SP     ({time},            :,:)
Q        = f_q->Q       (0:n_tim-1,{pB:pTp},:,:)
U        = f_u->U       (0:n_tim-1,{pB:pTp},:,:)
V        = f_v->V       (0:n_tim-1,{pB:pTp},:,:)
H        = f_z->Z       (0:n_tim-1,{pB:pTp},:,:)
OMEGA    = f_w->W       (0:n_tim-1,{pB:pTp},:,:)
print("Reading ends")

print("Calculation begins")

; Instantaneous steps
p        = dble2flt(lev)*100.
pB       = pB*100.
pT       = pT*100.
LAT      = dble2flt(lat)
LON      = dble2flt(lon)
QU       = Q ; copy attributes
QV       = Q ; copy attributes
VKE      = Q ; copy attributes
QU       = Q*U
QV       = Q*V
K        = (U*U+V*V)*0.5
VKE      = Q*Q*K
IQU      = vibeta(p,QU (time|:,latitude|:,longitude|:,level|::-1),1,ps,pB,pT) /g
IQV      = vibeta(p,QV (time|:,latitude|:,longitude|:,level|::-1),1,ps,pB,pT) /g
IVKE     = vibeta(p,VKE(time|:,latitude|:,longitude|:,level|::-1),1,ps,pB,pT) /g
delete(VKE)
IVT      = sqrt(IQU*IQU+IQV*IQV)
Q2U      = Q*QU
Q2V      = Q*QV
KQ2      = U*QU+V*QV

dUdp     = center_finite_diff_n(U,p,False,0,1)
dVdp     = center_finite_diff_n(V,p,False,0,1)
VAKE     = Q ; copy attributes
VAKE     = -Q2U*OMEGA*dUdp-Q2V*OMEGA*dVdp
delete(dUdp)
delete(dVdp)
IVAKE    = vibeta(p,VAKE (time|:,latitude|:,longitude|:,level|::-1),1,ps,pB,pT) /g
delete(VAKE)

dQdp     = center_finite_diff_n(Q,p,False,0,1)
VAV      = Q ; copy attributes
VAV      = -KQ2*OMEGA*dQdp
delete(dQdp)
IVAV     = vibeta(p,VAV  (time|:,latitude|:,longitude|:,level|:),1,ps,pB,pT) /g
delete(VAV)

IHAKE    = IVAKE
IHAKE    = IVAKE@_FillValue
IPEKE    = IHAKE
IHAV     = IHAKE

gradK    = grad_latlon_cfd(K  ,LAT,LON,True,False)
HAKE     = Q ; copy attributes
HAKE     = -  Q2U*gradK[1]  -  Q2V*gradK[0]
delete(gradK)
IHAKE(:,1:n_lat-2,:) = vibeta(p,HAKE (time|:,latitude|1:n_lat-2,longitude|:,level|::-1),1,ps(:,1:n_lat-2,:),pB,pT) /g
delete(HAKE)

gradH    = grad_latlon_cfd(H  ,LAT,LON,True,False)
PEKE     = Q ; copy attributes
PEKE     = -  Q2U*gradH[1]  -  Q2V*gradH[0]
delete(gradH)
IPEKE(:,1:n_lat-2,:) = vibeta(p,PEKE (time|:,latitude|1:n_lat-2,longitude|:,level|::-1),1,ps(:,1:n_lat-2,:),pB,pT) /g
delete(PEKE)

gradQ    = grad_latlon_cfd(Q  ,LAT,LON,True,False)
HAV      = Q ; copy attributes
HAV      = -KQ2*U*gradQ[1]  -KQ2*V*gradQ[0]
delete(gradQ)
IHAV(:,1:n_lat-2,:)  = vibeta(p,HAV  (time|:,latitude|1:n_lat-2,longitude|:,level|::-1),1,ps(:,1:n_lat-2,:),pB,pT) /g
delete(HAV)

KEpart   = Q(1:n_tim,:,:,:) ; copy attributes
Q2U      = Q(1:n_tim,:,:,:)*Q(1:n_tim,:,:,:)*U(1:n_tim,:,:,:)
Q2V      = Q(1:n_tim,:,:,:)*Q(1:n_tim,:,:,:)*V(1:n_tim,:,:,:)
Utend    = (U(2:n_tim+1,:,:,:)-U(0:n_tim-1,:,:,:))/7200.
Vtend    = (V(2:n_tim+1,:,:,:)-V(0:n_tim-1,:,:,:))/7200.
KEpart   = Q2U*Utend+Q2V*Vtend
delete(Utend)
delete(Vtend)
IKEpart  = vibeta(p,KEpart(time|:,latitude|:,longitude|:,level|::-1),1,ps,pB,pT) /g
delete(KEpart)

VPpart   = Q(1:n_tim,:,:,:) ; copy attributes
KQ2      = Q(1:n_tim,:,:,:)*(U(1:n_tim,:,:,:)*U(1:n_tim,:,:,:)+V(1:n_tim,:,:,:)*V(1:n_tim,:,:,:))
Qtend    = (Q(2:n_tim+1,:,:,:)-Q(0:n_tim-1,:,:,:))/7200.
VPpart   = KQ2*Qtend
delete(Qtend)
IVPpart  = vibeta(p,VPpart(time|:,latitude|:,longitude|:,level|::-1),1,ps,pB,pT) /g
delete(VPpart)

print("Calculation ends")

print("Writing begins")

outN = "./analysis_AR_era5_"+date+".nc"
outF = addfile(outN ,"c")

; Instantaneous steps

dim_names     = (/"time","latitude","longitude"/)
dim_sizes     = (/n_tim,n_lat,n_lon/)
dim_unlimited = (/True,False,False/)
filedimdef(outF,dim_names,dim_sizes,dim_unlimited)
filevardef(outF,"time",typeof(time),"time")
outF->time    = time
;filevardef(outF, "level" ,typeof(lev) ,"level")
;outF->level   = lev
filevardef(outF, "latitude" ,typeof(lat) ,"latitude")
outF->latitude= lat
filevardef(outF, "longitude" ,typeof(lon) ,"longitude")
outF->longitude=lon
VarAtts       = 0.0

filevardef(outF,"IVKE"    ,"float",(/"time","latitude","longitude"/))
VarAtts@long_name = "Integrated vapor kinetic energy"
VarAtts@units = "kg / s^2"
filevarattdef(outF, "IVKE"    ,VarAtts)
outF->IVKE     = (/IVKE/)

filevardef(outF,"IVT"     ,"float",(/"time","latitude","longitude"/))
VarAtts@long_name = "Integrated vapor transport"
VarAtts@units = "kg / m / s"
filevarattdef(outF, "IVT"     ,VarAtts)
outF->IVT      = (/IVT/)

;filevardef(outF,"VKE"     ,"float",(/"time","level","latitude","longitude"/))
;VarAtts@long_name = "Vapor kinetic energy"
;VarAtts@units = "m^2 / s^2"
;filevarattdef(outF, "VKE"     ,VarAtts)
;outF->VKE      = (/VKE/)

filevardef(outF,"IHAKE"   ,"float",(/"time","latitude","longitude"/))
VarAtts@long_name = "IVKE tendency due to horizontal advection of KE"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IHAKE"   ,VarAtts)
outF->IHAKE    = (/IHAKE/)

filevardef(outF,"IVAKE"   ,"float",(/"time","latitude","longitude"/))
VarAtts@long_name = "IVKE tendency due to vertical advection of KE"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IVAKE"   ,VarAtts)
outF->IVAKE    = (/IVAKE/)

filevardef(outF,"IPEKE"   ,"float",(/"time","latitude","longitude"/))
VarAtts@long_name = "IVKE tendency due to potential energy conversion to KE"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IPEKE"   ,VarAtts)
outF->IPEKE    = (/IPEKE/)

filevardef(outF,"IHAV"    ,"float",(/"time","latitude","longitude"/))
VarAtts@long_name = "IVKE tendency due to horizontal advection of vapor"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IHAV"    ,VarAtts)
outF->IHAV     = (/IHAV/)

filevardef(outF,"IVAV"    ,"float",(/"time","latitude","longitude"/))
VarAtts@long_name = "IVKE tendency due to vertical advection of vapor"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IVAV"    ,VarAtts)
outF->IVAV     = (/IVAV/)

filevardef(outF,"IKEpart","float",(/"time","latitude","longitude"/))
VarAtts@long_name = "KE part of Integrated VKE tendency"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IKEpart",VarAtts)
outF->IKEpart = (/IKEpart/)

filevardef(outF,"IVPpart","float",(/"time","latitude","longitude"/))
VarAtts@long_name = "Vapor part of Integrated VKE tendency"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IVPpart",VarAtts)
outF->IVPpart = (/IVPpart/)

;filevardef(outF,"HAKE"    ,"float",(/"time","level","latitude","longitude"/))
;VarAtts@long_name = "VKE tendency due to horizontal advection of KE"
;VarAtts@units = "m^2 / s^3"
;filevarattdef(outF, "HAKE"    ,VarAtts)
;outF->HAKE     = (/HAKE/)

;filevardef(outF,"VAKE"    ,"float",(/"time","level","latitude","longitude"/))
;VarAtts@long_name = "VKE tendency due to vertical advection of KE"
;VarAtts@units = "m^2 / s^3"
;filevarattdef(outF, "VAKE"    ,VarAtts)
;outF->VAKE     = (/VAKE/)

;filevardef(outF,"PEKE"    ,"float",(/"time","level","latitude","longitude"/))
;VarAtts@long_name = "VKE tendency due to potential energy conversion to KE"
;VarAtts@units = "m^2 / s^3"
;filevarattdef(outF, "PEKE"    ,VarAtts)
;outF->PEKE     = (/PEKE/)

;filevardef(outF,"HAV"     ,"float",(/"time","level","latitude","longitude"/))
;VarAtts@long_name = "VKE tendency due to horizontal advection of vapor"
;VarAtts@units = "m^2 / s^3"
;filevarattdef(outF, "HAV"     ,VarAtts)
;outF->HAV      = (/HAV/)

;filevardef(outF,"VAV"     ,"float",(/"time","level","latitude","longitude"/))
;VarAtts@long_name = "VKE tendency due to vertical advection of vapor"
;VarAtts@units = "m^2 / s^3"
;filevarattdef(outF, "VAV"     ,VarAtts)
;outF->VAV      = (/VAV/)

print("Writing ends")

