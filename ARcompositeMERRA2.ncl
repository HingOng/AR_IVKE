; AR composite using analysis data from ARanalysisMERRA2.ncl
; Hing Ong, 2 Oct 2024

; Parameter setting
siglvl =   0.05
lag    =   0
latW   =  34.5
latE   =  35.5
lonS   = 154.5
lonN   = 155.5

; Add input files
fn = systemfunc("find $SCRATCH/MERRA2_201*/* -name 'analysis_AR_inst*' | sort -k11")
f  = addfiles(fn,"r")
fn = systemfunc("find $SCRATCH/MERRA2_201*/* -name 'analysis_AR_tavg*' | sort -k11")
f1 = addfiles(fn,"r")

; Reading predictor
IVKE       = f[:]->IVKE(:,{20:50},{130:180})
lat        = f[0]->lat(   {20:50}          )
lon        = f[0]->lon(           {130:180})
dsize      = dimsizes(IVKE)
tsize      = dsize(0)
latsize    = dsize(1)
lonsize    = dsize(2)

; Masking
;AR1        = AR_(0:tsize-2,:,:)+AR_(1:tsize-1,:,:)
;AR         = runave(dim_avg(dim_avg(AR1)),17,1)
predictor_ = dim_avg(dim_avg(IVKE(:,{latW:latE},{lonS:lonN})))
predictor  = (predictor_(0:tsize-2)+predictor_(1:tsize-1))*0.5
;predictor  = mask(predictor,AR.eq.0,False)
predictor0 = stddev(predictor)

outN = "./regression_AR.nc"
outF = addfile(outN ,"c")

dim_names1     = (/"lat","lon"/)
dim_sizes1     = (/latsize,lonsize/)
dim_unlimited1 = (/False,False/)
filedimdef(outF,dim_names1,dim_sizes1,dim_unlimited1)
filevardef(outF, "lat" ,typeof(lat) ,"lat")
outF->lat     = lat
filevardef(outF, "lon" ,typeof(lon) ,"lon")
outF->lon     = lon
VarAtts       = 0.0

; To do: Write a function to make the following scripts more concise.

; Reading predictand
predictand = f1[:]->IVKETEND(0:tsize-2,{20:50},{130:180})
; Regressing
RegCoef    = regCoef_n(predictor,predictand,0,0)
TVal       = onedtond(RegCoef@tval,dimsizes(RegCoef))
DegFree    = equiv_sample_size(predictand(lat|:,lon|:,time|:),siglvl,0)-2
B          = TVal
B          = 0.5
Probability= betainc(DegFree/(DegFree+TVal^2),DegFree*0.5,B)
Composite  = RegCoef*predictor0
filevardef(outF,"IVKETEND","float",(/"lat","lon"/))
VarAtts@long_name = "IVKE tendency (regressed)"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IVKETEND",VarAtts)
outF->IVKETEND = (/Composite/)
filevardef(outF,"IVKETEND_prob","float",(/"lat","lon"/))
VarAtts@long_name = "Probability of IVKETEND (regressed)"
VarAtts@units = ""
filevarattdef(outF, "IVKETEND_prob",VarAtts)
outF->IVKETEND_prob = (/Probability/)

; Reading predictand
ivke       = predictand
ivke       = (IVKE(0:tsize-2,:,:)+IVKE(1:tsize-1,:,:))*0.5
predictand = ivke
; Regressing
RegCoef    = regCoef_n(predictor,predictand,0,0)
TVal       = onedtond(RegCoef@tval,dimsizes(RegCoef))
DegFree    = equiv_sample_size(predictand(lat|:,lon|:,time|:),siglvl,0)-2
B          = TVal
B          = 0.5
Probability= betainc(DegFree/(DegFree+TVal^2),DegFree*0.5,B)
Composite  = RegCoef*predictor0
filevardef(outF,"IVKE"     ,"float",(/"lat","lon"/))
VarAtts@long_name = "Integrated vapor kinetic energy (regressed)"
VarAtts@units = "kg / s^2"
filevarattdef(outF, "IVKE"     ,VarAtts)
outF->IVKE      = (/Composite/)
filevardef(outF,"IVKE_prob","float",(/"lat","lon"/))
VarAtts@long_name = "Probability of IVKE (regressed)"
VarAtts@units = ""
filevarattdef(outF, "IVKE_prob"     ,VarAtts)
outF->IVKE_prob = (/Probability/)

; Reading predictand
predictand = f1[:]->SPTE(0:tsize-2,{20:50},{130:180})
; Regressing
RegCoef    = regCoef_n(predictor,predictand,0,0)
TVal       = onedtond(RegCoef@tval,dimsizes(RegCoef))
DegFree    = equiv_sample_size(predictand(lat|:,lon|:,time|:),siglvl,0)-2
B          = TVal
B          = 0.5
Probability= betainc(DegFree/(DegFree+TVal^2),DegFree*0.5,B)
Composite  = RegCoef*predictor0
filevardef(outF,"SPTE","float",(/"lat","lon"/))
VarAtts@long_name = "Surface pressure tendency effect (regressed)"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "SPTE",VarAtts)
outF->SPTE = (/Composite/)
filevardef(outF,"SPTE_prob","float",(/"lat","lon"/))
VarAtts@long_name = "Probability of SPTE (regressed)"
VarAtts@units = ""
filevarattdef(outF, "SPTE_prob",VarAtts)
outF->SPTE_prob = (/Probability/)

; Reading predictand
predictand = f1[:]->I_VKETEND(0:tsize-2,{20:50},{130:180})
; Regressing
RegCoef    = regCoef_n(predictor,predictand,0,0)
TVal       = onedtond(RegCoef@tval,dimsizes(RegCoef))
DegFree    = equiv_sample_size(predictand(lat|:,lon|:,time|:),siglvl,0)-2
B          = TVal
B          = 0.5
Probability= betainc(DegFree/(DegFree+TVal^2),DegFree*0.5,B)
Composite  = RegCoef*predictor0
filevardef(outF,"I_VKETEND","float",(/"lat","lon"/))
VarAtts@long_name = "Integral of VKE tendency (regressed)"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "I_VKETEND",VarAtts)
outF->I_VKETEND = (/Composite/)
filevardef(outF,"I_VKETEND_prob","float",(/"lat","lon"/))
VarAtts@long_name = "Probability of I_VKETEND (regressed)"
VarAtts@units = ""
filevarattdef(outF, "I_VKETEND_prob",VarAtts)
outF->I_VKETEND_prob = (/Probability/)

; Reading predictand
predictand = f1[:]->IDOKE(0:tsize-2,{20:50},{130:180})
; Regressing
RegCoef    = regCoef_n(predictor,predictand,0,0)
TVal       = onedtond(RegCoef@tval,dimsizes(RegCoef))
DegFree    = equiv_sample_size(predictand(lat|:,lon|:,time|:),siglvl,0)-2
B          = TVal
B          = 0.5
Probability= betainc(DegFree/(DegFree+TVal^2),DegFree*0.5,B)
Composite  = RegCoef*predictor0
filevardef(outF,"IDOKE","float",(/"lat","lon"/))
VarAtts@long_name = "IVKE tendency due to dynamics on KE (regressed)"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IDOKE",VarAtts)
outF->IDOKE = (/Composite/)
filevardef(outF,"IDOKE_prob","float",(/"lat","lon"/))
VarAtts@long_name = "Probability of IDOKE (regressed)"
VarAtts@units = ""
filevarattdef(outF, "IDOKE_prob",VarAtts)
outF->IDOKE_prob = (/Probability/)

; Reading predictand
predictand = f1[:]->IDOV(0:tsize-2,{20:50},{130:180})
; Regressing
RegCoef    = regCoef_n(predictor,predictand,0,0)
TVal       = onedtond(RegCoef@tval,dimsizes(RegCoef))
DegFree    = equiv_sample_size(predictand(lat|:,lon|:,time|:),siglvl,0)-2
B          = TVal
B          = 0.5
Probability= betainc(DegFree/(DegFree+TVal^2),DegFree*0.5,B)
Composite  = RegCoef*predictor0
filevardef(outF,"IDOV","float",(/"lat","lon"/))
VarAtts@long_name = "IVKE tendency due to dynamics on vapor (regressed)"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IDOV",VarAtts)
outF->IDOV = (/Composite/)
filevardef(outF,"IDOV_prob","float",(/"lat","lon"/))
VarAtts@long_name = "Probability of IDOV (regressed)"
VarAtts@units = ""
filevarattdef(outF, "IDOV_prob",VarAtts)
outF->IDOV_prob = (/Probability/)

; Reading predictand
predictand = f1[:]->IHAKE(0:tsize-2,{20:50},{130:180})
; Regressing
RegCoef    = regCoef_n(predictor,predictand,0,0)
TVal       = onedtond(RegCoef@tval,dimsizes(RegCoef))
DegFree    = equiv_sample_size(predictand(lat|:,lon|:,time|:),siglvl,0)-2
B          = TVal
B          = 0.5
Probability= betainc(DegFree/(DegFree+TVal^2),DegFree*0.5,B)
Composite  = RegCoef*predictor0
filevardef(outF,"IHAKE","float",(/"lat","lon"/))
VarAtts@long_name = "IVKE tendency due to horizontal advection of KE (regressed)"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IHAKE",VarAtts)
outF->IHAKE = (/Composite/)
filevardef(outF,"IHAKE_prob","float",(/"lat","lon"/))
VarAtts@long_name = "Probability of IHAKE (regressed)"
VarAtts@units = ""
filevarattdef(outF, "IHAKE_prob",VarAtts)
outF->IHAKE_prob = (/Probability/)

; Reading predictand
predictand = f1[:]->IHAV(0:tsize-2,{20:50},{130:180})
; Regressing
RegCoef    = regCoef_n(predictor,predictand,0,0)
TVal       = onedtond(RegCoef@tval,dimsizes(RegCoef))
DegFree    = equiv_sample_size(predictand(lat|:,lon|:,time|:),siglvl,0)-2
B          = TVal
B          = 0.5
Probability= betainc(DegFree/(DegFree+TVal^2),DegFree*0.5,B)
Composite  = RegCoef*predictor0
filevardef(outF,"IHAV","float",(/"lat","lon"/))
VarAtts@long_name = "IVKE tendency due to horizontal advection of vapor (regressed)"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IHAV",VarAtts)
outF->IHAV = (/Composite/)
filevardef(outF,"IHAV_prob","float",(/"lat","lon"/))
VarAtts@long_name = "Probability of IHAV (regressed)"
VarAtts@units = ""
filevarattdef(outF, "IHAV_prob",VarAtts)
outF->IHAV_prob = (/Probability/)

; Reading predictand
predictand = f1[:]->IVAKE(0:tsize-2,{20:50},{130:180})
; Regressing
RegCoef    = regCoef_n(predictor,predictand,0,0)
TVal       = onedtond(RegCoef@tval,dimsizes(RegCoef))
DegFree    = equiv_sample_size(predictand(lat|:,lon|:,time|:),siglvl,0)-2
B          = TVal
B          = 0.5
Probability= betainc(DegFree/(DegFree+TVal^2),DegFree*0.5,B)
Composite  = RegCoef*predictor0
filevardef(outF,"IVAKE","float",(/"lat","lon"/))
VarAtts@long_name = "IVKE tendency due to vertical advection of KE (regressed)"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IVAKE",VarAtts)
outF->IVAKE = (/Composite/)
filevardef(outF,"IVAKE_prob","float",(/"lat","lon"/))
VarAtts@long_name = "Probability of IVAKE (regressed)"
VarAtts@units = ""
filevarattdef(outF, "IVAKE_prob",VarAtts)
outF->IVAKE_prob = (/Probability/)

; Reading predictand
predictand = f1[:]->IVAV(0:tsize-2,{20:50},{130:180})
; Regressing
RegCoef    = regCoef_n(predictor,predictand,0,0)
TVal       = onedtond(RegCoef@tval,dimsizes(RegCoef))
DegFree    = equiv_sample_size(predictand(lat|:,lon|:,time|:),siglvl,0)-2
B          = TVal
B          = 0.5
Probability= betainc(DegFree/(DegFree+TVal^2),DegFree*0.5,B)
Composite  = RegCoef*predictor0
filevardef(outF,"IVAV","float",(/"lat","lon"/))
VarAtts@long_name = "IVKE tendency due to vertical advection of vapor (regressed)"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IVAV",VarAtts)
outF->IVAV = (/Composite/)
filevardef(outF,"IVAV_prob","float",(/"lat","lon"/))
VarAtts@long_name = "Probability of IVAV (regressed)"
VarAtts@units = ""
filevarattdef(outF, "IVAV_prob",VarAtts)
outF->IVAV_prob = (/Probability/)

; Reading predictand
predictand = f1[:]->IHCV(0:tsize-2,{20:50},{130:180})
; Regressing
RegCoef    = regCoef_n(predictor,predictand,0,0)
TVal       = onedtond(RegCoef@tval,dimsizes(RegCoef))
DegFree    = equiv_sample_size(predictand(lat|:,lon|:,time|:),siglvl,0)-2
B          = TVal
B          = 0.5
Probability= betainc(DegFree/(DegFree+TVal^2),DegFree*0.5,B)
Composite  = RegCoef*predictor0
filevardef(outF,"IHCV","float",(/"lat","lon"/))
VarAtts@long_name = "IVKE tendency due to horizontal convergence of vapor flux (regressed)"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IHCV",VarAtts)
outF->IHCV = (/Composite/)
filevardef(outF,"IHCV_prob","float",(/"lat","lon"/))
VarAtts@long_name = "Probability of IHCV (regressed)"
VarAtts@units = ""
filevarattdef(outF, "IHCV_prob",VarAtts)
outF->IHCV_prob = (/Probability/)

; Reading predictand
predictand = f1[:]->IVCV(0:tsize-2,{20:50},{130:180})
; Regressing
RegCoef    = regCoef_n(predictor,predictand,0,0)
TVal       = onedtond(RegCoef@tval,dimsizes(RegCoef))
DegFree    = equiv_sample_size(predictand(lat|:,lon|:,time|:),siglvl,0)-2
B          = TVal
B          = 0.5
Probability= betainc(DegFree/(DegFree+TVal^2),DegFree*0.5,B)
Composite  = RegCoef*predictor0
filevardef(outF,"IVCV","float",(/"lat","lon"/))
VarAtts@long_name = "IVKE tendency due to vertical convergence of vapor flux (regressed)"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IVCV",VarAtts)
outF->IVCV = (/Composite/)
filevardef(outF,"IVCV_prob","float",(/"lat","lon"/))
VarAtts@long_name = "Probability of IVCV (regressed)"
VarAtts@units = ""
filevarattdef(outF, "IVCV_prob",VarAtts)
outF->IVCV_prob = (/Probability/)

; Reading predictand
predictand = f1[:]->IPEKE(0:tsize-2,{20:50},{130:180})
; Regressing
RegCoef    = regCoef_n(predictor,predictand,0,0)
TVal       = onedtond(RegCoef@tval,dimsizes(RegCoef))
DegFree    = equiv_sample_size(predictand(lat|:,lon|:,time|:),siglvl,0)-2
B          = TVal
B          = 0.5
Probability= betainc(DegFree/(DegFree+TVal^2),DegFree*0.5,B)
Composite  = RegCoef*predictor0
filevardef(outF,"IPEKE","float",(/"lat","lon"/))
VarAtts@long_name = "IVKE tendency due to potential energy conversion to KE (regressed)"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IPEKE",VarAtts)
outF->IPEKE = (/Composite/)
filevardef(outF,"IPEKE_prob","float",(/"lat","lon"/))
VarAtts@long_name = "Probability of IPEKE (regressed)"
VarAtts@units = ""
filevarattdef(outF, "IPEKE_prob",VarAtts)
outF->IPEKE_prob = (/Probability/)

; Reading predictand
predictand = f1[:]->IMOKE(0:tsize-2,{20:50},{130:180})
; Regressing
RegCoef    = regCoef_n(predictor,predictand,0,0)
TVal       = onedtond(RegCoef@tval,dimsizes(RegCoef))
DegFree    = equiv_sample_size(predictand(lat|:,lon|:,time|:),siglvl,0)-2
B          = TVal
B          = 0.5
Probability= betainc(DegFree/(DegFree+TVal^2),DegFree*0.5,B)
Composite  = RegCoef*predictor0
filevardef(outF,"IMOKE","float",(/"lat","lon"/))
VarAtts@long_name = "IVKE tendency due to moist convection on KE (regressed)"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IMOKE",VarAtts)
outF->IMOKE = (/Composite/)
filevardef(outF,"IMOKE_prob","float",(/"lat","lon"/))
VarAtts@long_name = "Probability of IMOKE (regressed)"
VarAtts@units = ""
filevarattdef(outF, "IMOKE_prob",VarAtts)
outF->IMOKE_prob = (/Probability/)

; Reading predictand
predictand = f1[:]->IPOV(0:tsize-2,{20:50},{130:180})
; Regressing
RegCoef    = regCoef_n(predictor,predictand,0,0)
TVal       = onedtond(RegCoef@tval,dimsizes(RegCoef))
DegFree    = equiv_sample_size(predictand(lat|:,lon|:,time|:),siglvl,0)-2
B          = TVal
B          = 0.5
Probability= betainc(DegFree/(DegFree+TVal^2),DegFree*0.5,B)
Composite  = RegCoef*predictor0
filevardef(outF,"IPOV","float",(/"lat","lon"/))
VarAtts@long_name = "IVKE tendency due to precipitation on vapor (regressed)"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IPOV",VarAtts)
outF->IPOV = (/Composite/)
filevardef(outF,"IPOV_prob","float",(/"lat","lon"/))
VarAtts@long_name = "Probability of IPOV (regressed)"
VarAtts@units = ""
filevarattdef(outF, "IPOV_prob",VarAtts)
outF->IPOV_prob = (/Probability/)

; Reading predictand
predictand = f1[:]->ITOKE(0:tsize-2,{20:50},{130:180})
; Regressing
RegCoef    = regCoef_n(predictor,predictand,0,0)
TVal       = onedtond(RegCoef@tval,dimsizes(RegCoef))
DegFree    = equiv_sample_size(predictand(lat|:,lon|:,time|:),siglvl,0)-2
B          = TVal
B          = 0.5
Probability= betainc(DegFree/(DegFree+TVal^2),DegFree*0.5,B)
Composite  = RegCoef*predictor0
filevardef(outF,"ITOKE","float",(/"lat","lon"/))
VarAtts@long_name = "IVKE tendency due to turbulence on KE (regressed)"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "ITOKE",VarAtts)
outF->ITOKE = (/Composite/)
filevardef(outF,"ITOKE_prob","float",(/"lat","lon"/))
VarAtts@long_name = "Probability of ITOKE (regressed)"
VarAtts@units = ""
filevarattdef(outF, "ITOKE_prob",VarAtts)
outF->ITOKE_prob = (/Probability/)

; Reading predictand
predictand = f1[:]->ITOV(0:tsize-2,{20:50},{130:180})
; Regressing
RegCoef    = regCoef_n(predictor,predictand,0,0)
TVal       = onedtond(RegCoef@tval,dimsizes(RegCoef))
DegFree    = equiv_sample_size(predictand(lat|:,lon|:,time|:),siglvl,0)-2
B          = TVal
B          = 0.5
Probability= betainc(DegFree/(DegFree+TVal^2),DegFree*0.5,B)
Composite  = RegCoef*predictor0
filevardef(outF,"ITOV","float",(/"lat","lon"/))
VarAtts@long_name = "IVKE tendency due to turbulence on vapor (regressed)"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "ITOV",VarAtts)
outF->ITOV = (/Composite/)
filevardef(outF,"ITOV_prob","float",(/"lat","lon"/))
VarAtts@long_name = "Probability of ITOV (regressed)"
VarAtts@units = ""
filevarattdef(outF, "ITOV_prob",VarAtts)
outF->ITOV_prob = (/Probability/)

; Reading predictand
predictand = f1[:]->IGOKE(0:tsize-2,{20:50},{130:180})
; Regressing
RegCoef    = regCoef_n(predictor,predictand,0,0)
TVal       = onedtond(RegCoef@tval,dimsizes(RegCoef))
DegFree    = equiv_sample_size(predictand(lat|:,lon|:,time|:),siglvl,0)-2
B          = TVal
B          = 0.5
Probability= betainc(DegFree/(DegFree+TVal^2),DegFree*0.5,B)
Composite  = RegCoef*predictor0
filevardef(outF,"IGOKE","float",(/"lat","lon"/))
VarAtts@long_name = "IVKE tendency due to gravity wave drag on KE (regressed)"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IGOKE",VarAtts)
outF->IGOKE = (/Composite/)
filevardef(outF,"IGOKE_prob","float",(/"lat","lon"/))
VarAtts@long_name = "Probability of IGOKE (regressed)"
VarAtts@units = ""
filevarattdef(outF, "IGOKE_prob",VarAtts)
outF->IGOKE_prob = (/Probability/)

; Reading predictand
predictand = f1[:]->ICOV(0:tsize-2,{20:50},{130:180})
; Regressing
RegCoef    = regCoef_n(predictor,predictand,0,0)
TVal       = onedtond(RegCoef@tval,dimsizes(RegCoef))
DegFree    = equiv_sample_size(predictand(lat|:,lon|:,time|:),siglvl,0)-2
B          = TVal
B          = 0.5
Probability= betainc(DegFree/(DegFree+TVal^2),DegFree*0.5,B)
Composite  = RegCoef*predictor0
filevardef(outF,"ICOV","float",(/"lat","lon"/))
VarAtts@long_name = "IVKE tendency due to chemistry on vapor (regressed)"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "ICOV",VarAtts)
outF->ICOV = (/Composite/)
filevardef(outF,"ICOV_prob","float",(/"lat","lon"/))
VarAtts@long_name = "Probability of ICOV (regressed)"
VarAtts@units = ""
filevarattdef(outF, "ICOV_prob",VarAtts)
outF->ICOV_prob = (/Probability/)

; Reading predictand
predictand = f1[:]->IAOVKE(0:tsize-2,{20:50},{130:180})
; Regressing
RegCoef    = regCoef_n(predictor,predictand,0,0)
TVal       = onedtond(RegCoef@tval,dimsizes(RegCoef))
DegFree    = equiv_sample_size(predictand(lat|:,lon|:,time|:),siglvl,0)-2
B          = TVal
B          = 0.5
Probability= betainc(DegFree/(DegFree+TVal^2),DegFree*0.5,B)
Composite  = RegCoef*predictor0
filevardef(outF,"IAOVKE","float",(/"lat","lon"/))
VarAtts@long_name = "IVKE tendency due to analysis on KE (regressed)"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IAOVKE",VarAtts)
outF->IAOVKE = (/Composite/)
filevardef(outF,"IAOVKE_prob","float",(/"lat","lon"/))
VarAtts@long_name = "Probability of IAOVKE (regressed)"
VarAtts@units = ""
filevarattdef(outF, "IAOVKE_prob",VarAtts)
outF->IAOVKE_prob = (/Probability/)

