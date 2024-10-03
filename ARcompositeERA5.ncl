; AR composite using analysis data from ARanalysisERA5.ncl
; Hing Ong, 2 Oct 2024

; Parameter setting
siglvl =   0.05
lag    =   0
latW   =  34.5
latE   =  35.5
lonS   = 154.5
lonN   = 155.5

; Add input files
fn = systemfunc("find $SCRATCH/ERA5/*/* -name 'analysis_AR_era5*' | sort -k11")
f  = addfiles(fn,"r")

; Reading predictor
IVKE       = f[:]->IVKE(:,{20:50},{130:180})
lat        = f[0]->latitude(   {20:50}          )
lon        = f[0]->longitude(           {130:180})
dsize      = dimsizes(IVKE)
tsize      = dsize(0)
latsize    = dsize(1)
lonsize    = dsize(2)

; Masking
predictor = dim_avg(dim_avg(IVKE(:,{latW:latE},{lonS:lonN})))
;predictor  = mask(predictor,AR.eq.0,False)
predictor0 = stddev(predictor)

outN = "./regression_AR.nc"
outF = addfile(outN ,"c")

dim_names1     = (/"latitude","longitude"/)
dim_sizes1     = (/latsize,lonsize/)
dim_unlimited1 = (/False,False/)
filedimdef(outF,dim_names1,dim_sizes1,dim_unlimited1)
filevardef(outF, "latitude" ,typeof(lat) ,"latitude")
outF->latitude = lat
filevardef(outF, "longitude" ,typeof(lon) ,"longitude")
outF->longitude= lon
VarAtts       = 0.0

; To do: Write a function to make the following scripts more concise.

; Reading predictand
predictand = IVKE
; Regressing
RegCoef    = regCoef_n(predictor,predictand,0,0)
TVal       = onedtond(RegCoef@tval,dimsizes(RegCoef))
DegFree    = equiv_sample_size(predictand(latitude|:,longitude|:,time|:),siglvl,0)-2
B          = TVal
B          = 0.5
Probability= betainc(DegFree/(DegFree+TVal^2),DegFree*0.5,B)
Composite  = RegCoef*predictor0
filevardef(outF,"IVKE","float",(/"latitude","longitude"/))
VarAtts@long_name = "Integrated vapor kinetic energy (regressed)"
VarAtts@units = "kg / s^2"
filevarattdef(outF, "IVKE",VarAtts)
outF->IVKE = (/Composite/)
filevardef(outF,"IVKE_prob","float",(/"latitude","longitude"/))
VarAtts@long_name = "Probability of IVKE (regressed)"
VarAtts@units = ""
filevarattdef(outF, "IVKE_prob",VarAtts)
outF->IVKE_prob = (/Probability/)

; Reading predictand
predictand = center_finite_diff_n(IVKE,3600.,False,0,0)
; Regressing
RegCoef    = regCoef_n(predictor,predictand,0,0)
TVal       = onedtond(RegCoef@tval,dimsizes(RegCoef))
DegFree    = equiv_sample_size(predictand(latitude|:,longitude|:,time|:),siglvl,0)-2
B          = TVal
B          = 0.5
Probability= betainc(DegFree/(DegFree+TVal^2),DegFree*0.5,B)
Composite  = RegCoef*predictor0
filevardef(outF,"IVKETEND","float",(/"latitude","longitude"/))
VarAtts@long_name = "IVKE tendency (regressed)"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IVKETEND",VarAtts)
outF->IVKETEND = (/Composite/)
filevardef(outF,"IVKETEND_prob","float",(/"latitude","longitude"/))
VarAtts@long_name = "Probability of IVKETEND (regressed)"
VarAtts@units = ""
filevarattdef(outF, "IVKETEND_prob",VarAtts)
outF->IVKETEND_prob = (/Probability/)

; Reading predictand
predictand = f[:]->IPEKE(:,{20:50},{130:180})
; Regressing
RegCoef    = regCoef_n(predictor,predictand,0,0)
TVal       = onedtond(RegCoef@tval,dimsizes(RegCoef))
DegFree    = equiv_sample_size(predictand(latitude|:,longitude|:,time|:),siglvl,0)-2
B          = TVal
B          = 0.5
Probability= betainc(DegFree/(DegFree+TVal^2),DegFree*0.5,B)
Composite  = RegCoef*predictor0
filevardef(outF,"IPEKE","float",(/"latitude","longitude"/))
VarAtts@long_name = "IVKE tendency due to potential energy conversion (regressed)"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IPEKE",VarAtts)
outF->IPEKE = (/Composite/)
filevardef(outF,"IPEKE_prob","float",(/"latitude","longitude"/))
VarAtts@long_name = "Probability of IPEKE (regressed)"
VarAtts@units = ""
filevarattdef(outF, "IPEKE_prob",VarAtts)
outF->IPEKE_prob = (/Probability/)

; Reading predictand
predictand = f[:]->IHAKE(:,{20:50},{130:180})
; Regressing
RegCoef    = regCoef_n(predictor,predictand,0,0)
TVal       = onedtond(RegCoef@tval,dimsizes(RegCoef))
DegFree    = equiv_sample_size(predictand(latitude|:,longitude|:,time|:),siglvl,0)-2
B          = TVal
B          = 0.5
Probability= betainc(DegFree/(DegFree+TVal^2),DegFree*0.5,B)
Composite  = RegCoef*predictor0
filevardef(outF,"IHAKE","float",(/"latitude","longitude"/))
VarAtts@long_name = "IVKE tendency due to horizontal advection of KE (regressed)"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IHAKE",VarAtts)
outF->IHAKE = (/Composite/)
filevardef(outF,"IHAKE_prob","float",(/"latitude","longitude"/))
VarAtts@long_name = "Probability of IHAKE (regressed)"
VarAtts@units = ""
filevarattdef(outF, "IHAKE_prob",VarAtts)
outF->IHAKE_prob = (/Probability/)

; Reading predictand
predictand = f[:]->IHAV(:,{20:50},{130:180})
; Regressing
RegCoef    = regCoef_n(predictor,predictand,0,0)
TVal       = onedtond(RegCoef@tval,dimsizes(RegCoef))
DegFree    = equiv_sample_size(predictand(latitude|:,longitude|:,time|:),siglvl,0)-2
B          = TVal
B          = 0.5
Probability= betainc(DegFree/(DegFree+TVal^2),DegFree*0.5,B)
Composite  = RegCoef*predictor0
filevardef(outF,"IHAV","float",(/"latitude","longitude"/))
VarAtts@long_name = "IVKE tendency due to horizontal advection of vapor (regressed)"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IHAV",VarAtts)
outF->IHAV = (/Composite/)
filevardef(outF,"IHAV_prob","float",(/"latitude","longitude"/))
VarAtts@long_name = "Probability of IHAV (regressed)"
VarAtts@units = ""
filevarattdef(outF, "IHAV_prob",VarAtts)
outF->IHAV_prob = (/Probability/)

; Reading predictand
predictand = f[:]->IVAKE(:,{20:50},{130:180})
; Regressing
RegCoef    = regCoef_n(predictor,predictand,0,0)
TVal       = onedtond(RegCoef@tval,dimsizes(RegCoef))
DegFree    = equiv_sample_size(predictand(latitude|:,longitude|:,time|:),siglvl,0)-2
B          = TVal
B          = 0.5
Probability= betainc(DegFree/(DegFree+TVal^2),DegFree*0.5,B)
Composite  = RegCoef*predictor0
filevardef(outF,"IVAKE","float",(/"latitude","longitude"/))
VarAtts@long_name = "IVKE tendency due to vertical advection of KE (regressed)"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IVAKE",VarAtts)
outF->IVAKE = (/Composite/)
filevardef(outF,"IVAKE_prob","float",(/"latitude","longitude"/))
VarAtts@long_name = "Probability of IVAKE (regressed)"
VarAtts@units = ""
filevarattdef(outF, "IVAKE_prob",VarAtts)
outF->IVAKE_prob = (/Probability/)

; Reading predictand
predictand = f[:]->IVAV(:,{20:50},{130:180})
; Regressing
RegCoef    = regCoef_n(predictor,predictand,0,0)
TVal       = onedtond(RegCoef@tval,dimsizes(RegCoef))
DegFree    = equiv_sample_size(predictand(latitude|:,longitude|:,time|:),siglvl,0)-2
B          = TVal
B          = 0.5
Probability= betainc(DegFree/(DegFree+TVal^2),DegFree*0.5,B)
Composite  = RegCoef*predictor0
filevardef(outF,"IVAV","float",(/"latitude","longitude"/))
VarAtts@long_name = "IVKE tendency due to vertical advection of vapor (regressed)"
VarAtts@units = "kg / s^3"
filevarattdef(outF, "IVAV",VarAtts)
outF->IVAV = (/Composite/)
filevardef(outF,"IVAV_prob","float",(/"latitude","longitude"/))
VarAtts@long_name = "Probability of IVAV (regressed)"
VarAtts@units = ""
filevarattdef(outF, "IVAV_prob",VarAtts)
outF->IVAV_prob = (/Probability/)

