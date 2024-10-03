#!/bin/bash -l

#PBS -J 20190101-20191231

### Template for running TempestExetremes in a PBS job array
./DetectBlobs --in_data "./analysis_AR_inst"$PBS_ARRAY_INDEX".nc" --out "./detect_AR_IVT"$PBS_ARRAY_INDEX".nc"  --thresholdcmd "_LAPLACIAN{9,7.0}(IVT),<,-40000.0,0.0"  --filtercmd "_LAPLACIAN{9,7.0}(IVT),<,-40000.0,50"  --minabslat "20.0"
./DetectBlobs --in_data "./analysis_AR_inst"$PBS_ARRAY_INDEX".nc" --out "./detect_AR_IVKE"$PBS_ARRAY_INDEX".nc" --thresholdcmd "_LAPLACIAN{9,7.0}(IVKE),<,-2000.0 ,0.0" --filtercmd "_LAPLACIAN{9,7.0}(IVKE),<,-2000.0 ,50" --minabslat "20.0"
