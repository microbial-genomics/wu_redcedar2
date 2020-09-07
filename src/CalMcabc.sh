#!/bin/csh

#SBATCH --job-name=MCABC_calibration
#SBATCH --partition=ord
#SBATCH --time=50:00:00
#SBACTH --ntasks=64
#SBATCH --mem=16g
#SBARCH --mail-user=wu.huiyun@epa.gov
#SBATCH --mail-type=BEGIN,END

setenv TMPDIR /work/OVERFLOW/RCR/calibration/MSU

module load intel/19.0.5
module load R/3.6.2
module load geos/3.8.0
module load gdal-2.4.3/intel-19.0
module load proj-5.2.0/intel-19.0
module load udunits-2.2.26/intel-19.0



Rscript  bac_cal.R 
