#!/bin/csh

#SBATCH --job-name=MSU_MCABC
#SBATCH --partition=compute
#SBATCH --time=99:99:99
#SBATCH --constraint=cascadelake
#SBATCH --ntasks=48
#SBATCH --mail-user=wu.huiyun@epa.gov
#SBATCH --mail-type=BEGIN,END
#SBATCH --exclusive

setenv TMPDIR /work/OVERFLOW/RCR/sim56

module load intel/19.0.5
module load R/3.6.2
module load geos/3.8.0
module load gdal-2.4.3/intel-19.0
module load proj-5.2.0/intel-19.0
module load udunits-2.2.26/intel-19.0



Rscript  abc_loop.R
