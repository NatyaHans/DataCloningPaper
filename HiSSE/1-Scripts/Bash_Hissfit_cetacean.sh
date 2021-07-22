#!/bin/sh
# here if we change sim numbers in each R script, we will get different simulation results for MLE
# So run this first and then run the second Bash file for running it all again to avoid overpopulating 
# the number of bash files. So the last run of the scripts. sim should be about 100 (if i run 100 simulations) or 100 if i run 100 simulations

for file in 2 4 8 16 32 64 128 256 512 1024 ;
   do
   echo $file;
      cp Hissefit_cetacean_forcluster.R $file".hisse_cetaceanTAGNH.R";
      sed -i -e "s/k.vec<-2/k.vec<-$file/g" $file".hisse_cetaceanTAGNH.R";
     
      # cp "Bisse_taxa"$file"_forcluster.R"  $file.$s".TAGNH.R";
      cp "sbatch_Cetaceanfit.sh"  $file".hisse_cetaceanTAGNH.sh";
      sed -i -e "s/Hissefit_cetacean_forcluster.R/$file.hisse_cetaceanTAGNH.R/g" $file".hisse_cetaceanTAGNH.sh";
      sbatch $file".hisse_cetaceanTAGNH.sh";
      done

# Dependencies: 
# Hisse_cetacean_forcluster.R
# sbatch_Cetaceanfit.sh
