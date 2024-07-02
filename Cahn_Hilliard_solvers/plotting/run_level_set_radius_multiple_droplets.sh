indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/CPC_alpha_-0.5"
epsilon=0.14
for cohesin in 4 12 16 20 
do
   for CPC in 10 20 24 28 32
    do
        echo $cohesin
        echo $CPC
        /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "level_set_radius_multiple_droplets($CPC, $cohesin, $epsilon, '$indir', false);quit;"
    done
done


