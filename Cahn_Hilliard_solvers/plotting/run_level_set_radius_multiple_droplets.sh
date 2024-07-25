indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/CPC_alpha_-0.2"
epsilon=0.0125
for cohesin in 0.2 0.3
do
   for CPC in  0.12 0.173 0.22 0.25 0.3
    do
        echo $cohesin
        echo $CPC
        echo $epsilon
        /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "level_set_radius_multiple_droplets($CPC, $cohesin, $epsilon, '$indir', '-0.2');quit;"
    done
done


# for epsilon in [0.01, 0.02]
#     for c in [2, 6, 10] #number of grid points in 128 case
#         for CPC in [5, 10, 14, 16]

indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/CPC_alpha_0"
epsilon=0.0031
for cohesin in 0.1 0.2 0.3
do
   for CPC in 0.12 0.173 0.22 0.25 0.3
    do
        echo $cohesin
        echo $CPC
        echo $epsilon
        /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "level_set_radius_multiple_droplets($CPC, $cohesin, $epsilon, '$indir', '0');quit;"
    done
done


indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/CPC_alpha_0"
epsilon=0.0096
for cohesin in 0.1 
do
   for CPC in 0.125 
    do
        echo $cohesin
        echo $CPC
        echo $epsilon
        /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "level_set_radius_multiple_droplets($CPC, $cohesin, '$epsilon', '$indir', '0', 512, 0.000001525878906, 0.03);quit;"
    done
done


indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/CPC_alpha_0"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "level_set_radius_multiple_droplets(20, 8, '0.007504684956431058', '$indir', '0', 512, 0.00000152595, 0.03);quit;"
