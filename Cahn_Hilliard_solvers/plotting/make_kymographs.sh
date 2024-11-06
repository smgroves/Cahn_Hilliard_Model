
for time in 0 100 200 300 400
do
    sim="10_24_23_CPC_tensed_RefModel_128x64_post_transition_10_25_23_400s_post_transition_base_20Pac_${time}_256x256_${time}s_8.4max"
    indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC/$sim"
    dtout=10
    dt=2.5e-5
    for eps in 0.015009
    do
        for alpha in -0.5 -0.2 0.0 0.2
        do
            name="phi_256_2000_1.0e-5__eps_${eps}_alpha_${alpha}"
            echo $name 
            outdir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/radii_over_time_level_set_plots/$sim/$name"
            echo $outdir
            /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "kymograph_central_droplets('$indir', '$outdir','$name', $dt, $dtout, true);quit;"
        done
    done
    echo "Done."
done

for time in 100 200 300 400
do
    sim="03_25_24_CPC_relaxed_RefModel_128x64_03_25_24_relaxed_RefModel_${time}_256x256_${time}s_8.4max"
    indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC/$sim"
    dtout=10
    dt=2.5e-5
    for eps in 0.0075 0.015009
    do
        for alpha in -0.5 -0.2 0.0 0.2
        do
            name="phi_256_2000_1.0e-5__eps_${eps}_alpha_${alpha}"
            echo $name 
            outdir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/radii_over_time_level_set_plots/$sim/$name"
            echo $outdir
            /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "kymograph_central_droplets('$indir', '$outdir','$name', $dt, $dtout, true);quit;"
        done
    done
    echo "Done."
done

for time in 100 200 300 400
do
    sim="03_25_24_CPC_tensed_RefModel_128x64_04_02_24_tensed_RefModel_${time}_256x256_${time}s_8.4max"
    indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC/$sim"
    dtout=10
    dt=2.5e-5
    for eps in 0.00375 0.0075 0.015009
    do
        for alpha in -0.5 -0.2 0.0 0.2
        do
            name="phi_256_2000_1.0e-5__eps_${eps}_alpha_${alpha}"
            echo $name 
            outdir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/radii_over_time_level_set_plots/$sim/$name"
            echo $outdir
            /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "kymograph_central_droplets('$indir', '$outdir','$name', $dt, $dtout, true);quit;"
        done
    done
    echo "Done."
done

## error: '/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC/03_25_24_CPC_relaxed_RefModel_128x64_03_25_24_relaxed_RefModel_100_256x256_100s_8.4max/phi_256_2000_1.0e-5__eps_0.030019_alpha_-0.5.txt'.


for time in 100 200 300 400
do
    sim="10_16_23_CPC_relaxed_RefModel_128x64_10_16_23_relaxed_RefModel_Mps1_phos_Plk1a_20Pac_transactiv_${time}_256x256_${time}s_8.4max"
    indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC/$sim"
    dtout=10
    dt=2.5e-5
    for eps in 0.00375 0.030019 0.006 0.0075 0.015009
    do
        for alpha in -0.5 -0.2 0.0 0.2
        do
            name="phi_256_2000_1.0e-5__eps_${eps}_alpha_${alpha}"
            echo $name 
            outdir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/radii_over_time_level_set_plots/$sim/$name"
            echo $outdir
            /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "kymograph_central_droplets('$indir', '$outdir','$name', $dt, $dtout, true);quit;"
        done
    done
    echo "Done."
done


for time in 0 10 20 50

do
    sim="10_24_23_CPC_tensed_RefModel_128x64_post_transition_07_14_24_500s_post_transition_base_20Pac_${time}_256x256_${time}s_8.4max"
    indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC/$sim"
    dtout=10
    dt=2.5e-5
    for eps in 0.00375 0.030019 0.006 0.0075 0.015009 0.0125 0.04
    do
        for alpha in -0.5 -0.2 0.0 0.2
        do
            name="phi_256_2000_1.0e-5__eps_${eps}_alpha_${alpha}"
            echo $name 
            outdir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/radii_over_time_level_set_plots/$sim/$name"
            echo $outdir
            /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "kymograph_central_droplets('$indir', '$outdir','$name', $dt, $dtout, true);quit;"
        done
    done
    echo "Done."
done

sim="10_16_23_CPC_relaxed_RefModel_128x64_10_16_23_relaxed_RefModel_Mps1_phos_Plk1a_20Pac_transactiv_70_256x256_70s_8.4max"
indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC/$sim"
dtout=10
dt=2.5e-5
for eps in 0.00375 0.030019 0.006 0.0075 0.015009 0.0125 0.04
do
    for alpha in -0.5 -0.2 0.0 0.2
    do
        name="phi_256_2000_1.0e-5__eps_${eps}_alpha_${alpha}"
        echo $name 
        outdir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/radii_over_time_level_set_plots/$sim/$name"
        echo $outdir
        /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "kymograph_central_droplets('$indir', '$outdir','$name', $dt, $dtout, true);quit;"
    done
done
echo "Done."

###
sim="10_16_23_CPC_relaxed_RefModel_128x64_10_16_23_relaxed_RefModel_Mps1_phos_Plk1a_20Pac_transactiv_70_256x256_70s_8.4max"
indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC/$sim"
dtout=10
dt=2.5e-5
for eps in 0.016 0.018 0.02 0.025
do
    for alpha in -0.2
    do
        name="phi_256_2000_1.0e-5__eps_${eps}_alpha_${alpha}"
        echo $name 
        outdir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/radii_over_time_level_set_plots/$sim/$name"
        echo $outdir
        /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "kymograph_central_droplets('$indir', '$outdir','$name', $dt, $dtout, true);quit;"
    done
done
echo "Done."

sim="10_24_23_CPC_tensed_RefModel_128x64_post_transition_07_14_24_500s_post_transition_base_20Pac_0_256x256_0s_8.4max"
indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC/$sim"
dtout=10
dt=2.5e-5
for eps in 0.024 0.023
do
    for alpha in -0.2
    do
        name="phi_256_2000_1.0e-5__eps_${eps}_alpha_${alpha}"
        echo $name 
        outdir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/radii_over_time_level_set_plots/$sim/$name"
        echo $outdir
        /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "kymograph_central_droplets('$indir', '$outdir','$name', $dt, $dtout, true);quit;"
    done
done
echo "Done."




# for SASCO presentation
sim="10_16_23_CPC_relaxed_RefModel_128x64_10_16_23_relaxed_RefModel_Mps1_phos_Plk1a_20Pac_transactiv_100_256x256_100s_8.4max"
indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC/$sim"
dtout=10
dt=2.5e-6
name="phi_256_8000_1.0e-5__eps_0.030019_alpha_-0.2"
echo $name 
outdir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/radii_over_time_level_set_plots/$sim/$name"
echo $outdir
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "kymograph_central_droplets('$indir', '$outdir','$name', $dt, $dtout, true);quit;"
echo "Done."

sim="10_24_23_CPC_tensed_RefModel_128x64_post_transition_10_25_23_400s_post_transition_base_20Pac_0_256x256_0s_8.4max"
indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC/$sim"
dtout=10
dt=2.5e-6
name="phi_256_8000_1.0e-5__eps_0.030019_alpha_-0.2"
echo $name 
outdir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/radii_over_time_level_set_plots/$sim/$name"
echo $outdir
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "kymograph_central_droplets('$indir', '$outdir','$name', $dt, $dtout, true);quit;"
echo "Done."


indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/CPC_alpha_0"
dtout=10
dt=0.000001525878906
name="phi_512_19661_1.0e-5__CPC_0.125_cohesin_0.1_eps_0.0096_alpha_0"
echo $name 
outdir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/radii_over_time_level_set_plots/$name"
echo $outdir
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "kymograph_central_droplets('$indir', '$outdir','$name', $dt, $dtout, false);quit;"
echo "Done."

######################################################
# Compare domain sizes for CPC droplets
######################################################

indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/CPC_alpha_0"
dtout=10
dt=0.000001525878906
name="phi_256_19661_1.0e-5__CPC_0.125_cohesin_0.1_eps_0.0096_alpha_0"
echo $name 
outdir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/radii_over_time_level_set_plots/$name"
echo $outdir
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "kymograph_central_droplets('$indir', '$outdir','$name', $dt, $dtout, false);quit;"
echo "Done."

indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/CPC_domain_0_2"
dtout=10
dt=0.000025
name="phi_512_2000_1.0e-5__CPC_0.125_cohesin_0.1_eps_0.0192_alpha_0_domain_0_2"
echo $name 
outdir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/radii_over_time_level_set_plots/$name"
echo $outdir
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "kymograph_central_droplets('$indir', '$outdir','$name', $dt, $dtout, false);quit;"
echo "Done."


indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/CPC_domain_0_2"
dtout=10
dt=0.000025
name="phi_512_2000_1.0e-5__CPC_0.125_cohesin_0.1_eps_0.0096_alpha_0"
echo $name 
outdir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/radii_over_time_level_set_plots/$name"
echo $outdir
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "kymograph_central_droplets('$indir', '$outdir','$name', $dt, $dtout, false);quit;"
echo "Done."


indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/CPC_domain_0_2"
dtout=10
dt=0.000001525878906
name="phi_512_10000_1.0e-5__CPC_0.125_cohesin_0.1_eps_0.0192_alpha_0_domain_0_2"
echo $name 
outdir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/radii_over_time_level_set_plots/$name"
echo $outdir
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "kymograph_central_droplets('$indir', '$outdir','$name', $dt, $dtout, false);quit;"
echo "Done."


# TODO

indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/CPC_domain_0_2"
dtout=10
dt=0.000001525878906
name="phi_512_9830_1.0e-5__CPC_0.125_cohesin_0.1_eps_0.0096_alpha_0_domain_0_2"
echo $name 
outdir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/radii_over_time_level_set_plots/$name"
echo $outdir
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "kymograph_central_droplets('$indir', '$outdir','$name', $dt, $dtout, false);quit;"
echo "Done."

# for direct comparison between domain 0-1 and 0-2, need same time length
indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/CPC_alpha_0"
dtout=10
dt=0.000001525878906
name="phi_256_19661_1.0e-5__CPC_0.125_cohesin_0.1_eps_0.0096_alpha_0"
echo $name 
outdir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/radii_over_time_level_set_plots/$name"
echo $outdir
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "kymograph_central_droplets('$indir', '$outdir','$name', $dt, $dtout, false, 983);quit;"
echo "Done."

# domain change in axis labels
indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/CPC_domain_0_2"
dtout=10
dt=0.000001525878906
nx=512
name="phi_${nx}_9830_1.0e-5__CPC_0.125_cohesin_0.1_eps_0.0096_alpha_0_domain_0_2"
echo $name 
outdir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/radii_over_time_level_set_plots/$name"
echo $outdir
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "kymograph_central_droplets_domain('$indir', '$outdir','$name', $dt, $dtout, false, 0, 6.4);quit;"
echo "Done."

indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/CPC_alpha_0"
dtout=10
dt=0.000001525878906
nx=256
name="phi_${nx}_19661_1.0e-5__CPC_0.125_cohesin_0.1_eps_0.0096_alpha_0"
echo $name 
outdir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/radii_over_time_level_set_plots/$name"
echo $outdir
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "kymograph_central_droplets_domain('$indir', '$outdir','$name', $dt, $dtout, false, 983, 3.2);quit;"
echo "Done."

# job array of double domain from rivanna
indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/CPC_domain_0_2"
dtout=10
dt=0.000001525878906
nx=512
for CPC in .1 .12 0.125 0.15 0.173 0.22
do
    for cohesin in .1 0.05 0.15 0.2
    do
        name="phi_${nx}_9830_1.0e-5__CPC_${CPC}_cohesin_${cohesin}_eps_0.0096_alpha_0_domain_0_2"
        echo $name 
        outdir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/radii_over_time_level_set_plots/$name"
        echo $outdir
        /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "kymograph_central_droplets_domain('$indir', '$outdir','$name', $dt, $dtout, false, 0, 6.4);quit;"
        echo "Done."
    done
done

# job array of double domain from rivanna
indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/CPC_domain_0_2_alpha_-0.5_e_0.0075"
dtout=10
dt=0.000001525878906
nx=512
for CPC in 0.125 
do
    for cohesin in  0.05 0.12 0.13 0.14 0.15 0.2
    do
        name="phi_${nx}_19661_1.0e-5__CPC_${CPC}_cohesin_${cohesin}_eps_0.0075_alpha_-0.5_domain_0_2"
        echo $name 
        outdir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/radii_over_time_level_set_plots/$name"
        echo $outdir
        /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "kymograph_central_droplets_domain('$indir', '$outdir','$name', $dt, $dtout, false, 0, 6.4);quit;"
        echo "Done."
    done
done


