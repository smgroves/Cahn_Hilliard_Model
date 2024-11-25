
# indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/droplet"
# dtout=10
# frame_rate=10
# timesteps=10000
# for dt in 6.25e-6 5.0e-5 2.5e-5 1.25e-5 0.0001
# do
#     name="phi_128_${timesteps}_1.0e-6__dt_$dt"
#     /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'contourf', 'b');quit;"
# done

# echo "Done."


# indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/CPC_geometry"
# dtout=10
# frame_rate=1
# # total_time=.25
# dt=0.000001525878906
# dtdec=$(printf "%.14f" $dt)
# # timesteps=$(echo " $total_time / $dtdec" | bc)
# timsteps=19660
# # for cohesin in 2 4
# # do
# name="phi_512_19660_1.0e-5__CPC_20_cohesin_8_eps_0.007504684956431058"
# /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', 10,'heatmap','red');quit;"

# /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, '', $frame_rate,'heatmap','red');quit;"


# # done
# echo "Done."


# indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/critical_radius_alpha_v2"
# dtout=10
# frame_rate=10
# timesteps=10000
# dt=2.5e-5
# for alpha in -0.2 -0.002 -0.4 0.0 0.2 0.002 0.4
# do
#     name="phi_128_${timesteps}_1.0e-6__alpha_${alpha}"
#     /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate,'contourf', 'blue');quit;"
# done
# echo "Done."

# indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smoothed/t_5"
# dtout=10
# frame_rate=1
# dt=6.25e-6
# name="phi_128_3200_1.0e-6_dt_6.25e-6_eps_0.015009369912862116"
# /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, '', $frame_rate,'contourf', 'blue');quit;"
# echo "Done."

# indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/CPC_geometry"
# dtout=10
# frame_rate=10
# dt=0.000001525878906
# name="incomplete_phi_1024_39320_1.0e-5__CPC_40_cohesin_16_eps_0.007504684956431058"
# /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'heatmap', 'red');quit;"
# echo "Done."

# from rivanna
# indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/CPC_geometry"
# dtout=10
# frame_rate=10
# dt=0.000001525878906
# name="phi_1024_19660_0.0001__CPC_40_cohesin_16_eps_0.007504684956431058"
# /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'heatmap', 'red');quit;"
# echo "Done."

# indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smoothed/smoothed_with_multigrid"
# dtout=10
# frame_rate=1
# dt=6.2e-6
# name="phi_128_3200_1.0e-6_dt_6.25e-6_initial_500"
# /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, '', $frame_rate, 'contourf', 'b');quit;"
# echo "Done."

# indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smoothed/t_500"
# dtout=10
# frame_rate=1
# dt=6.25e-6
# name="phi_128_3200_1.0e-6_dt_6.25e-6_periodic"
# /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, '', $frame_rate, 'contourf', 'blue');quit;"
# echo "Done."

# indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/CPC_geometry"
# dtout=10
# frame_rate=1
# dt=2.5e-5
# name="phi_256_1000_1.0e-5__CPC_10_cohesin_4_eps_0.007504684956431058"
# /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, '', $frame_rate, 'contourf', 'blue');quit;"
# echo "Done."

# indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/CPC_geometry"
# dtout=10
# frame_rate=16
# dt=0.000001525878906
# name="phi_512_19660_1.0e-5__CPC_20_cohesin_8_eps_0.007504684956431058"
# /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'contourf', 'blue');quit;"
# echo "Done."

# indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/CPC_geometry"
# dtout=10
# frame_rate=16
# dt=0.000001525878906
# name="phi_256_19660_1.0e-5__CPC_10_cohesin_4_eps_0.007504684956431058"
# /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'contourf', 'blue');quit;"
# echo "Done."

# indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/CPC_geometry"
# dtout=10
# frame_rate=16
# dt=0.000001525878906
# name="phi_128_19660_1.0e-5__CPC_5_cohesin_2_eps_0.007504684956431058"
# /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'contourf', 'blue');quit;"
# echo "Done."

indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/CPC_geometry/CPC_alpha_-0.5"
dtout=10
frame_rate=4
dt=2.5e-5
name="phi_256_2000_1.0e-5__CPC_28_cohesin_16_eps_0.14_alpha_-0.5"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'contourf', 'blue');quit;"
echo "Done."

indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/CPC_geometry/CPC_alpha_-0.5"
dtout=10
frame_rate=1
dt=2.5e-5
for c in 12 16 
do
    for CPC in 36 40 44 48 52
    do
        name="phi_256_2000_1.0e-5__CPC_${CPC}_cohesin_${c}_eps_0.08_alpha_-0.5"
        echo $name 
        /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, '', $frame_rate, 'contourf', 'blue');quit;"
    done
done
echo "Done."

#4 12 16 20
# 10 20 24 28 32 next do 4 and 32, then do numbers above for c

# 10 20 24 28 32

indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/CPC_geometry/CPC_alpha_0"
dtout=10
frame_rate=4
dt=2.5e-5
name="phi_256_2000_1.0e-5__CPC_0.173_cohesin_0.1_eps_0.0031_alpha_0"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'contourf', 'blue');quit;"
echo "Done."

indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/CPC_geometry/CPC_alpha_-0.2"
dtout=10
frame_rate=4
dt=2.5e-5
name="phi_256_2000_1.0e-5__CPC_0.173_cohesin_0.1_eps_0.0125_alpha_-0.2"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'contourf', 'blue');quit;"
echo "Done."

indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/CPC_geometry/CPC_domain_0_2"
dtout=10
frame_rate=1
dt=2.5e-5
name="phi_256_2000_1.0e-5__CPC_0.22_cohesin_0.2_eps_0.0075_alpha_0_domain_0_2"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, '', $frame_rate, 'contourf', 'blue');quit;"
echo "Done."

indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC/10_24_23_CPC_tensed_RefModel_128x64_post_transition_10_25_23_400s_post_transition_base_20Pac_0_256x256_0s_8.4max"
dtout=10
frame_rate=4
dt=2.5e-5
for eps in 0.015009
do
    for alpha in -0.5 -0.2 0.0 0.2
    do
        name="phi_256_2000_1.0e-5__eps_${eps}_alpha_${alpha}"
        echo $name 
        /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'contourf', 'blue');quit;"
    done
done
echo "Done."

indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/CPC_geometry/CPC_alpha_-0.2"
dtout=10
frame_rate=10
dt=2.5e-5
name="phi_256_2000_1.0e-5__CPC_0.173_cohesin_0.1_eps_0.0125_alpha_-0.2"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'changing_cbar', $frame_rate, 'contourf', 'blue');quit;"
echo "Done."

indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC/10_24_23_CPC_tensed_RefModel_128x64_post_transition_10_25_23_400s_post_transition_base_20Pac_0_256x256_0s_8.4max"
dtout=10
frame_rate=1
dt=2.5e-5
name="phi_256_2000_1.0e-5__eps_0.030019_alpha_-0.2"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, '', $frame_rate, 'contourf', 'blue');quit;"
echo "Done."

indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC/10_16_23_CPC_relaxed_RefModel_128x64_10_16_23_relaxed_RefModel_Mps1_phos_Plk1a_20Pac_transactiv_100_256x256_100s_8.4max"
dtout=10
frame_rate=1
dt=2.5e-5
name="phi_256_2000_1.0e-5__eps_0.030019_alpha_-0.2"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, '', $frame_rate, 'contourf', 'blue');quit;"
echo "Done."

indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/CPC_geometry/CPC_alpha_0"
dtout=10
frame_rate=1
dt=2.5e-5
name="phi_256_4000_1.0e-5__CPC_0.125_cohesin_0.1_eps_0.0096_alpha_0"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, '', $frame_rate, 'contourf', 'blue');quit;"
echo "Done."

## to do
indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC/10_16_23_CPC_relaxed_RefModel_128x64_10_16_23_relaxed_RefModel_Mps1_phos_Plk1a_20Pac_transactiv_100_256x256_100s_8.4max"
dtout=10
frame_rate=1
dt=2.5e-5
name="phi_256_2000_1.0e-5__eps_0.030019_alpha_0.2"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'shorter', $frame_rate, 'contourf', 'blue', 100);quit;"
echo "Done."

indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC/10_24_23_CPC_tensed_RefModel_128x64_post_transition_10_25_23_400s_post_transition_base_20Pac_0_256x256_0s_8.4max"
dtout=10
frame_rate=1
dt=2.5e-6
name="phi_256_4000_1.0e-5__eps_0.030019_alpha_-0.2"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, '', $frame_rate, 'contourf', 'blue', 0);quit;"
echo "Done."


indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC/10_16_23_CPC_relaxed_RefModel_128x64_10_16_23_relaxed_RefModel_Mps1_phos_Plk1a_20Pac_transactiv_100_256x256_100s_8.4max"
dtout=10
frame_rate=4
dt=2.5e-5
name="phi_256_8000_1.0e-5__eps_0.030019_alpha_-0.2"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'contourf', 'blue', 0);quit;"
echo "Done."

#todo
indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/CPC_geometry/VCell_IC/10_24_23_CPC_tensed_RefModel_128x64_post_transition_10_25_23_400s_post_transition_base_20Pac_0_256x256_0s_8.4max"
dtout=10
frame_rate=4
dt=2.5e-6
name="phi_256_8000_1.0e-5__eps_0.030019_alpha_-0.2"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'contourf', 'blue', 0);quit;"
echo "Done."

# alpha = 0 eps = 0.0096
indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/CPC_geometry/CPC_alpha_0"
dtout=10
frame_rate=10
dt=0.000001525878906
name="phi_512_19661_1.0e-5__CPC_0.125_cohesin_0.1_eps_0.0096_alpha_0"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'contourf', 'blue', 0);quit;"
echo "Done."

# alpha = 0 eps = 0.0096
indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/CPC_geometry/CPC_alpha_0"
dtout=10
frame_rate=10
dt=0.000001525878906
name="phi_512_19661_1.0e-5__CPC_0.1_cohesin_0.1_eps_0.0096_alpha_0"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'contourf', 'blue', 0);quit;"
echo "Done."

# alpha = 0 eps = 0.0096
indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/CPC_geometry/CPC_alpha_0"
dtout=10
frame_rate=10
dt=0.000001525878906
name="phi_512_19661_1.0e-5__CPC_0.12_cohesin_0.1_eps_0.0096_alpha_0"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'contourf', 'blue', 0);quit;"
echo "Done."

indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/CPC_geometry/CPC_alpha_0"
dtout=10
frame_rate=10
dt=0.000001525878906
name="phi_512_19661_1.0e-5__CPC_0.15_cohesin_0.1_eps_0.0096_alpha_0"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'contourf', 'blue', 0);quit;"
echo "Done."

# domain 0 to 2
#not done
indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/CPC_geometry/CPC_domain_0_2"
dtout=10
frame_rate=10
dt=0.000001525878906
name="phi_512_9830_1.0e-5__CPC_0.125_cohesin_0.1_eps_0.0096_alpha_0_domain_0_2"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'contourf', 'blue', 0);quit;"
echo "Done."

#done
indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/CPC_geometry/CPC_domain_0_2"
dtout=10
frame_rate=10
dt=0.000001525878906
name="phi_512_10000_1.0e-5__CPC_0.125_cohesin_0.1_eps_0.0192_alpha_0_domain_0_2"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'contourf', 'blue', 0);quit;"
echo "Done."

indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_normal_IC"
dtout=10
frame_rate=1
dt=1.25e-5
name="phi_128_2400_1.0e-6_dt_1.25e-5_mean_0_sd_0.2"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'heatmap', 'red', 0);quit;"
echo "Done."


indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_normal_IC/output"
dtout=10
frame_rate=1
dt=3e-9
name="phi_128_10000_1.0e-5_dt_3.0e-9_mean_0_sd_0.2"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'heatmap', 'red', 0);quit;"
echo "Done."


indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_normal_IC/output"
dtout=1000
frame_rate=1
dt=3e-9
name="phi_128_1000000_1.0e-5_dt_3.0e-9_mean_0_sd_0.2"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'heatmap', 'red', 0);quit;"
echo "Done."


indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_normal_IC/output"
dtout=10
frame_rate=10
dt=6.25e-6
name="phi_128_4800_1.0e-6_dt_6.25e-6_mean_0_sd_0.2"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'heatmap', 'red', 0);quit;"
echo "Done."


indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_+1_-1_IC/output"
dtout=10
frame_rate=10
dt=6.25e-6
name="phi_128_4800_1.0e-6_dt_6.25e-6"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'heatmap', 'red', 0);quit;"
echo "Done."


indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_+1_-1_IC/output"
dtout=10
frame_rate=10
dt=6.25e-6
name="phi_128_9600_1.0e-6_dt_6.25e-6"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'heatmap', 'red', 0);quit;"
echo "Done."

indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_smooth_relax_function/output"
dtout=10
frame_rate=10
dt=6.25e-6
name="phi_128_9600_1.0e-6_dt_6.25e-6"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'heatmap', 'red', 0);quit;"
echo "Done."

indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/checkerboard_IC/output"
dtout=10
frame_rate=10
dt=6.103515625e-6
name="MG_20000_dt_6.103515625e-6_Nx_128_eps_0.015009369912862116_grid_8_phi"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, '', $frame_rate, 'heatmap', 'red', 0);quit;"
echo "Done."


indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/large_and_small_droplets/output"
dtout=10
frame_rate=1
dt=6.103515625e-5
name="MG_2000_dt_6.103515625e-5_Nx_128_eps_0.015009369912862116_r1_5_r2_10_space_2_phi"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, '', $frame_rate, 'heatmap', 'red', 0);quit;"
echo "Done."

indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/large_and_small_droplets/output"
dtout=10
frame_rate=1
dt=6.103515625e-5
name="MG_2000_dt_6.103515625e-5_Nx_128_eps_0.015009369912862116_r1_20_r2_30_space_10_phi"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, '', $frame_rate, 'heatmap', 'red', 0);quit;"
echo "Done."


indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/large_and_small_droplets/output"
dtout=10
frame_rate=1
dt=0.001
name="MG_5000_dt_0.001_Nx_128_eps_0.015009369912862116_r1_20_r2_30_space_10_phi"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, '', $frame_rate, 'heatmap', 'red', 0);quit;"
echo "Done."

indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal_MG_timepoint_IC/output"
dtout=10
frame_rate=10
dt=6.25e-6
name="MG_9600_dt_6.25e-6_Nx_128_eps_0.015009369912862116_phi"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, '', $frame_rate, 'heatmap', 'red', 0);quit;"
echo "Done."

indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/matlab_multigrid/output/large_and_small_droplets"
dtout=10
frame_rate=10
dt=0.001
name="MG_MATLAB_5000_dt_1.00e-03_Nx_128_r1_20_r2_30_space_10_phi"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, '', $frame_rate, 'heatmap', 'red', 0);quit;"
echo "Done."