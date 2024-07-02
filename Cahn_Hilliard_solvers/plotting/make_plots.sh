
# indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/droplet"
# dtout=10
# frame_rate=10
# timesteps=10000
# for dt in 6.25e-6 5.0e-5 2.5e-5 1.25e-5 0.0001
# do
#     name="phi_128_${timesteps}_1.0e-6__dt_$dt"
#     /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'contourf', 'b');quit;"
# done

# echo "Done."


# indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry"
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


# indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/critical_radius_alpha_v2"
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

# indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/spinodal_smoothed/t_5"
# dtout=10
# frame_rate=1
# dt=6.25e-6
# name="phi_128_3200_1.0e-6_dt_6.25e-6_eps_0.015009369912862116"
# /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, '', $frame_rate,'contourf', 'blue');quit;"
# echo "Done."

# indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry"
# dtout=10
# frame_rate=10
# dt=0.000001525878906
# name="incomplete_phi_1024_39320_1.0e-5__CPC_40_cohesin_16_eps_0.007504684956431058"
# /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'heatmap', 'red');quit;"
# echo "Done."

# from rivanna
# indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry"
# dtout=10
# frame_rate=10
# dt=0.000001525878906
# name="phi_1024_19660_0.0001__CPC_40_cohesin_16_eps_0.007504684956431058"
# /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'heatmap', 'red');quit;"
# echo "Done."

# indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/spinodal_smoothed/smoothed_with_multigrid"
# dtout=10
# frame_rate=1
# dt=6.2e-6
# name="phi_128_3200_1.0e-6_dt_6.25e-6_initial_500"
# /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, '', $frame_rate, 'contourf', 'b');quit;"
# echo "Done."

# indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/spinodal_smoothed/t_500"
# dtout=10
# frame_rate=1
# dt=6.25e-6
# name="phi_128_3200_1.0e-6_dt_6.25e-6_periodic"
# /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, '', $frame_rate, 'contourf', 'blue');quit;"
# echo "Done."

# indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry"
# dtout=10
# frame_rate=1
# dt=2.5e-5
# name="phi_256_1000_1.0e-5__CPC_10_cohesin_4_eps_0.007504684956431058"
# /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, '', $frame_rate, 'contourf', 'blue');quit;"
# echo "Done."

# indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry"
# dtout=10
# frame_rate=16
# dt=0.000001525878906
# name="phi_512_19660_1.0e-5__CPC_20_cohesin_8_eps_0.007504684956431058"
# /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'contourf', 'blue');quit;"
# echo "Done."

# indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry"
# dtout=10
# frame_rate=16
# dt=0.000001525878906
# name="phi_256_19660_1.0e-5__CPC_10_cohesin_4_eps_0.007504684956431058"
# /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'contourf', 'blue');quit;"
# echo "Done."

# indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry"
# dtout=10
# frame_rate=16
# dt=0.000001525878906
# name="phi_128_19660_1.0e-5__CPC_5_cohesin_2_eps_0.007504684956431058"
# /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'contourf', 'blue');quit;"
# echo "Done."

indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/CPC_alpha_-0.5"
dtout=10
frame_rate=4
dt=2.5e-5
name="phi_256_2000_1.0e-5__CPC_28_cohesin_16_eps_0.14_alpha_-0.5"
/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate, 'contourf', 'blue');quit;"
echo "Done."

indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/CPC_alpha_-0.5"
dtout=10
frame_rate=1
dt=2.5e-5
for c in 4 12 16 20 
do
    for CPC in 10 20 24 28 32
    do
        name="phi_256_2000_1.0e-5__CPC_${CPC}_cohesin_${c}_eps_0.14_alpha_-0.5"
        echo $name
        /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, '', $frame_rate, 'contourf', 'blue');quit;"
    done
done
echo "Done."

#4 12 16 20
# 10 20 24 28 32 next do 4 and 32, then do numbers above for c