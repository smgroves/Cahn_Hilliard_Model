
indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/droplet"
dtout=10
frame_rate=10
timesteps=10000
for dt in 6.25e-6 5.0e-5 2.5e-5 1.25e-5 0.0001
do
    name="phi_128_${timesteps}_1.0e-6__dt_$dt"
    /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate);quit;"
done

echo "Done."


# indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/droplet_2"
# dtout=10
# frame_rate=100
# total_time=1
# for dt in 6.25e-6 5.0e-5 2.5e-5 1.25e-5 0.0001
# do
#     dtdec=$(printf "%.14f" $dt)
#     timesteps=$(echo " $total_time / $dtdec" | bc)
#     name="phi_128_${timesteps}_1.0e-6__dt_$dt"
#     /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate);quit;"
# done

# echo "Done."