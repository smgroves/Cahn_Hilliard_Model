
indir="/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/spinodal"
dtout=10
frame_rate=5
for dt in 6.25e-6 5.0e-5 2.5e-5 1.25e-5 0.0001
do
    name="phi_128_5000_1.0e-6_dt_$dt"
    /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "CHplotting_function('$indir', '$name', $dt, $dtout, 'fast', $frame_rate)";exit()
done

echo "Done."