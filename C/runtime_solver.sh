# gcc -o solver_16 CHsolver_CPC_compare_to_py.c #with gnx = gny = 256; for other sizes, change in the c file and then rerun this script

#  max_it, max_it_CH
arg1=("10000")
arg2=("1000" "10000" "100000")

# Iterate over the matrices and run the C program
for a1 in "${arg1[@]}"; do
    for a2 in "${arg2[@]}"; do
        ./solver_256 "$a1" "$a2" 
    done
done