import numpy as np
file = "/nonlinear_multigrid/julia_multigrid/manuscript_output/spinodal__same_total_time/phi_128_200_1.0e-6_dt_0.0001.txt"
def LastNlines(file_name, N):
    # opening file using with() method
    # so that file get closed
    # after completing work
    lines = np.empty((0,128))
    
    with open(file_name) as file:
        num_lines = sum(1 for _ in file)

        # loop to read iterate 
        # last n lines and print it
        for line in (file.readlines() [-N:]):
            line = np.array([float(i) for i in line[0:-2].split(" ")]).reshape(1,-1)
            print(line)
            lines = np.append(lines, line, axis = 0)
    return lines, num_lines

lines, num_lines = (LastNlines(file, 12))
print(lines.shape, num_lines)