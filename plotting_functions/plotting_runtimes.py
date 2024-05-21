import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

C_in_dir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/C"
py_in_dir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/"
#
# jobs_py = pd.read_csv(f"{py_in_dir}/Job_specs.csv", index_col = 0, header = 0).reset_index()
# jobs_c = pd.read_csv(f"{C_in_dir}/Job_specs.csv", index_col = None, header = 0)
# vars_to_scan = ['max_it','max_it_CH','nx']
# vars_to_set = ['ny','ns','Cahn','dt','c_relax', 'n_level']
#
# #double check below if we want to merge on just the max_it and max_it_CH columns, or all
# compare_pd = jobs_py.merge(jobs_c, how = "outer", on = ['max_it','max_it_CH','nx','ny','c_relax'], suffixes = ["_py",'_c']).drop(['brcd','index'], axis = 1)
#
# print(compare_pd)
#
# # sns.lineplot(data = compare_pd, x = 'nx', y = 'time_py', hue = 'max_it_CH')
# # sns.lineplot(data = compare_pd, x = 'nx', y = 'time_c', hue = 'max_it_CH')
# # plt.show()
#
# # or put on the same graph colored by coding type
# #maybe need to add in other variables into id_vars to preserve them as well
# long_compare = pd.melt(compare_pd, id_vars = ['max_it','max_it_CH','nx'], value_vars = ['time_py','time_c'], var_name = 'solver', value_name = "time")
# long_compare.solver = [{'time_py':'Python', 'time_c':"C"}[i] for i in long_compare.solver]

long_compare = pd.read_csv(f"{py_in_dir}/Job_specs_all_py_c_julia.csv", index_col = None, header = 0)
long_compare['time_min'] = long_compare['time']/60

long_compare = long_compare[long_compare['tol'] == 1e-6]
long_compare = long_compare[long_compare['max_it'] == 10000]
print(long_compare)
x = 'nx'
plt.figure(dpi = 300)
long_compare_py_c = long_compare.loc[long_compare['solver'].isin(['Python','C'])]
p = sns.lineplot(data = long_compare_py_c, x = x, y = 'time_min', hue = 'solver', style = 'max_it', markers = True, dashes = True, palette = ['royalblue','green'])
# p.set(yscale='log')
plt.xscale('log', base=2)
plt.title("Runtime for Python and C CHSolvers")
plt.xlabel('Grid size (nx)')
plt.ylabel('Time (min)')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.savefig("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/plotting_functions/plots/min_py_c_10000.pdf")

if True:
    x = 'nx'
    plt.figure(dpi = 300)
    p = sns.lineplot(data = long_compare, x = x, y = 'time', hue = 'solver', style = 'max_it', markers = True, dashes = True)
    p.set(yscale='log')
    plt.xscale('log', base=2)
    plt.title("Runtime for All CHSolvers")
    plt.xlabel('Grid size (nx)')
    plt.ylabel('Time (sec)')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.savefig("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/plotting_functions/plots/log_sec_all_10000.pdf")

    x = 'nx'
    plt.figure(dpi = 300)
    long_compare_jl_c = long_compare.loc[long_compare['solver'].isin(['Julia','C'])]

    p = sns.lineplot(data = long_compare_jl_c, x = x, y = 'time', hue = 'solver', style = 'max_it', markers = True, dashes = True)
    # p.set(yscale='log')
    plt.xscale('log', base=2)
    plt.title("Runtime for Julia and C CHSolvers")
    plt.xlabel('Grid size (nx)')
    plt.ylabel('Time (sec)')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.savefig("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/plotting_functions/plots/sec_c_jl_10000.pdf")
