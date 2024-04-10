import numpy as np
import random
import os
import time
import pandas as pd
######################
# Global variables
######################
nx = 64  # max number of grid points in x-direction defined as a global variable
ny = 64  # max number of grid points in y-direction defined as a global variable

n_level = int(np.log(nx) / np.log(2.0) + 0.1)  # original c code uses natural log too
c_relax = int(2)  # number of SMOOTH relaxation operations defined as a global variable
xleft = 0.0  # left x-coordinate defined as a global variable
xright = 1.0  # right x-coordinate defined as a global variable
yleft = 0.0  # left y-coordinate defined as a global variable
yright = 1.0  # right y-coordinate defined as a global variable

# todo: check if these are needed; it appears that only ht2 (temp h^2) is used in the code
h = xright / float(nx)  # space step size defined as a global variable
h2 = h**2 #space step size squared defined as a global variable
dt = 0.1 * h2  # ∆t defined as a global variable
gam = 4 * h / (2 * np.sqrt(2) * np.arctanh(0.9))
Cahn = gam ** 2  # ϵ^2 defined as a global variable


def dmatrix(nrows, ncols):
    """
    Create an empty matrix of size nrows x ncols
    :param nrows: Number of rows
    :param ncols: Number of columns
    :return: Matrix of size nrows x ncols
    """
    return np.empty((nrows, ncols), dtype=float)


def print_data(filename, a):
    """
    Write data to file in space-delimited format
    :param filename: Name of file to write data to
    :param a: Data that is written to file
    :return: None
    """
    np.savetxt(filename, a, fmt='%16.15f', delimiter=' ')


def laplace(a, nxt, nyt, xright=xright, xleft=xleft):
    """
    Compute the discrete Laplacian of a
    :param xright:
    :param a: matrix
    :param nxt: nx temp (number of grid points in x-direction, locally defined)
    :param nyt: ny temp (number of grid points in y-direction, locally defined)
    :return: lap_a, the discrete laplacian of a
    """
    lap_a = dmatrix(nxt, nyt)
    h2 = ((xright - xleft) / nxt) ** 2
    for i in range(nxt):
        for j in range(nyt):
            if i > 0:
                dadx_L = (a[i, j] - a[i - 1, j])
            else:
                dadx_L = 0
            if i < nxt - 1:
                dadx_R = (a[i + 1, j] - a[i, j])
            else:
                dadx_R = 0
            if j > 0:
                dady_B = (a[i, j] - a[i, j - 1])
            else:
                dady_B = 0
            if j < nyt - 1:
                dady_T = (a[i, j + 1] - a[i, j])
            else:
                dady_T = 0
            lap_a[i, j] = (dadx_R - dadx_L + dady_T - dady_B) / h2
    return lap_a


# note: the C code updates ct, but it doesn't appear to be used anywhere. We'll need to check if it should be returned.
def source(c_old, nx=nx, ny=ny, dt=dt):
    """
    Compute the source term for phi and mu
    :param c_old: phi at a time step
    :return: src_c, the source term for phi, and src_mu, the source term for mu
    """

    src_mu = dmatrix(nx, ny)
    src_c = dmatrix(nx, ny)
    ct = laplace(c_old, nx, ny)
    for i in range(nx):
        for j in range(ny):
            src_c[i, j] = c_old[i, j] / dt - ct[i, j]  # update source term of phi
            src_mu[i, j] = 0  # set source term for mu to zero
    return src_c, src_mu


def df(c):
    return c ** 3


def d2f(c):
    return 3 * c ** 2


def relax(c_new, mu_new, su, sw, nxt, nyt,
          c_relax=c_relax, xright=xright, xleft=xleft, dt=dt, Cahn=Cahn):
    """
    SMOOTH Relaxation operator. This is just solving x =b*A-1 for the system of equations c_new and mu_new, where A is
    the LHS of equations 22 and 23, and b is the RHS.
    :param c_new: c to be smoothed
    :param mu_new: mu to be smoothed
    :param su: sc, locally defined
    :param sw: smu, locally defined
    :param nxt: temp (number of grid points in x-direction, locally defined)
    :param nyt: temp (number of grid points in y-direction, locally defined)
    :param c_relax: number of relaxation operations
    :param xright: right x-coordinate
    :param xleft: left x-coordinate
    :param dt: time step
    :param Cahn: ϵ^2
    :return: c_new, mu_new
    """
    ht2 = ((xright - xleft) / nxt) ** 2  # h2 temp, defined locally
    a = np.empty(4)
    f = np.empty(2)
    # print("c_new before relaxation: \n", c_new)
    # print("mu_new before relaxation: \n", mu_new)
    # print("su before relaxation: \n", su)
    # print("sw before relaxation: \n", sw)
    for iter in range(c_relax):  # c_relax is defined to be 2 in CHsolver.c
        for i in range(nxt):
            for j in range(nyt):
                if i > 0 and i < nxt - 1:
                    x_fac = 2.0
                else:
                    x_fac = 1.0
                if j > 0 and j < nyt - 1:
                    y_fac = 2.0
                else:
                    y_fac = 1.0
                a[0] = 1 / dt
                a[1] = (x_fac + y_fac) / ht2
                a[2] = -(x_fac + y_fac) * Cahn / ht2 - d2f(c_new[i][j])
                a[3] = 1.0

                f[0] = su[i][j]
                f[1] = sw[i][j] - 2 * (c_new[i][j] ** 3)  # replaced from c code with a more condensed version
                if i > 0:  # boundary cases are slightly different because i-1 doesn't exist for i = 0, for example (same for j)
                    f[0] += mu_new[i - 1][j] / ht2
                    f[1] -= Cahn * c_new[i - 1][j] / ht2
                if i < nxt - 1:
                    f[0] += mu_new[i + 1][j] / ht2
                    f[1] -= Cahn * c_new[i + 1][j] / ht2
                if j > 0:
                    f[0] += mu_new[i][j - 1] / ht2
                    f[1] -= Cahn * c_new[i][j - 1] / ht2
                if j < nyt - 1:
                    f[0] += mu_new[i][j + 1] / ht2
                    f[1] -= Cahn * c_new[i][j + 1] / ht2
                det = a[0] * a[3] - a[1] * a[2]
                c_new[i][j] = (a[3] * f[0] - a[1] * f[1]) / det
                mu_new[i][j] = (-a[2] * f[0] + a[0] * f[1]) / det
        # print("f: \n", f)
        # print("a: \n", a)
        # print("c_new: \n", c_new)
        # print("mu_new: \n", mu_new)
    return c_new, mu_new


def restrict_ch(uf, vf, nxc, nyc):
    """
    Restrict the defect twofold in each direction
    uf and vf get compressed to uc and vc with dimensions nxc and nyc
    Note that changing from C to Python requires adding 1 instead of subtracting in formulas
    :param uf: uf matrix to be restricted
    :param vf: vf matrix to be restricted
    :param nxc: number of grid points in x-direction of uc
    :param nyc: number of grid points in y-direction of vc
    :return: uc, vc
    """
    uc = dmatrix(nxc, nyc)
    vc = dmatrix(nxc, nyc)
    for i in range(nxc):
        for j in range(nyc):
            uc[i][j] = 0.25 * (
                    uf[2 * i][2 * j] + uf[2 * i + 1][2 * j] + uf[2 * i][2 * j + 1] + uf[2 * i + 1][2 * j + 1])
            vc[i][j] = 0.25 * (
                    vf[2 * i][2 * j] + vf[2 * i + 1][2 * j] + vf[2 * i][2 * j + 1] + vf[2 * i + 1][2 * j + 1])
    return uc, vc


def nonL(c_new, mu_new, nxt, nyt, dt=dt, Cahn=Cahn):
    """
    NSO operator
    :param c_new: c at a time step
    :param mu_new: mu at a time step
    :param nxt: temp (number of grid points in x-direction, locally defined)
    :param nyt: temp (number of grid points in y-direction, locally defined)
    :return: ru, rw
    """
    ru = dmatrix(nxt, nyt)
    rw = dmatrix(nxt, nyt)

    lap_c = laplace(c_new, nxt, nyt)
    lap_mu = laplace(mu_new, nxt, nyt)
    for i in range(nxt):
        for j in range(nyt):
            ru[i][j] = c_new[i][j] / dt - lap_mu[i][j]
            rw[i][j] = mu_new[i][j] - df(c_new[i][j]) + Cahn * lap_c[i][j]
    return ru, rw


def defect(uf_new, wf_new, suf, swf, nxf, nyf, uc_new, wc_new, nxc, nyc):
    ruc, rwc = nonL(uc_new, wc_new, nxc, nyc)
    ruf, rwf = nonL(uf_new, wf_new, nxf, nyf)
    ruf = suf - ruf
    rwf = swf - rwf
    rruf, rrwf = restrict_ch(ruf, rwf, nxc, nyc)
    duc = ruc + rruf
    dwc = rwc + rrwf
    return duc, dwc


def prolong_ch(uc, vc, nxc, nyc):
    uf = np.zeros((2 * nxc, 2 * nyc))
    vf = np.zeros((2 * nxc, 2 * nyc))
    for i in range(nxc):
        for j in range(nyc):
            uf[2 * i][2 * j] = uf[2 * i + 1][2 * j] = uf[2 * i][2 * j + 1] = uf[2 * i + 1][2 * j + 1] = uc[i][j]
            vf[2 * i][2 * j] = vf[2 * i + 1][2 * j] = vf[2 * i][2 * j + 1] = vf[2 * i + 1][2 * j + 1] = vc[i][j]
    return uf, vf



def vcycle(uf_new, wf_new, su, sw, nxf, nyf, ilevel):
    """
    FAS multigrid cycle
    """
    #relax the input data
    # print("before relaxing IN VCYCLE")
    # print("nxf: ", nxf) #same
    # print("uf_new: \n", uf_new) #same
    # print("wf_new: \n", wf_new)
    # print("su: \n", su)
    # print("sw: \n", sw)
    uf_new, wf_new = relax(c_new=uf_new, mu_new=wf_new, su=su, sw=sw, nxt=nxf, nyt=nyf, c_relax=c_relax, xright=xright,
                           xleft=xleft, dt=dt, Cahn=Cahn)
    # print("after relaxing IN VCYCLE")
    # print("uf_new: \n", uf_new) #same
    # print("wf_new: \n", wf_new)
    # print("su: \n", su)
    # print("sw: \n", sw)
    # If the number of multigrid levels has not been reached
    if ilevel < n_level:
        # print("ilevel: ", ilevel)
        nxc = int(nxf / 2)
        nyc = int(nyf / 2)
        uc_new, wc_new = restrict_ch(uf = uf_new, vf = wf_new, nxc = nxc, nyc = nyc)

        duc, dwc = defect(uf_new = uf_new, wf_new = wf_new, suf = su, swf = sw, nxf = nxf, nyf = nyf, uc_new=uc_new, wc_new = wc_new, nxc = nxc, nyc = nyc)

        uc_def = uc_new.copy()
        wc_def = wc_new.copy()
        # print("before recursion")
        # print("nxc: \n", nxc)
        # print("duc: \n", duc)
        # print("dwc: \n", dwc)
        # print("ilevel \n", ilevel)
        # print("uc_def: \n", uc_def)
        # print("wc_def: \n", wc_def)
        uc_def, wc_def = vcycle(uf_new = uc_def, wf_new = wc_def, su = duc, sw = dwc, nxf = nxc, nyf = nyc, ilevel = ilevel + 1)
        # print("uc_def right after vcycle")
        # print("uc_def: \n", uc_def) #different
        # print("uc_new: \n", uc_new) #same
        uc_def = uc_def - uc_new
        wc_def = wc_def - wc_new
        # print("uc_def after matrix subtraction")
        # print("uc_def: \n", uc_def) #different

        uf_def, wf_def = prolong_ch(uc = uc_def, vc = wc_def, nxc = nxc, nyc = nyc)
        # print('after prolonging \n')
        # print("nxc: \n", nxc)
        # print("uf_new: \n", uf_new)
        # print("uf_def: \n", uf_def)
        uf_new = uf_new + uf_def
        wf_new = wf_new + wf_def

        uf_new, wf_new = relax(c_new=uf_new, mu_new=wf_new, su=su, sw=sw, nxt=nxf, nyt=nyf, c_relax=c_relax, xright=xright, xleft=xleft,
              dt=dt, Cahn=Cahn)
        # print("after final relaxing IN VCYCLE")
        # print("uf_new: \n", uf_new)
        # print("wf_new: \n", wf_new)
    # print("Done with VCycle")
    return uf_new, wf_new

def error2(c_old, c_new, mu, nxt, nyt, dt = dt):
    """
    Calculate the residual for phi
    :param c_old: old phi
    :param c_new: updated phi
    :param mu: updated mu
    :param nxt: temp (number of grid points in x-direction, locally defined)
    :param nyt: temp (number of grid points in y-direction, locally defined)
    :return: res2, Frobenius norm (residual), calculated after each vcycle update to c_new and mu
    """
    rr = dmatrix(nxt, nyt)
    for i in range(nxt):
        for j in range(nyt):
            rr[i][j] = mu[i][j] - c_old[i][j]

    sor = laplace(rr,nxt,nyt)
    for i in range(nxt):
        for j in range(nyt):
            rr[i][j] = sor[i][j] - (c_new[i][j] - c_old[i][j]) / dt
    res2 = np.sqrt(np.sum(rr ** 2)/(nxt*nyt))
    return res2

def initialization(nx, ny):
    """
    Example initialization scheme for phi; this initializes phi for spinodal decomposition
    :return: phi
    """
    phi = dmatrix(nx, ny)
    CPC_width = 5
    cohesin_width = 1
    # for i in range(nx):
    #     for j in range(ny):
    #         phi[i][j] = 0.1 * (1 - 2 * random.randint(0,1))
    for i in range(nx):
        for j in range(ny):
            if i > int(nx/2 - CPC_width-1) and i < int(nx/2 + CPC_width-1):
                if j > int(ny/2 - CPC_width-1) and j < int(ny/2 + CPC_width-1):
                    phi[i][j] = 1
                elif i > int(nx/2 - cohesin_width-1) and i < int(nx/2 + cohesin_width-1):
                    phi[i][j] = 1
                else:
                    phi[i][j] = -1
            else:
                phi[i][j] = -1
    return phi

#todo: test how quickly this function runs without printing to a file
def cahn(c_old, c_new, mu, nx = nx, ny = ny, dt = dt, max_it_CH = 10000, tol = 1e-10):
    it_mg2 = 0
    resid2 = 1
    sc, smu = source(c_old, nx= nx, ny = ny, dt = dt)
    # print("sc after sourcing \n", sc)

    # with open(os.path.join(out_dir, resid_file)) as f:
    while it_mg2 < max_it_CH and resid2 > tol:
        # print("c_new before vcycle \n", c_new)
        # print("mu before vcycle \n", mu)
        # print("sc before vcycle \n", sc)
        # print("smu before vcycle \n", smu)
        c_new, mu = vcycle(uf_new = c_new, wf_new = mu, su = sc, sw = smu, nxf = nx, nyf = ny, ilevel = 1) # TODO why does this give any error when assigned to c_new, mu?
        resid2 = error2(c_old = c_old, c_new = c_new, mu = mu, nxt = nx, nyt = ny, dt = dt)
            # f.write(f"{it_mg2} error: {resid2}\n")
        it_mg2 += 1
        # print("c_new after vcycle \n", c_new)
        # print("mu after vcycle \n", mu)
        # print("sc after vcycle \n", sc)
        # print("smu after vcycle \n", smu)
    # f.close()
    return c_new

if __name__ == "__main__":
    for max_it in [1, 1000, 10000]:
        for max_it_CH in [1000,10000,100000]:
            brcd = random.randint(0,1000)
            start = time.time()
            # max_it_CH = 1
            tol = 1e-6
            # max_it = 100 # number of time steps
            ns = 10 # frequency that time step is printed to file (1 = every time step, 10 = every 10 steps, etc)
            print(f"nx = {nx}, ny = {ny}, dt = {dt}, Cahn = {Cahn}, max_it = {max_it}, ns = {ns}, n_level = {n_level}")
            mu = np.zeros((nx, ny))  # µ defined as a global variable
            oc = initialization(nx, ny)
            # print(oc)
            nc = oc.copy()

            # for it in range(max_it):
            #     nc = cahn(oc, nc, mu, max_it_CH=max_it_CH, tol = tol) # update phi
            #     oc = nc.copy() # copy updated phi to old phi

            with open(f'../outputs/runtime_tests/phi_{brcd}.txt', "w") as f:
                # write initial state to file
                for i in range(nx):
                    for j in range(ny):
                        f.write(f"{oc[i][j]} ")
                    f.write("\n")
                for it in range(max_it):
                    nc = cahn(oc, nc, mu, max_it_CH=max_it_CH, tol = tol) # update phi
                    oc = nc.copy() # copy updated phi to old phi
                    if it % ns == 0:
                        # write state to file every ns iterations
                        for i in range(nx):
                            for j in range(ny):
                                f.write(f"{oc[i][j]} ")
                            f.write("\n")
                        print(it)
            end = time.time()
            print(f"Time elapsed: {end - start} seconds")
            T = {}
            # T['brcd'] = brcd
            T['c_relax'] = c_relax
            T['Cahn'] = Cahn
            T['solver'] = "Python"
            T['dt'] = dt
            T['max_it'] = max_it
            T['max_it_CH'] = max_it_CH
            T['n_level'] = n_level
            T['ns'] = ns
            T['nx'] = nx
            T['ny'] = ny
            T['time'] = end - start
            T['tol'] = tol
            T = pd.DataFrame([T])
            file = '../Job_specs_all_py_c.csv'
            if max_it != 1:
                if not os.path.isfile(file):
                    T.to_csv(file)
                else:
                    with open(file, 'a') as f:
                        T.to_csv(f, header=False, index = False)


