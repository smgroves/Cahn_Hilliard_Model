
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
/* #include <malloc.h> <malloc.h> is apparently in <stdlib.h> and is Linux specific */
#include <time.h>
#define gnx 4                             /* Number of grid points in x- direction defined as a global variable */
#define gny 4                             /* Number of grid points in y- direction defined as a global variable */
#define PI 4.0 * atan(1.0)                /* Defines π as a global variable */
#define iloop for (i = 1; i <= gnx; i++)  /* Increments i from 1 to gnx */
#define jloop for (j = 1; j <= gny; j++)  /* Increments j from 1 to gny */
#define ijloop iloop jloop                /* Double for loop sweeping through i and j */
#define iloopt for (i = 1; i <= nxt; i++) /* Increments i from 1 to nxt (nx temp defined locally) */
#define jloopt for (j = 1; j <= nyt; j++) /* Increments j from 1 to nyt (ny temp defined locally) */
#define ijloopt iloopt jloopt

int nx, ny, n_level, c_relax; /* Declare integers as global variables:
   max number of grid points in x- direction, max number of grid points in y- direction
   number of multigrid levels, number of SMOOTH relaxation operations */
double **ct, **sc, **smu, **sor, h, h2, dt, xleft, xright, yleft, yright, gam, Cahn, **mu, **mi;

double **dmatrix(long nrl, long nrh, long ncl, long nch)
{                                                       /* Discrete matrix passed:
 number of rows low (nrl), number of rows high (nrh),
 number of columns low (ncl), number of columns high (nch) */
    double **m;                                         /* Declare discrete matrix */
    long i, nrow = nrh - nrl + 2, ncol = nch - ncl + 2; /* Declare dummy variable, number of rows-columns (including flank) */
    m = (double **)malloc((nrow) * sizeof(double *));
    m += 1;
    m -= nrl; /* Allocate matrix rows starting at nrl */
    m[nrl] = (double *)malloc((nrow * ncol) * sizeof(double));
    m[nrl] += 1;
    m[nrl] -= ncl; /* Allocate matrix columns starting at ncl */
    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol; /* Add columns to each row from nrl+1 to nrh */
    return m;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
    free(m[nrl] + ncl - 1);
    free(m + nrl - 1); /* Frees up pointers created by malloc */
}

void zero_matrix(double **a, int xl, int xr, int yl, int yr)
{ /* Creates zero matrix from xl –> xr and yl –> yr */
    int i, j;
    for (i = xl; i <= xr; i++)
    {
        for (j = yl; j <= yr; j++)
        {
            a[i][j] = 0.0;
        }
    }
}

double df(double c) { return pow(c, 3); } /* Used in Eyre's unconditionally gradient stable scheme and the smoothing sweeps of phi */

double d2f(double c) { return 3.0 * c * c; } /* Used in Eyre's unconditionally gradient stable scheme and the smoothing sweeps of phi */

void laplace(double **a, double **lap_a, int nxt, int nyt)
{ /* Calculate discrete Laplace operator of matrix a,
passing nx temp and ny temp */
    int i, j;
    double ht2, dadx_L, dadx_R, dady_B, dady_T;   /* Declare local differentials */
    ht2 = pow((xright - xleft) / (double)nxt, 2); /* Space step size squared [~THIS ASSUMES A SQUARE MESH GRID... NOTE THAT nyt IS NOT USED~] */
    ijloopt
    {
        if (i > 1)
            dadx_L = a[i][j] - a[i - 1][j]; /* Calculate left differential */
        else
            dadx_L = 0.0; /* Neumann boundary conditions */
        if (i < nxt)
            dadx_R = a[i + 1][j] - a[i][j]; /* Calculate right differential */
        else
            dadx_R = 0.0; /* Neumann boundary conditions */
        if (j > 1)
            dady_B = a[i][j] - a[i][j - 1]; /* Calculate bottom differential */
        else
            dady_B = 0.0; /* Neumann boundary conditions */
        if (j < nyt)
            dady_T = a[i][j + 1] - a[i][j]; /* Calculate top differential */
        else
            dady_T = 0.0; /* Neumann boundary conditions */
        lap_a[i][j] = (dadx_R - dadx_L + dady_T - dady_B) / ht2;
    } /* Defined just below Equation 6 in Mathematics 8:97 (2020) */
}

void make_oc_random(double **oc, int nxt, int nyt)
{ /* Make initial condition for oc */
    int i, j;
    ijloopt
    {
        oc[i][j] = rand() % 10;
    }
}

void make_oc_identity(double **oc, int nxt, int nyt)
{ /* Make initial condition for oc */
    int i, j;
    ijloopt
    {
        oc[i][j] = 0.0;
    }
    ijloopt
    {
        if (i == j)
            oc[i][j] = 1.0;
    }
}

void make_oc_ones(double **oc, int nxt, int nyt)
{ /* Make initial condition for oc */
    int i, j;
    ijloopt
    {
        oc[i][j] = 1.0;
    }
    ijloopt
    {
        if (i == j)
            oc[i][j] = -1.0;
    }
}

// make oc a random matrix of 1 and -1
void make_oc_random_ones(double **oc, int nxt, int nyt)
{ /* Make initial condition for oc */
    int i, j;
    ijloopt
    {
        int random = rand() % 2;
        oc[i][j] = (random == 0) ? -1 : 1;
    }
}

void print_mat(double **a, int nxt, int nyt)
{ /* Print matrix a */
    for (int i = 1; i < nxt + 1; i++)
    {
        for (int j = 1; j < nyt + 1; j++)
        {
            printf("%f ", a[i][j]);
        }
        printf("\n");
    }
}

void nonL(double **ru, double **rw, double **c_new, double **mu_new, int nxt, int nyt)
{ /* NSO operator just above Equation 19 of Mathematics 8:97 (2020) */
    int i, j;
    double **lap_mu, **lap_c;
    lap_mu = dmatrix(1, nxt, 1, nyt);
    lap_c = dmatrix(1, nxt, 1, nyt); /* Create discrete matrices of ones */
    laplace(c_new, lap_c, nxt, nyt);
    laplace(mu_new, lap_mu, nxt, nyt); /* Initialize discrete Laplace operator */
    ijloopt
    {
        ru[i][j] = c_new[i][j] / dt - lap_mu[i][j]; /* Apply NSO to phi */
        rw[i][j] = mu_new[i][j] - df(c_new[i][j]) + Cahn * lap_c[i][j];
    } /* Apply NSO to mu */
    free_dmatrix(lap_mu, 1, nxt, 1, nyt);
    free_dmatrix(lap_c, 1, nxt, 1, nyt); /* Free up pointers */
}

void restrictCH(double **uf, double **uc, double **vf, double **vc, int nxc, int nyc)
{ /* Restrict the defect twofold in each direction,
passing uf and vf matrices compressed to uc and vc with original dimensions nxc and nyc; Step 3 of Page 7 of Mathematics 8:97 (2020) */
    int i, j;
    for (i = 1; i <= nxc; i++)
        for (j = 1; j <= nyc; j++)
        {
            uc[i][j] = 0.25 * (uf[2 * i][2 * j] + uf[2 * i - 1][2 * j] + uf[2 * i][2 * j - 1] + uf[2 * i - 1][2 * j - 1]);
            vc[i][j] = 0.25 * (vf[2 * i][2 * j] + vf[2 * i - 1][2 * j] + vf[2 * i][2 * j - 1] + vf[2 * i - 1][2 * j - 1]);
        }
}

void mat_add2(double **a, double **b, double **c, double **a2,
              double **b2, double **c2, int xl, int xr, int yl, int yr)
{ /* Adds two matrices a = b + c and a2 = b2 + c2
from xl –> xr and yl –> yr */
    int i, j;
    for (i = xl; i <= xr; i++)
        for (j = yl; j <= yr; j++)
        {
            a[i][j] = b[i][j] + c[i][j];
            a2[i][j] = b2[i][j] + c2[i][j];
        }
}

void mat_sub2(double **a, double **b, double **c, double **a2,
              double **b2, double **c2, int nrl, int nrh, int ncl, int nch)
{ /* Subtracts two matrices a = b - c and a2 = b2 - c2
from nrl –> nrh and ncl –> nch */
    int i, j;
    for (i = nrl; i <= nrh; i++)
        for (j = ncl; j <= nch; j++)
        {
            a[i][j] = b[i][j] - c[i][j];
            a2[i][j] = b2[i][j] - c2[i][j];
        }
}

void defect(double **duc, double **dwc, double **uf_new, double **wf_new, double **suf,
            double **swf, int nxf, int nyf, double **uc_new, double **wc_new, int nxc, int nyc)
{ /* Compute defect */
    double **ruf, **rwf, **rruf, **rrwf, **ruc, **rwc;
    ruc = dmatrix(1, nxc, 1, nyc);
    rwc = dmatrix(1, nxc, 1, nyc);
    ruf = dmatrix(1, nxf, 1, nyf); /* Initialized discrete matrices */
    rwf = dmatrix(1, nxf, 1, nyf);
    rruf = dmatrix(1, nxc, 1, nyc);
    rrwf = dmatrix(1, nxc, 1, nyc);
    nonL(ruc, rwc, uc_new, wc_new, nxc, nyc);
    nonL(ruf, rwf, uf_new, wf_new, nxf, nyf);
    mat_sub2(ruf, suf, ruf, rwf, swf, rwf, 1, nxf, 1, nyf);   /* ruf = suf - ruf and rwf = swf - rwf */
    restrictCH(ruf, rruf, rwf, rrwf, nxc, nyc);               /* Restrict the defect and assign to rruf and rrwf */
    mat_add2(duc, ruc, rruf, dwc, rwc, rrwf, 1, nxc, 1, nyc); /* duc = ruc + rruf and dwc = rwc + rrwf */
    free_dmatrix(ruc, 1, nxc, 1, nyc);
    free_dmatrix(rwc, 1, nxc, 1, nyc); /* Free up pointers */
    free_dmatrix(ruf, 1, nxf, 1, nyf);
    free_dmatrix(rwf, 1, nxf, 1, nyf);
    free_dmatrix(rruf, 1, nxc, 1, nyc);
    free_dmatrix(rrwf, 1, nxc, 1, nyc);
}

int main()
{
    int it = 1, max_it, ns, count = 1, it_mg = 1;
    double **ruf, **rwf, **rruf, **rrwf, **ruc, **rwc;

    double **oc, resid2 = 1.0;
    c_relax = 2;
    nx = gnx;
    ny = gny;
    xleft = 0.0;
    xright = 1.0;
    yleft = 0.0;
    yright = 1.0;
    dt = 0.01; /* Set x-y dimenison, max iterations, and number of steps before printing results */
    h = xright / (double)nx;
    h2 = pow(h, 2);
    gam = 0.06;
    Cahn = pow(gam, 2);
    oc = dmatrix(0, nx + 1, 0, ny + 1);
    mu = dmatrix(1, nx, 1, ny); /* Initialize matrices; note that oc and nc are larger than the other */
    int nxc, nyc;
    nxc = nx / 2;
    nyc = ny / 2; /* Coarsen grid by twofold */
    sc = dmatrix(1, nx, 1, ny);
    smu = dmatrix(1, nx, 1, ny);
    double **duc, **dwc, **uc_new, **wc_new, **uc_def, **wc_def, **uf_def, **wf_def;
    duc = dmatrix(1, nxc, 1, nyc);
    dwc = dmatrix(1, nxc, 1, nyc); /* Initialize matrices */
    uc_new = dmatrix(1, nxc, 1, nyc);
    wc_new = dmatrix(1, nxc, 1, nyc);
    uf_def = dmatrix(1, nx, 1, ny);
    wf_def = dmatrix(1, nx, 1, ny);
    uc_def = dmatrix(1, nxc, 1, nyc);
    wc_def = dmatrix(1, nxc, 1, nyc);

    // make_oc_identity(oc, nxc, nyc);
    // print_mat(oc, nxc, nyc);
    // print_mat(mu, nxc, nyc);

    // nonL(ruc, rwc, oc, mu, nxc, nyc);
    // print_mat(ruc, nxc, nyc);
    // print_mat(rwc, nxc, nyc);

    make_oc_random_ones(oc, nx, ny);
    make_oc_random_ones(mu, nx, ny);

    print_mat(oc, nx, ny);
    print_mat(mu, nx, ny);
    restrictCH(oc, uc_new, mu, wc_new, nxc, nyc);                        /* Restrict the defect upon initialization */
    defect(duc, dwc, oc, mu, sc, smu, nx, ny, uc_new, wc_new, nxc, nyc); /* Compute the defect */

    print_mat(duc, nxc, nyc);
    print_mat(dwc, nxc, nyc);

    // zero_matrix(oc, 1, nx, 1, ny);
    // print_mat(oc, nx, ny);

    // print_mat(sc, nx, ny);
    // print_mat(smu, nx, ny);
}
