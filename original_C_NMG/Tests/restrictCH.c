
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
/* #include <malloc.h> <malloc.h> is apparently in <stdlib.h> and is Linux specific */
#include <time.h>
#define gnx 8                             /* Number of grid points in x- direction defined as a global variable */
#define gny 8                             /* Number of grid points in y- direction defined as a global variable */
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

int main()
{
    int it = 1, max_it, ns, count = 1, it_mg = 1;

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
    double **uc_new, **wc_new;
    nxc = nx / 2;
    nyc = ny / 2; /* Coarsen grid by twofold */
    uc_new = dmatrix(1, nxc, 1, nyc);
    wc_new = dmatrix(1, nxc, 1, nyc);

    make_oc_identity(oc, nx, ny);
    // print_mat(oc, nx, ny);
    restrictCH(oc, uc_new, mu, wc_new, nxc, nyc);
    print_mat(uc_new, nxc, nyc);
    print_mat(wc_new, nxc, nyc);

    make_oc_random(oc, nx, ny);
    print_mat(oc, nx, ny);
    restrictCH(oc, uc_new, mu, wc_new, nxc, nyc);
    print_mat(uc_new, nxc, nyc);
    print_mat(wc_new, nxc, nyc);

    // zero_matrix(oc, 1, nx, 1, ny);
    // print_mat(oc, nx, ny);

    // print_mat(sc, nx, ny);
    // print_mat(smu, nx, ny);
}
