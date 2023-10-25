
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
/* #include <malloc.h> <malloc.h> is apparently in <stdlib.h> and is Linux specific */
#include <time.h>
#define gnx 3                             /* Number of grid points in x- direction defined as a global variable */
#define gny 3                             /* Number of grid points in y- direction defined as a global variable */
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
void mat_copy(double **a, double **b, int xl, int xr, int yl, int yr)
{ /* Copies matrix b to matrix a
from xl –> xr and yl –> yr */
    int i, j;
    for (i = xl; i <= xr; i++)
    {
        for (j = yl; j <= yr; j++)
        {
            a[i][j] = b[i][j];
        }
    }
}

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

void source(double **c_old, double **src_c, double **src_mu)
{ /* Calculates the spatial residuals of phi and mu as source terms to the smoothing operator
passing old phi, the spatial residual of phi, and the spatial residual of mu */
    int i, j;
    laplace(c_old, ct, nx, ny); /* Grab the Laplacian of the c_old and assign to ct given the spatial grid nx by ny */
    ijloop
    {
        src_c[i][j] = c_old[i][j] / dt - ct[i][j];
        src_mu[i][j] = 0.0;
    } /* Update spatial residual of phi and set spatial residual of mu to zero */
}

double df(double c) { return pow(c, 3); } /* Used in Eyre's unconditionally gradient stable scheme and the smoothing sweeps of phi */

double d2f(double c) { return 3.0 * c * c; } /* Used in Eyre's unconditionally gradient stable scheme and the smoothing sweeps of phi */

void relax(double **c_new, double **mu_new, double **su, double **sw, int ilevel,
           int nxt, int nyt)
{ /* Relaxation operator passing new phi, new mu, 𝝵 source term, ψ source term, [ilevel multigrid index appears to be
passed for no reason] */
    int i, j, iter;
    double ht2, x_fac, y_fac, a[4], f[2], det;
    ht2 = pow((xright - xleft) / (double)nxt, 2); /* Space step size squared [~THIS ASSUMES A SQUARE MESH GRID... NOTE THAT nyt IS NOT USED~] */
    for (iter = 1; iter <= c_relax; iter++)
    { /* For each relaxation cycle */
        ijloopt
        {
            if (i > 1 && i < nxt)
                x_fac = 2.0;
            else
                x_fac = 1.0;
            if (j > 1 && j < nyt)
                y_fac = 2.0;
            else
                y_fac = 1.0;
            a[0] = 1.0 / dt;
            a[1] = (x_fac + y_fac) / ht2; /* First- and second-term coefficients of Equation 22 in Mathematics 8:97 (2020) */
            a[2] = -(x_fac + y_fac) * Cahn / ht2 - d2f(c_new[i][j]);
            a[3] = 1.0; /* First- and second-term coefficients of Equation 23 in Mathematics 8:97 (2020) */
            f[0] = su[i][j];
            f[1] = sw[i][j] + df(c_new[i][j]) - d2f(c_new[i][j]) * c_new[i][j]; /* Third term of Equation 22 and third and fourth terms
     of Equation 23 in Mathematics 8:97 (2020) [~NOT CLEAR WHY THE SECOND TERM OF f[1] IS SO CONVOLUTED JUST TO RETURN -2*c_new[i][j]^2 */
            /* Update mu and phi according to the last terms of Equations 22 and 23 in Mathematics 8:97 (2020) [except on the very edges] */
            if (i > 1)
            {
                f[0] += mu_new[i - 1][j] / ht2;
                f[1] -= Cahn * c_new[i - 1][j] / ht2;
            }
            if (i < nxt)
            {
                f[0] += mu_new[i + 1][j] / ht2;
                f[1] -= Cahn * c_new[i + 1][j] / ht2;
            }
            if (j > 1)
            {
                f[0] += mu_new[i][j - 1] / ht2;
                f[1] -= Cahn * c_new[i][j - 1] / ht2;
            }
            if (j < nyt)
            {
                f[0] += mu_new[i][j + 1] / ht2;
                f[1] -= Cahn * c_new[i][j + 1] / ht2;
            }
            det = a[0] * a[3] - a[1] * a[2];                 /* Calculate determinant */
            c_new[i][j] = (a[3] * f[0] - a[1] * f[1]) / det; /* Solve for the next phi */
            mu_new[i][j] = (-a[2] * f[0] + a[0] * f[1]) / det;
        }
        printf("F and A matrices \n");
        printf("%f \n", f[0]);
        printf("%f \n", f[1]);
        printf("%f \n", a[0]);
        printf("%f \n", a[1]);
        printf("%f \n", a[2]);
        printf("%f \n", a[3]);
        printf("c_new and mu_new \n");
        print_mat(c_new, nxt, nyt);
        print_mat(mu_new, nxt, nyt);
    } /* Solve for the next mu */
}

int main()
{
    int it = 1, max_it, ns, count = 1, it_mg = 1;

    double **oc, **nc, resid2 = 1.0;
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
    nc = dmatrix(0, nx + 1, 0, ny + 1);
    sc = dmatrix(1, nx, 1, ny);
    smu = dmatrix(1, nx, 1, ny);
    ct = dmatrix(1, nx, 1, ny);
    mu = dmatrix(1, nx, 1, ny); /* Initialize matrices; note that oc and nc are larger than the other */

    // iterate 100 times
    // for (int c = 0; c < 100; c++)
    // {
    // printf("%d \n", c);
    make_oc_identity(oc, nx, ny);
    mat_copy(nc, oc, 1, nx, 1, ny); /* Initialize oc and copy oc to nc */

    /* print sc and smu to command line; matching print_mat by starting at 1s */
    // print_mat(oc, nx, ny);
    /* in original code, source gets called in the cahn function with c_old = oc, src_c = sc and src_mu = smu */
    source(oc, sc, smu);
    printf("sc and smu \n");
    print_mat(sc, nx, ny);
    print_mat(smu, nx, ny);
    printf("c_new, mu_new before relax \n");
    print_mat(nc, nx, ny);
    print_mat(mu, nx, ny);
    relax(nc, mu, sc, smu, 1, nx, ny);
    printf("final c_new and mu_new \n");
    print_mat(nc, nx, ny);
    print_mat(mu, nx, ny);
    // }
    // make_oc_random(oc, nx, ny);
    // print_mat(oc, nx, ny);
    // /* print sc and smu to command line; matching print_mat by starting at 1s */
    // source(oc, sc, smu);

    // print_mat(sc, nx, ny);
    // print_mat(smu, nx, ny);

    // zero_matrix(oc, 1, nx, 1, ny);
    // print_mat(oc, nx, ny);
    // source(oc, sc, smu);
    // print_mat(sc, nx, ny);
    // print_mat(smu, nx, ny);
}
