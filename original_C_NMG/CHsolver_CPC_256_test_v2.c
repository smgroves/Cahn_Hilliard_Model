#include <stdio.h>
#include <math.h>
#include <stdlib.h>
/* #include <malloc.h> <malloc.h> is apparently in <stdlib.h> and is Linux specific */
#include <time.h>
#define gnx 256                           /* Number of grid points in x- direction defined as a global variable */
#define gny 256                           /* Number of grid points in y- direction defined as a global variable */
#define PI 4.0 * atan(1.0)                /* Defines π as a global variable */
#define iloop for (i = 1; i <= gnx; i++)  /* Increments i from 1 to gnx */
#define jloop for (j = 1; j <= gny; j++)  /* Increments j from 1 to gny */
#define ijloop iloop jloop                /* Double for loop sweeping through i and j */
#define iloopt for (i = 1; i <= nxt; i++) /* Increments i from 1 to nxt (nx temp defined locally) */
#define jloopt for (j = 1; j <= nyt; j++) /* Increments j from 1 to nyt (ny temp defined locally) */
#define ijloopt iloopt jloopt

int version = 2;
int nx, ny, n_level, c_relax;                                                                    /* Declare integers as global variables:
                                                                      max number of grid points in x- direction, max number of grid points in y- direction
                                                                      number of multigrid levels, number of SMOOTH relaxation operations */
double **ct, **sc, **smu, **sor, h, h2, dt, xleft, xright, yleft, yright, gam, Cahn, **mu, **mi; /* Declare doubles as global variables:
              Pointer-to-pointers phi at a time step, source term for phi,
              source term for µ, start of residual,
              space step size, space step size squared, ∆t,
              left and right x- coordinates, left and right y- coordinates,
              gam = ϵ, Cahn = ϵ^2 */

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

void mat_copy2(double **a, double **b, double **a2, double **b2, int xl, int xr, int yl, int yr)
{ /* Copies matrix b to matrix a
and matrix b2 to matrix a2 from xl –> xr and yl –> yr */
    int i, j;
    for (i = xl; i <= xr; i++)
        for (j = yl; j <= yr; j++)
        {
            a[i][j] = b[i][j];
            a2[i][j] = b2[i][j];
        }
}

void print_mat(FILE *fptr, double **a, int nrl, int nrh, int ncl, int nch)
{ /* Prints matrix a in a space-delimited format
from nrl –> nrh and ncl –> nch */
    int i, j;
    for (i = nrl; i <= nrh; i++)
    {
        for (j = ncl; j <= nch; j++)
            fprintf(fptr, " %16.15f", a[i][j]);
        fprintf(fptr, "\n");
    } /* Newline used to separate rows */
}

void print_data(double **phi)
{ /* Prints the conserved scalar field phi as a space-delimited file */
    FILE *fphi;
    fphi = fopen("phi_CPC_test256_tan_IC.txt", "a");
    print_mat(fphi, phi, 0, nx + 1, 0, ny + 1);
    fclose(fphi); /* Appends to an existing file */
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
    FILE *fphi_a, *fphi_f, *fphi_d2f, *fphi_cnew;
    double ht2, x_fac, y_fac, a[4], f[2], det, d2f_val, cnew_val;

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
            cnew_val = c_new[i][j];
            d2f_val = -d2f(c_new[i][j]);
            a[3] = 1.0; /* First- and second-term coefficients of Equation 23 in Mathematics 8:97 (2020) */
            f[0] = su[i][j];
            f[1] = sw[i][j] + df(c_new[i][j]) - d2f(c_new[i][j]) * c_new[i][j]; /* Third term of Equation 22 and third and fourth terms
     of Equation 23 in Mathematics 8:97 (2020) [~NOT CLEAR WHY THE SECOND TERM OF f[1] IS SO CONVOLUTED JUST TO RETURN -2*c_new[i][j]^2 */
            /* Update mu and phi according to the last terms of Equations 22 and 23 in Mathematics 8:97 (2020) [except on the very edges] */
            fphi_f = fopen("f_full_initial.csv", "a");
            fprintf(fphi_f, "%16.15f %16.15f \n", f[0], f[1]);
            fclose(fphi_f);

            fphi_cnew = fopen("su.csv", "a");
            fprintf(fphi_cnew, "%16.15f ", su[i][j]);
            if (j == nxt)
            {
                fprintf(fphi_cnew, "\n");
            }
            fclose(fphi_cnew);

            fphi_cnew = fopen("sw.csv", "a");
            fprintf(fphi_cnew, "%16.15f ", sw[i][j]); /* Outputs error */
            if (j == nxt)
            {
                fprintf(fphi_cnew, "\n");
            }
            fclose(fphi_cnew);
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
            fphi_cnew = fopen("c_new_before_update.txt", "a");
            fprintf(fphi_cnew, "%16.15f ", c_new[i][j]); /* Outputs error */
            if (j == nxt)
            {
                fprintf(fphi_cnew, "\n");
            }
            fclose(fphi_cnew);

            fphi_d2f = fopen("mu_new_before_update.txt", "a");
            fprintf(fphi_d2f, "%16.15f ", mu_new[i][j]);
            if (j == nxt)
            {
                fprintf(fphi_d2f, "\n");
            }
            fclose(fphi_d2f);
            det = a[0] * a[3] - a[1] * a[2];                 /* Calculate determinant */
            c_new[i][j] = (a[3] * f[0] - a[1] * f[1]) / det; /* Solve for the next phi */
            mu_new[i][j] = (-a[2] * f[0] + a[0] * f[1]) / det;

            -d2f(c_new[i][j]);
            fphi_f = fopen("f_full.csv", "a");
            fprintf(fphi_f, "%16.15f %16.15f \n", f[0], f[1]); /* Outputs error */
            fclose(fphi_f);

            fphi_a = fopen("a_full.csv", "a");
            fprintf(fphi_a, "%16.15f %16.15f %16.15f %16.15f  \n", a[0], a[1], a[2], a[3]); /* Outputs error */
            fclose(fphi_a);

            fphi_cnew = fopen("cnew.csv", "a");
            fprintf(fphi_cnew, "%16.15f ", cnew_val); /* Outputs error */
            if (j == nxt)
            {
                fprintf(fphi_cnew, "\n");
            }
            fclose(fphi_cnew);

            fphi_d2f = fopen("d2f.csv", "a");
            fprintf(fphi_d2f, "%16.15f ", d2f_val); /* Outputs error */
            if (j == nxt)
            {
                fprintf(fphi_d2f, "\n");
            }
            fclose(fphi_d2f);
            fphi_cnew = fopen("c_new_after_update.txt", "a");
            fprintf(fphi_cnew, "%16.15f ", c_new[i][j]); /* Outputs error */
            if (j == nxt)
            {
                fprintf(fphi_cnew, "\n");
            }
            fclose(fphi_cnew);

            fphi_d2f = fopen("mu_new_after_update.txt", "a");
            fprintf(fphi_d2f, "%16.15f ", mu_new[i][j]); /* Outputs error */
            if (j == nxt)
            {
                fprintf(fphi_d2f, "\n");
            }
            fclose(fphi_d2f);
        }
    } /* Solve for the next mu */
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

void prolong_ch(double **uc, double **uf, double **vc, double **vf, int nxc, int nyc)
{ /* Expand twofold in each direction,
passing uf and vf matrices expanded to uc and vc with original dimensions nxc and nyc; Step 7 of Page 7 of Mathematics 8:97 (2020) */
    int i, j;
    for (i = 1; i <= nxc; i++)
        for (j = 1; j <= nyc; j++)
        {
            uf[2 * i][2 * j] = uf[2 * i - 1][2 * j] = uf[2 * i][2 * j - 1] = uf[2 * i - 1][2 * j - 1] = uc[i][j];
            vf[2 * i][2 * j] = vf[2 * i - 1][2 * j] = vf[2 * i][2 * j - 1] = vf[2 * i - 1][2 * j - 1] = vc[i][j];
        }
}

void vcycle(double **uf_new, double **wf_new, double **su, double **sw, int nxf, int nyf,
            int ilevel)
{
    // FILE *wf_new_1, *wf_new_2, *wf_new_3, *wf_new_4, *wf_new_5, *wf_new_6;
    // wf_new_1 = fopen("Test_256/uf_new_1.csv", "a");
    // print_mat(wf_new_1, uf_new, 1, nxf, 1, nyf);
    // fclose(wf_new_1);
    // /* FAS multigrid cycle */
    relax(uf_new, wf_new, su, sw, ilevel, nxf, nyf); /* Relax the input data */
    // wf_new_2 = fopen("Test_256/uf_new_2.csv", "a");
    // print_mat(wf_new_2, uf_new, 1, nxf, 1, nyf);
    // fclose(wf_new_2);
    if (ilevel < n_level)
    { /* If the number of multigrid levels has not been reached */
        int nxc, nyc;
        double **duc, **dwc, **uc_new, **wc_new, **uc_def, **wc_def, **uf_def, **wf_def;
        nxc = nxf / 2;
        nyc = nyf / 2; /* Coarsen grid by twofold */
        duc = dmatrix(1, nxc, 1, nyc);
        dwc = dmatrix(1, nxc, 1, nyc); /* Initialize matrices */
        uc_new = dmatrix(1, nxc, 1, nyc);
        wc_new = dmatrix(1, nxc, 1, nyc);
        uf_def = dmatrix(1, nxf, 1, nyf);
        wf_def = dmatrix(1, nxf, 1, nyf);
        uc_def = dmatrix(1, nxc, 1, nyc);
        wc_def = dmatrix(1, nxc, 1, nyc);
        restrictCH(uf_new, uc_new, wf_new, wc_new, nxc, nyc); /* Restrict the defect upon initialization */

        // wf_new_3 = fopen("Test_256/uf_new_3.csv", "a");
        // print_mat(wf_new_3, uf_new, 1, nxf, 1, nyf);
        // fclose(wf_new_3);

        defect(duc, dwc, uf_new, wf_new, su, sw, nxf, nyf, uc_new, wc_new, nxc, nyc); /* Compute the defect */
        mat_copy2(uc_def, uc_new, wc_def, wc_new, 1, nxc, 1, nyc);                    /* Copy uc_new to uc_def and wc_new to wc_def */
        vcycle(uc_def, wc_def, duc, dwc, nxc, nyc, ilevel + 1);                       /* FAS multigrid cycle at the next coarser multigrid level */
        mat_sub2(uc_def, uc_def, uc_new, wc_def, wc_def, wc_new, 1, nxc, 1, nyc);     /* uc_def = uc_def - uc_new and wc_def = wc_def - wc_new */
        prolong_ch(uc_def, uf_def, wc_def, wf_def, nxc, nyc);                         /* Expand grid twofold; Step 7 of Page 7 of Mathematics 8:97 (2020) */
        mat_add2(uf_new, uf_new, uf_def, wf_new, wf_new, wf_def, 1, nxf, 1, nyf);     /* uf_new = uf_new + uf_def and wf_new = wf_new + wf_def; Step 8 of Page 7 of Mathematics 8:97 (2020) */

        // wf_new_4 = fopen("Test_256/uf_new_4.csv", "a");
        // print_mat(wf_new_4, uf_new, 1, nxf, 1, nyf);
        // fclose(wf_new_4);

        relax(uf_new, wf_new, su, sw, ilevel, nxf, nyf); /* Post-smoothing; Step 9 of Page 7 of Mathematics 8:97 (2020) */

        // wf_new_5 = fopen("Test_256/uf_new_5.csv", "a");
        // print_mat(wf_new_5, uf_new, 1, nxf, 1, nyf);
        // fclose(wf_new_5);

        free_dmatrix(duc, 1, nxc, 1, nyc);
        free_dmatrix(dwc, 1, nxc, 1, nyc); /* Free up pointers */
        free_dmatrix(uc_new, 1, nxc, 1, nyc);
        free_dmatrix(wc_new, 1, nxc, 1, nyc);
        free_dmatrix(uf_def, 1, nxf, 1, nyf);
        free_dmatrix(wf_def, 1, nxf, 1, nyf);
        free_dmatrix(uc_def, 1, nxc, 1, nyc);
        free_dmatrix(wc_def, 1, nxc, 1, nyc);
    }
    // wf_new_6 = fopen("Test_256/uf_new_6.csv", "a");
    // print_mat(wf_new_6, uf_new, 1, nxf, 1, nyf);
    // fclose(wf_new_6);
}

double error2(double **c_old, double **c_new, double **mu, int nxt, int nyt)
{ /* Calculates 2D residual for phi */
    int i, j;
    double **rr, res2, x = 0.0;
    FILE *fphi_rr, *fphi_rr2, *fphi_mu, *fphi_c_old;
    rr = dmatrix(1, nxt, 1, nyt);                  /* Initialize matrix */
    ijloopt { rr[i][j] = mu[i][j] - c_old[i][j]; } /* Calculate starting residual */
    // fphi_mu = fopen("Test_256/mu.csv", "a");
    // print_mat(fphi_mu, mu, 1, nxt, 1, nyt);
    // fclose(fphi_mu);

    // fphi_c_old = fopen("Test_256/c_old.csv", "a");
    // print_mat(fphi_c_old, c_old, 1, nxt, 1, nyt);
    // fclose(fphi_c_old);

    // fphi_rr = fopen("Test_256/rr_1.csv", "a");
    // print_mat(fphi_rr, rr, 1, nxt, 1, nyt);
    // fclose(fphi_rr);

    laplace(rr, sor, nx, ny);                                            /* Update with Laplace operator from rr applied to sor [~WHY IS SOR INITIALIZED OUTSIDE THIS
                                                  FUNCTION AS A GLOBAL VARIABLE INSTEAD OF LOCALLY LIKE rr? ~] */
    ijloopt { rr[i][j] = sor[i][j] - (c_new[i][j] - c_old[i][j]) / dt; } /* Update residual to reflect Laplacian step; Equation 27 of Mathematics 8:97 (2020) */
    // fphi_rr2 = fopen("Test_256/rr_2.csv", "a");
    // print_mat(fphi_rr2, rr, 1, nxt, 1, nyt);
    // fclose(fphi_rr2);
    ijloopt { x = (rr[i][j]) * (rr[i][j]) + x; } /* Calculate sum of squares */
    res2 = sqrt(x / (nx * ny));                  /* Calculate Frobenius norm */
    free_dmatrix(rr, 1, nxt, 1, nyt);            /* Free up pointers */
    return res2;
}

void initialization(double **phi)
{ /* Initialize the domain for convergence test; Equation 26 of Mathematics 8:97 (2020) */
    int i, j;
    double x, y;
    srand(time(NULL)); /* Initialize pseudorandom number generator seed */
    ijloop
    {
        x = (i - 0.5) * h;
        y = (j - 0.5) * h;
        /*phi[i][j]=cos(PI*x)*cos(PI*y);*/ /* Sinusoidal input for convergence test */
        /*if (i > round(nx/2)-5 && i < round(nx/2)+5  && j > 5 && j < ny-5) {
            phi[i][j]=-0.5*cos(0.5*(i-round(nx/2)/nx))-0.5*cos(0.5*(j-round(ny/2)/ny));}
        else {phi[i][j]=-1;}*/
        /* Setting to zero gives ringing artifacts */
        if (i > round(nx / 2) - 5 && i < round(nx / 2) + 5)
        {
            if (j > round(ny / 2) - 5 && j < round(ny / 2) + 5)
            {
                phi[i][j] = 1;
            }
            else if (i > round(nx / 2) - 1 && i < round(nx / 2) + 1)
            {
                phi[i][j] = 1;
            }
            else
            {
                phi[i][j] = -1;
            }
        }
        else
        {
            phi[i][j] = -1;
        }
    }
}

/* initialize from tanh function */
double **meshgrid(double *x, int nx, double *y, int ny)
{
    int i, j;
    double **grid = (double **)malloc(nx * sizeof(double *));
    iloop
    {
        grid[i] = (double *)malloc(ny * sizeof(double));
        jloop
        {
            grid[i][j] = x[i];
        }
    }
    jloop
    {
        iloop
        {
            grid[i][j] = y[j];
        }
    }
    return grid;
}

void free_grid(double **grid, int nx)
{
    int i;
    iloop
    {
        free(grid[i]);
    }
    free(grid);
}

double **transpose(double **xx)
{
    double **xx_t = dmatrix(1, nx, 1, ny);
    int i, j;
    ijloop
    {
        xx_t[j][i] = xx[i][j];
    }
    return xx_t;
}

double **initialize_from_function(double **phi)
{
    int i, j;
    double *x = (double *)malloc(nx * sizeof(double));
    double *y = (double *)malloc(ny * sizeof(double));
    iloop
    {
        x[i] = h * (i - 1);
        // printf("x[i]=%f\n", x[i]); /* Output model parameters */
    }
    jloop
    {
        y[j] = h * (j - 1);
        // printf("y[i]=%f\n", y[j]); /* Output model parameters */
    }
    double **xx = meshgrid(x, nx, y, ny);
    double **yy_tmp = meshgrid(x, nx, y, ny);
    double **yy = transpose(yy_tmp);
    // ijloop
    // {
    //     // printf("xx: %f, ", xx[i][j]);
    //     // printf("yy: %f (i = %d, j = %d) \n", yy[i][j], i, j);
    // }
    double R0 = 0.1;
    double eps_c = 0.01;
    ijloop
    {
        // printf("x: %f, y: %f", pow(xx[i][j] - 0.5, 2), pow(yy[i][j] - 0.5, 2));
        double R = sqrt(pow(xx[i][j] - 0.5, 2) + pow(yy[i][j] - 0.5, 2));
        double psi0 = 0.5 * (1 + tanh((R0 - R) / (2 * eps_c)));
        phi[i][j] = 2 * psi0 - 1;
        // printf("phi = %f, psi0 = %f, R=%f (i=%d, j = %d)\n", phi[i][j], psi0, R, i, j); /* Output model parameters */
    }

    free(x);
    free(y);
    free_grid(xx, nx);
    free_dmatrix(yy, 1, nx, 1, ny);
    free_grid(yy_tmp, nx);
    return phi;
}

void cahn(double **c_old, double **c_new)
{ /* Cahn-Hilliard solver */
    FILE *fphi2, *fphi_mu, *fphi_c_new;
    int i, j, max_it_CH = 10, it_mg2 = 1;
    double tol = 1.0e-10, resid2 = 1.0;
    source(c_old, sc, smu); /* Initialize source terms for phi and mu from c_old */
    while (it_mg2 <= max_it_CH && resid2 > tol)
    { /* While max iterations have not been hit and residual is greater than tolerance */
        it_mg2++;
        vcycle(c_new, mu, sc, smu, nx, ny, 1); /* Update counter and run vcycle */

        fphi_mu = fopen("mu_tan_IC.csv", "a");
        print_mat(fphi_mu, mu, 1, nx, 1, ny);
        fclose(fphi_mu);

        fphi_c_new = fopen("c_new_tan_IC.csv", "a");
        print_mat(fphi_c_new, c_new, 1, nx, 1, ny);
        fclose(fphi_c_new);

        resid2 = error2(c_old, c_new, mu, nx, ny);          /* Calculate residual error */
        printf("error2 %16.15f %d \n", resid2, it_mg2 - 1); /* Outputs error */
        fphi2 = fopen("res2_tan_IC.csv", "a");
        fprintf(fphi2, "%16.15f %d \n", resid2, it_mg2 - 1);
        fclose(fphi2);
    } /* [~NOT CLEAR WHY THIS STEP IS IN THE WHILE LOOP~] */
}

int main()
{
    int it = 1, max_it, ns, count = 1, it_mg = 1;
    double **oc, **nc, resid2 = 1.0;
    FILE *fphi, *fphi2;
    c_relax = 2;
    nx = gnx;
    ny = gny;
    n_level = (int)(log(nx) / log(2.0) + 0.1); /* Set number of relaxation cycles and number of muligrid levels */
    xleft = 0.0;
    xright = 1.0;
    yleft = 0.0;
    yright = 1.0;
    max_it = 1;
    ns = 1; /* Set x-y dimenison, max iterations, and number of steps before printing results */
    h = xright / (double)nx;
    h2 = pow(h, 2);
    gam = 4 * h / (2 * sqrt(2) * atanh(0.9));
    Cahn = pow(gam, 2);
    double c = 0.898430778152558;
    printf("%16.15f ", d2f(c));
    dt = 0.1 * h2; /* gam is gradient interfacial energy; Page 1 of Mathematics 8:97 (2020) */
    printf("nx=%d,ny=%d\n", nx, ny);
    printf("dt=%f\n", dt); /* Output model parameters */
    printf("max_it=%d\n", max_it);
    printf("ns=%d\n", ns);
    printf("n_level=%d\n\n", n_level);
    oc = dmatrix(0, nx + 1, 0, ny + 1);
    nc = dmatrix(0, nx + 1, 0, ny + 1);
    mu = dmatrix(1, nx, 1, ny); /* Initialize matrices */
    sor = dmatrix(1, nx, 1, ny);
    ct = dmatrix(1, nx, 1, ny);
    sc = dmatrix(1, nx, 1, ny);
    mi = dmatrix(1, nx, 1, ny);
    smu = dmatrix(1, nx, 1, ny);
    zero_matrix(mu, 1, nx, 1, ny);
    // initialization(oc);
    initialize_from_function(oc);
    mat_copy(nc, oc, 1, nx, 1, ny); /* Initialize oc and copy oc to nc */
    fphi = fopen("phi_CPC_test256_tan_IC.txt", "w");
    print_data(oc); /* Save initial conditions by opening a new writeable file and appending with print_data */
    fclose(fphi);

    for (it = 1; it <= max_it; it++)
    { /* Run Cahn-Hilliard solver */
        cahn(oc, nc);
        mat_copy(oc, nc, 1, nx, 1, ny); /* Run one iteration of Cahn-Hilliard and copy oc to nc */
        if (it % ns == 0)
        {
            count++;
            print_data(oc);
            printf("print out counts %d \n", count);
        } /* Every ns iterations print out counts to phi_CPC.m*/
        printf(" %d \n", it);
    }
    return 0;
}