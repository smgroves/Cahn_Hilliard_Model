#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include <time.h>
#define gnx 32
#define gny 32
#define PI 4.0 * atan(1.0)
#define iloop                  \
    for (i = 1; i <= gnx; i++) \
    #define jloop for (j = 1; j <= gny; j++) #define ijloop iloop jloop
#define iloopt                 \
    for (i = 1; i <= nxt; i++) \
    #define jloopt for (j = 1; j <= nyt; j++) #define ijloopt iloopt jloopt
int nx, ny, n_level, c_relax;
double **ct, **sc, **smu, **sor, h, h2, dt, xleft, xright, yleft, yright, gam, Cahn, **mu, **mi;
double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
    double **m;
    long i, nrow = nrh - nrl + 2, ncol = nch - ncl + 2;
    m = (double **)malloc((nrow) * sizeof(double *));
    m += 1;
    m -= nrl;
    m[nrl] = (double *)malloc((nrow * ncol) * sizeof(double));
    m[nrl] += 1;
    m[nrl] -= ncl;
    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;
    return m;
}
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
    free(m[nrl] + ncl - 1);
    free(m + nrl - 1);
}
void zero_matrix(double **a, int xl, int xr, int yl, int yr)
{
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
{
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
{
    int i, j;
    for (i = nrl; i <= nrh; i++)
        for (j = ncl; j <= nch; j++)
        {
            a[i][j] = b[i][j] - c[i][j];
            a2[i][j] = b2[i][j] - c2[i][j];
        }
}
void mat_copy(double **a, double **b, int xl, int xr, int yl, int yr)
{
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
{
    int i, j;
    for (i = xl; i <= xr; i++)
        for (j = yl; j <= yr; j++)
        {
            a[i][j] = b[i][j];
            a2[i][j] = b2[i][j];
        }
}
void print_mat(FILE *fptr, double **a, int nrl, int nrh, int ncl, int nch)
{
    int i, j;
    for (i = nrl; i <= nrh; i++)
    {
        for (j = ncl; j <= nch; j++)
            fprintf(fptr, " %16.15f", a[i][j]);
        fprintf(fptr, "\n");
    }
}
void print_data(double **phi)
{
    FILE *fphi;
    fphi = fopen("phi.m", "a");
    print_mat(fphi, phi, 1, nx, 1, ny);
    fclose(fphi);
}
void laplace(double **a, double **lap_a, int nxt, int nyt)
{
    int i, j;
    double ht2, dadx_L, dadx_R, dady_B, dady_T;
    ht2 = pow((xright - xleft) / (double)nxt, 2);
    ijloopt
    {
        if (i > 1)
            dadx_L = a[i][j] - a[i - 1][j];
        else
               dadx_L = 0.0;
        if (i < nxt)
            dadx_R = a[i + 1][j] - a[i][j];
        else
               dadx_R = 0.0;
        if (j > 1)
            dady_B = a[i][j] - a[i][j - 1];
        else
               dady_B = 0.0;
        if (j < nyt)
            dady_T = a[i][j + 1] - a[i][j];
        else
               dady_T = 0.0;
        lap_a[i][j] = (dadx_R - dadx_L + dady_T - dady_B) / ht2;
    }
}
void source(double **c_old, double **src_c, double **src_mu)
{
    int i, j;
    laplace(c_old, ct, nx, ny);
    ijloop
    {
        src_c[i][j] = c_old[i][j] / dt - ct[i][j];
        src_mu[i][j] = 0.0;
    }
}
double df(double c) { return pow(c, 3); }
double d2f(double c) { return 3.0 * c * c; }
void relax(double **c_new, double **mu_new, double **su, double **sw, int ilevel,
           int nxt, int nyt)
{
    int i, j, iter;
    double ht2, x_fac, y_fac, a[4], f[2], det;
    ht2 = pow((xright - xleft) / (double)nxt, 2);
    for (iter = 1; iter <= c_relax; iter++)
    {
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
            a[1] = (x_fac + y_fac) / ht2;
            a[2] = -(x_fac + y_fac) * Cahn / ht2 - d2f(c_new[i][j]);
            a[3] = 1.0;
            f[0] = su[i][j];
            f[1] = sw[i][j] + df(c_new[i][j]) - d2f(c_new[i][j]) * c_new[i][j];
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
            det = a[0] * a[3] - a[1] * a[2];
            c_new[i][j] = (a[3] * f[0] - a[1] * f[1]) / det;
            mu_new[i][j] = (-a[2] * f[0] + a[0] * f[1]) / det;
        }
    }
}
void restrictCH(double **uf, double **uc, double **vf, double **vc, int nxc, int nyc)
{
    int i, j;
    for (i = 1; i <= nxc; i++)
        for (j = 1; j <= nyc; j++)
        {
            uc[i][j] = 0.25 * (uf[2 * i][2 * j] + uf[2 * i - 1][2 * j] + uf[2 * i][2 * j - 1] + uf[2 * i - 1][2 * j - 1]);
            vc[i][j] = 0.25 * (vf[2 * i][2 * j] + vf[2 * i - 1][2 * j] + vf[2 * i][2 * j - 1] + vf[2 * i - 1][2 * j - 1]);
        }
}
void nonL(double **ru, double **rw, double **c_new, double **mu_new, int nxt, int nyt)
{
    int i, j;
    double **lap_mu, **lap_c;
    lap_mu = dmatrix(1, nxt, 1, nyt);
    lap_c = dmatrix(1, nxt, 1, nyt);
    laplace(c_new, lap_c, nxt, nyt);
    laplace(mu_new, lap_mu, nxt, nyt);
    ijloopt
    {
        ru[i][j] = c_new[i][j] / dt - lap_mu[i][j];
        rw[i][j] = mu_new[i][j] - df(c_new[i][j]) + Cahn * lap_c[i][j];
    }
    free_dmatrix(lap_mu, 1, nxt, 1, nyt);
    free_dmatrix(lap_c, 1, nxt, 1, nyt);
}
void defect(double **duc, double **dwc, double **uf_new, double **wf_new, double **suf,
            double **swf, int nxf, int nyf, double **uc_new, double **wc_new, int nxc, int nyc)
{
    double **ruf, **rwf, **rruf, **rrwf, **ruc, **rwc;
    ruc = dmatrix(1, nxc, 1, nyc);
    rwc = dmatrix(1, nxc, 1, nyc);
    ruf = dmatrix(1, nxf, 1, nyf);
    rwf = dmatrix(1, nxf, 1, nyf);
    rruf = dmatrix(1, nxc, 1, nyc);
    rrwf = dmatrix(1, nxc, 1, nyc);
    nonL(ruc, rwc, uc_new, wc_new, nxc, nyc);
    nonL(ruf, rwf, uf_new, wf_new, nxf, nyf);
    mat_sub2(ruf, suf, ruf, rwf, swf, rwf, 1, nxf, 1, nyf);
    restrictCH(ruf, rruf, rwf, rrwf, nxc, nyc);
    mat_add2(duc, ruc, rruf, dwc, rwc, rrwf, 1, nxc, 1, nyc);
    free_dmatrix(ruc, 1, nxc, 1, nyc);
    free_dmatrix(rwc, 1, nxc, 1, nyc);
    free_dmatrix(ruf, 1, nxf, 1, nyf);
    free_dmatrix(rwf, 1, nxf, 1, nyf);
    free_dmatrix(rruf, 1, nxc, 1, nyc);
    free_dmatrix(rrwf, 1, nxc, 1, nyc);
}

void prolong_ch(double **uc, double **uf, double **vc, double **vf, int nxc, int nyc)
{
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
    relax(uf_new, wf_new, su, sw, ilevel, nxf, nyf);
    if (ilevel < n_level)
    {
        int nxc, nyc;
        double **duc, **dwc, **uc_new, **wc_new, **uc_def, **wc_def, **uf_def, **wf_def;
        nxc = nxf / 2;
        nyc = nyf / 2;
        duc = dmatrix(1, nxc, 1, nyc);
        dwc = dmatrix(1, nxc, 1, nyc);
        uc_new = dmatrix(1, nxc, 1, nyc);
        wc_new = dmatrix(1, nxc, 1, nyc);
        uf_def = dmatrix(1, nxf, 1, nyf);
        wf_def = dmatrix(1, nxf, 1, nyf);
        uc_def = dmatrix(1, nxc, 1, nyc);
        wc_def = dmatrix(1, nxc, 1, nyc);
        restrictCH(uf_new, uc_new, wf_new, wc_new, nxc, nyc);
        defect(duc, dwc, uf_new, wf_new, su, sw, nxf, nyf, uc_new, wc_new, nxc, nyc);
        mat_copy2(uc_def, uc_new, wc_def, wc_new, 1, nxc, 1, nyc);
        vcycle(uc_def, wc_def, duc, dwc, nxc, nyc, ilevel + 1);
        mat_sub2(uc_def, uc_def, uc_new, wc_def, wc_def, wc_new, 1, nxc, 1, nyc);
        prolong_ch(uc_def, uf_def, wc_def, wf_def, nxc, nyc);
        mat_add2(uf_new, uf_new, uf_def, wf_new, wf_new, wf_def, 1, nxf, 1, nyf);
        relax(uf_new, wf_new, su, sw, ilevel, nxf, nyf);
        free_dmatrix(duc, 1, nxc, 1, nyc);
        free_dmatrix(dwc, 1, nxc, 1, nyc);
        free_dmatrix(uc_new, 1, nxc, 1, nyc);
        free_dmatrix(wc_new, 1, nxc, 1, nyc);
        free_dmatrix(uf_def, 1, nxf, 1, nyf);
        free_dmatrix(wf_def, 1, nxf, 1, nyf);
        free_dmatrix(uc_def, 1, nxc, 1, nyc);
        free_dmatrix(wc_def, 1, nxc, 1, nyc);
    }
}
double error2(double **c_old, double **c_new, double **mu, int nxt, int nyt)
{
    int i, j;
    double **rr, res2, x = 0.0;
    rr = dmatrix(1, nxt, 1, nyt);
    ijloopt { rr[i][j] = mu[i][j] - c_old[i][j]; }
    laplace(rr, sor, nx, ny);
    ijloopt { rr[i][j] = sor[i][j] - (c_new[i][j] - c_old[i][j]) / dt; }
    ijloopt { x = (rr[i][j]) * (rr[i][j]) + x; }
    res2 = sqrt(x / (nx * ny));
    free_dmatrix(rr, 1, nxt, 1, nyt);
    return res2;
}
void initialization(double **phi)
{
    int i, j;
    double x, y;
    ijloop
    {
        x = (i - 0.5) * h;
        y = (j - 0.5) * h;
        phi[i][j] = cos(PI * x) * cos(PI * y);
    }
}
void cahn(double **c_old, double **c_new)
{
    FILE *fphi2;
    int i, j, max_it_CH = 10000, it_mg2 = 1;
    double tol = 1.0e-10, resid2 = 1.0;
    source(c_old, sc, smu);
    while (it_mg2 <= max_it_CH && resid2 > tol)
    {
        it_mg2++;
        vcycle(c_new, mu, sc, smu, nx, ny, 1);
        resid2 = error2(c_old, c_new, mu, nx, ny);
        printf("error2 %16.15f %d \n", resid2, it_mg2 - 1);
        fphi2 = fopen("phi2.m", "a");
        fprintf(fphi2, "%16.15f %d \n", resid2, it_mg2 - 1);
        fclose(fphi2);
    }
}
int main()
{
    int it = 1, max_it, ns, count = 1, it_mg = 1;
    double **oc, **nc, resid2 = 1.0;
    FILE *fphi, *fphi2;
    c_relax = 2;
    nx = gnx;
    ny = gny;
    n_level = (int)(log(nx) / log(2.0) + 0.1);
    xleft = 0.0;
    xright = 1.0;
    yleft = 0.0;
    yright = 1.0;
    max_it = 100;
    ns = 10;
    dt = 0.01;
    h = xright / (double)nx;
    h2 = pow(h, 2);
    gam = 0.06;
    Cahn = pow(gam, 2);
    printf("nx=%d,ny=%d\n", nx, ny);
    printf("dt=%f\n", dt);
    printf("max_it=%d\n", max_it);
    printf("ns=%d\n", ns);
    printf("n_level=%d\n\n", n_level);
    oc = dmatrix(0, nx + 1, 0, ny + 1);
    nc = dmatrix(0, nx + 1, 0, ny + 1);
    mu = dmatrix(1, nx, 1, ny);
    sor = dmatrix(1, nx, 1, ny);
    ct = dmatrix(1, nx, 1, ny);
    sc = dmatrix(1, nx, 1, ny);
    mi = dmatrix(1, nx, 1, ny);
    smu = dmatrix(1, nx, 1, ny);
    zero_matrix(mu, 1, nx, 1, ny);
    initialization(oc);
    mat_copy(nc, oc, 1, nx, 1, ny);
    fphi = fopen("phi.m", "w");
    fclose(fphi);
    print_data(oc);
    for (it = 1; it <= max_it; it++)
    {
        cahn(oc, nc);
        mat_copy(oc, nc, 1, nx, 1, ny);
        if (it % ns == 0)
        {
            count++;
            print_data(oc);
            printf("print out counts %d \n", count);
        }
        printf(" %d \n", it);
    }
    return 0;
}