/*
 * ============================================================================
 *  Méthode 3 (répertoire 04_VF_TPFA) : Volumes Finis classiques à deux points
 *  (VF TPFA sur cellules de Voronoï) pour -div(K grad u) = f
 *  sur le domaine NON CONVEXE en L = [0,1]x[0,2] U [1,2]x[0,1]
 *
 *  Objectif scientifique :
 *    - Vérifier que TPFA converge en O(h^2) sur tenseur ISOTROPE (Ex. 0 : K = I)
 *    - Mettre en évidence l'INCONSISTANCE INTRINSÈQUE de TPFA dès que K est un
 *      tenseur plein (Ex. 1, Ex. 2, Ex. 3) : la matrice reste symétrique et
 *      conservative, mais l'erreur stagne (ne converge pas), même en adaptatif !
 *
 *  Options CLI :
 *      ./vf_tpfa_lshape                  (temps réel Gnuplot, Ex. 0..3)
 *      ./vf_tpfa_lshape ex=3 -o ex3.gif  (enregistrement animation GIF)
 *      ./vf_tpfa_lshape gui=0            (mode console seul)
 *
 *  Compilation :
 *      gcc -O3 -std=c99 -Wall -Wextra -Wpedantic -o vf_tpfa_lshape vf_tpfa_lshape.c -lm
 * ============================================================================
 */

#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifdef _WIN32
  #include <windows.h>
  #define SLEEP_MS(ms) Sleep(ms)
  #define POPEN _popen
  #define PCLOSE _pclose
#else
  #include <unistd.h>
  #define SLEEP_MS(ms) do { struct timespec _ts = {(ms)/1000, ((ms)%1000)*1000000L}; nanosleep(&_ts, NULL); } while (0)
  #define POPEN popen
  #define PCLOSE pclose
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * 1. Utilitaires mémoire, temps, géométrie L-shape
 * ============================================================================ */

static void *xmalloc(size_t n) {
    void *p = malloc(n ? n : 1);
    if (!p) { fprintf(stderr, "Erreur: memoire insuffisante\n"); exit(1); }
    return p;
}

static void *xcalloc(size_t n, size_t s) {
    void *p = calloc(n ? n : 1, s);
    if (!p) { fprintf(stderr, "Erreur: memoire insuffisante\n"); exit(1); }
    return p;
}

static void *xrealloc(void *p, size_t n) {
    void *q = realloc(p, n ? n : 1);
    if (!q) { fprintf(stderr, "Erreur: memoire insuffisante\n"); exit(1); }
    return q;
}

static double now_sec(void) { return (double)clock() / (double)CLOCKS_PER_SEC; }

static double orient2d(double ax, double ay, double bx, double by, double cx, double cy) {
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
}

static int in_lshape(double x, double y, double tol) {
    if (x < -tol || x > 2.0 + tol || y < -tol || y > 2.0 + tol) return 0;
    if (x > 1.0 + tol && y > 1.0 + tol) return 0;
    return 1;
}

static int on_lshape_bdry(double x, double y) {
    const double e = 1e-10;
    if (fabs(x) < e || fabs(y) < e) return 1;
    if (fabs(x - 2.0) < e && y <= 1.0 + e) return 1;
    if (fabs(y - 2.0) < e && x <= 1.0 + e) return 1;
    if (fabs(x - 1.0) < e && y >= 1.0 - e) return 1;
    if (fabs(y - 1.0) < e && x >= 1.0 - e) return 1;
    return 0;
}

/* ============================================================================
 * 2. Cas tests sur domaine non convexe en L (Ex. 0 isotrope + Ex. 1..3 plein)
 * ============================================================================ */

static int EX_ID = 1;

static void eval_K(double x, double K[2][2]) {
    if (EX_ID == 0) {
        /* Témoin isotrope : K = I (TPFA est consistant et converge en O(h^2)) */
        K[0][0] = 1.0; K[0][1] = 0.0;
        K[1][0] = 0.0; K[1][1] = 1.0;
    } else if (EX_ID == 1) {
        K[0][0] =  2.0; K[0][1] = -1.0;
        K[1][0] = -1.0; K[1][1] =  2.0;
    } else if (EX_ID == 2) {
        if (x < 1.0) { K[0][0] = 1.0; K[0][1] = 0.0; K[1][0] = 0.0; K[1][1] = 1.0; }
        else         { K[0][0] = 8.0; K[0][1] = -7.0; K[1][0] = -7.0; K[1][1] = 20.0; }
    } else {
        double th = 30.0 * (M_PI / 180.0);
        double c = cos(th), s = sin(th), l1 = 50.0, l2 = 1.0;
        K[0][0] = l1 * c * c + l2 * s * s;
        K[0][1] = (l1 - l2) * c * s;
        K[1][0] = K[0][1];
        K[1][1] = l1 * s * s + l2 * c * c;
    }
}

static void uex_d(double x, double y, double d[6]) {
    if (EX_ID == 0 || EX_ID == 1) {
        double k = 0.5 * M_PI;
        double sx = sin(k * x), cx = cos(k * x), sy = sin(k * y), cy = cos(k * y);
        d[0] = sx * sy + 0.25 * x * y;
        d[1] = k * cx * sy + 0.25 * y;
        d[2] = k * sx * cy + 0.25 * x;
        d[3] = -k * k * sx * sy;
        d[4] =  k * k * cx * cy + 0.25;
        d[5] = -k * k * sx * sy;
    } else if (EX_ID == 2) {
        double a = 0.5, e = exp(a * (x + y));
        d[0] = e; d[1] = a * e; d[2] = a * e;
        d[3] = a * a * e; d[4] = a * a * e; d[5] = a * a * e;
    } else {
        double dx = x - 0.92, dy = y - 0.92, a = 35.0;
        double e = exp(-a * (dx * dx + dy * dy));
        d[0] = e + 0.2 * (x + y);
        d[1] = -2.0 * a * dx * e + 0.2;
        d[2] = -2.0 * a * dy * e + 0.2;
        d[3] = (4.0 * a * a * dx * dx - 2.0 * a) * e;
        d[4] =  4.0 * a * a * dx * dy * e;
        d[5] = (4.0 * a * a * dy * dy - 2.0 * a) * e;
    }
}

static double u_ex(double x, double y) { double d[6]; uex_d(x, y, d); return d[0]; }

static double f_src(double x, double y) {
    double K[2][2], d[6];
    eval_K(x, K);
    uex_d(x, y, d);
    return -(K[0][0] * d[3] + 2.0 * K[0][1] * d[4] + K[1][1] * d[5]);
}

/* ============================================================================
 * 3. Maillage Delaunay (Bowyer-Watson) filtré sur le domaine en L + Voronoï
 * ============================================================================ */

typedef struct { int n, cap; double *x, *y; } Pts;

static void pts_init(Pts *P, int cap) {
    P->n = 0; P->cap = cap > 32 ? cap : 32;
    P->x = (double *)xmalloc((size_t)P->cap * sizeof(double));
    P->y = (double *)xmalloc((size_t)P->cap * sizeof(double));
}

static void pts_free(Pts *P) { free(P->x); free(P->y); P->x = P->y = NULL; P->n = P->cap = 0; }

static void pts_add(Pts *P, double x, double y) {
    if (P->n >= P->cap) {
        P->cap *= 2;
        P->x = (double *)xrealloc(P->x, (size_t)P->cap * sizeof(double));
        P->y = (double *)xrealloc(P->y, (size_t)P->cap * sizeof(double));
    }
    P->x[P->n] = x; P->y[P->n] = y; P->n++;
}

static void pts_lshape_grid(Pts *P, int n1) {
    pts_init(P, 3 * (n1 + 1) * (n1 + 1));
    double h = 1.0 / (double)n1;
    for (int i = 0; i <= 2 * n1; i++) {
        for (int j = 0; j <= 2 * n1; j++) {
            double x = i * h, y = j * h;
            if (x > 1.0 + 1e-11 && y > 1.0 + 1e-11) continue;
            pts_add(P, x, y);
        }
    }
}

typedef struct { int v[3]; double cx, cy, r2; int alive; } BWTri;

static void bw_circum(const double *x, const double *y, BWTri *t) {
    double ax = x[t->v[0]], ay = y[t->v[0]];
    double bx = x[t->v[1]] - ax, by = y[t->v[1]] - ay;
    double cx = x[t->v[2]] - ax, cy = y[t->v[2]] - ay;
    double d = 2.0 * (bx * cy - by * cx);
    if (fabs(d) < 1e-24) { t->cx = ax; t->cy = ay; t->r2 = 1e30; return; }
    double b2 = bx * bx + by * by, c2 = cx * cx + cy * cy;
    double ux = (cy * b2 - by * c2) / d, uy = (bx * c2 - cx * b2) / d;
    t->cx = ax + ux; t->cy = ay + uy; t->r2 = ux * ux + uy * uy;
}

static void build_delaunay_lshape(int np, const double *px, const double *py,
                                  int **tri_out, int *nt_out) {
    double *x = (double *)xmalloc((size_t)(np + 3) * sizeof(double));
    double *y = (double *)xmalloc((size_t)(np + 3) * sizeof(double));
    memcpy(x, px, (size_t)np * sizeof(double));
    memcpy(y, py, (size_t)np * sizeof(double));

    double M = 40.0;
    x[np] = -M; y[np] = -M;
    x[np + 1] = 2 * M; y[np + 1] = -M;
    x[np + 2] = 1.0; y[np + 2] = 2 * M;

    int cap = 4 * np + 32, nt = 1;
    BWTri *T = (BWTri *)xmalloc((size_t)cap * sizeof(BWTri));
    T[0].v[0] = np; T[0].v[1] = np + 1; T[0].v[2] = np + 2; T[0].alive = 1;
    bw_circum(x, y, &T[0]);

    int *bad = (int *)xmalloc((size_t)cap * sizeof(int));
    int (*poly)[2] = (int (*)[2])xmalloc((size_t)(3 * cap) * sizeof(int[2]));

    for (int p = 0; p < np; p++) {
        int nb = 0;
        for (int t = 0; t < nt; t++) {
            if (!T[t].alive) continue;
            double dx = x[p] - T[t].cx, dy = y[p] - T[t].cy;
            if (dx * dx + dy * dy <= T[t].r2 * (1.0 + 1e-11)) bad[nb++] = t;
        }
        int ne = 0;
        for (int i = 0; i < nb; i++) {
            int t = bad[i];
            for (int e = 0; e < 3; e++) {
                int a = T[t].v[e], b = T[t].v[(e + 1) % 3], shared = 0;
                for (int j = 0; j < nb; j++) {
                    if (j == i) continue;
                    int u = bad[j];
                    for (int f = 0; f < 3; f++) {
                        if (T[u].v[f] == b && T[u].v[(f + 1) % 3] == a) { shared = 1; break; }
                    }
                    if (shared) break;
                }
                if (!shared) { poly[ne][0] = a; poly[ne][1] = b; ne++; }
            }
        }
        for (int i = 0; i < nb; i++) T[bad[i]].alive = 0;
        for (int e = 0; e < ne; e++) {
            if (nt >= cap) {
                cap *= 2;
                T = (BWTri *)xrealloc(T, (size_t)cap * sizeof(BWTri));
                bad = (int *)xrealloc(bad, (size_t)cap * sizeof(int));
                poly = (int (*)[2])xrealloc(poly, (size_t)(3 * cap) * sizeof(int[2]));
            }
            int a = poly[e][0], b = poly[e][1];
            if (orient2d(x[a], y[a], x[b], y[b], x[p], y[p]) < 0.0) { int tmp = a; a = b; b = tmp; }
            T[nt].v[0] = a; T[nt].v[1] = b; T[nt].v[2] = p; T[nt].alive = 1;
            bw_circum(x, y, &T[nt]);
            nt++;
        }
    }

    int out_n = 0;
    int *tri = (int *)xmalloc((size_t)(3 * nt) * sizeof(int));
    for (int t = 0; t < nt; t++) {
        if (!T[t].alive) continue;
        int a = T[t].v[0], b = T[t].v[1], c = T[t].v[2];
        if (a >= np || b >= np || c >= np) continue;
        if (fabs(orient2d(x[a], y[a], x[b], y[b], x[c], y[c])) < 1e-14) continue;
        /* Filtrage strict : rejeter tout triangle dans le coin exclu (1,2)x(1,2) */
        double mx = (x[a] + x[b] + x[c]) / 3.0;
        double my = (y[a] + y[b] + y[c]) / 3.0;
        if (!in_lshape(mx, my, -1e-9)) continue;
        tri[3 * out_n] = a; tri[3 * out_n + 1] = b; tri[3 * out_n + 2] = c;
        out_n++;
    }
    free(x); free(y); free(T); free(bad); free(poly);
    *tri_out = tri; *nt_out = out_n;
}

typedef struct {
    int     nn, nt;
    double *x, *y;
    int    *bdry;
    int    *tri;
    int    *et;       /* adjacence triangle-triangle */
    double *area_v;   /* aire |V_i| du volume de Voronoï dual */
} Mesh;


/* ---- qualite du maillage : angle minimal (deg) et nombre de triangles obtus (> 90 deg) ---- */
static int PRINT_PTS = 0;      /* pts=1 : affiche le tableau u_h - u_exact en chaque noeud */
static void mesh_quality(const Mesh *M, double *amin_deg, int *n_obtus) {
    double amin = 180.0; int nob = 0;
    for (int t = 0; t < M->nt; t++) {
        int v[3] = { M->tri[3 * t], M->tri[3 * t + 1], M->tri[3 * t + 2] }, obt = 0;
        for (int k = 0; k < 3; k++) {
            int a = v[k], b = v[(k + 1) % 3], c = v[(k + 2) % 3];
            double ux = M->x[b] - M->x[a], uy = M->y[b] - M->y[a], wx = M->x[c] - M->x[a], wy = M->y[c] - M->y[a];
            double ang = acos(fmax(-1.0, fmin(1.0, (ux * wx + uy * wy) / (hypot(ux, uy) * hypot(wx, wy))))) * 180.0 / M_PI;
            if (ang < amin) amin = ang;
            if (ang > 90.0 + 1e-9) obt = 1;
        }
        nob += obt;
    }
    *amin_deg = amin; *n_obtus = nob;
}
/* erreur ponctuelle u_h - u_exact en chaque noeud (fichier + affichage optionnel) */
static double point_errors(const Mesh *M, const double *u, int ex_id, int final) {
    double emax = 0.0;
    FILE *f = NULL;
    if (final) { char fn[64]; snprintf(fn, sizeof(fn), "ex%d_erreur_points.dat", ex_id); f = fopen(fn, "w"); }
    if (f) fprintf(f, "# i x y u_h u_exact u_h-u_exact\n");
    if (final && PRINT_PTS) printf("   %6s %10s %10s %14s %14s %12s\n", "i", "x", "y", "u_h", "u_exact", "u_h-u_ex");
    for (int i = 0; i < M->nn; i++) {
        double ue = u_ex(M->x[i], M->y[i]), d = u[i] - ue;
        if (fabs(d) > emax) emax = fabs(d);
        if (f) fprintf(f, "%d %.8f %.8f %.10e %.10e %.4e\n", i, M->x[i], M->y[i], u[i], ue, d);
        if (final && PRINT_PTS) printf("   %6d %10.5f %10.5f %14.8f %14.8f %12.3e\n", i, M->x[i], M->y[i], u[i], ue, d);
    }
    if (f) fclose(f);
    return emax;
}
typedef struct { int a, b, t, e; } EdgeKey;
static int cmp_ek(const void *p, const void *q) {
    const EdgeKey *u = (const EdgeKey *)p, *v = (const EdgeKey *)q;
    return (u->a != v->a) ? (u->a - v->a) : (u->b - v->b);
}

static void mesh_build(Mesh *M, const Pts *P) {
    int nn = P->n;
    M->nn = nn;
    M->x = (double *)xmalloc((size_t)nn * sizeof(double));
    M->y = (double *)xmalloc((size_t)nn * sizeof(double));
    M->bdry = (int *)xcalloc((size_t)nn, sizeof(int));
    memcpy(M->x, P->x, (size_t)nn * sizeof(double));
    memcpy(M->y, P->y, (size_t)nn * sizeof(double));
    for (int i = 0; i < nn; i++) M->bdry[i] = on_lshape_bdry(M->x[i], M->y[i]);

    build_delaunay_lshape(nn, M->x, M->y, &M->tri, &M->nt);
    int nt = M->nt;

    EdgeKey *ek = (EdgeKey *)xmalloc((size_t)(3 * nt) * sizeof(EdgeKey));
    for (int t = 0; t < nt; t++) {
        for (int e = 0; e < 3; e++) {
            int a = M->tri[3 * t + e], b = M->tri[3 * t + (e + 1) % 3];
            ek[3 * t + e].a = a < b ? a : b;
            ek[3 * t + e].b = a < b ? b : a;
            ek[3 * t + e].t = t; ek[3 * t + e].e = e;
        }
    }
    qsort(ek, (size_t)(3 * nt), sizeof(EdgeKey), cmp_ek);
    M->et = (int *)xmalloc((size_t)(3 * nt) * sizeof(int));
    for (int k = 0; k < 3 * nt; k++) M->et[k] = -1;
    for (int k = 0; k + 1 < 3 * nt; ) {
        if (ek[k].a == ek[k + 1].a && ek[k].b == ek[k + 1].b) {
            M->et[3 * ek[k].t + ek[k].e] = ek[k + 1].t;
            M->et[3 * ek[k + 1].t + ek[k + 1].e] = ek[k].t;
            k += 2;
        } else k += 1;
    }
    free(ek);

    /* Aire duale (barycentrique / Voronoï exacte sur triangles) */
    M->area_v = (double *)xcalloc((size_t)nn, sizeof(double));
    for (int t = 0; t < nt; t++) {
        int a = M->tri[3 * t], b = M->tri[3 * t + 1], c = M->tri[3 * t + 2];
        double A = 0.5 * fabs(orient2d(M->x[a], M->y[a], M->x[b], M->y[b], M->x[c], M->y[c]));
        M->area_v[a] += A / 3.0;
        M->area_v[b] += A / 3.0;
        M->area_v[c] += A / 3.0;
    }
}

static void mesh_free(Mesh *M) {
    free(M->x); free(M->y); free(M->bdry); free(M->tri); free(M->et); free(M->area_v);
}

/* ============================================================================
 * 4. CSR + BiCGSTAB préconditionné ILU(0)
 * ============================================================================ */

typedef struct { int n; long nz, cap; int *I, *J; double *V; } Trip;
typedef struct { int n, nz; int *rp, *ci; double *v; } CSR;

static void trip_init(Trip *T, int n) {
    T->n = n; T->nz = 0; T->cap = 10L * n + 64;
    T->I = (int *)xmalloc((size_t)T->cap * sizeof(int));
    T->J = (int *)xmalloc((size_t)T->cap * sizeof(int));
    T->V = (double *)xmalloc((size_t)T->cap * sizeof(double));
}

static void trip_add(Trip *T, int i, int j, double v) {
    if (fabs(v) < 1e-16) return;
    if (T->nz >= T->cap) {
        T->cap *= 2;
        T->I = (int *)xrealloc(T->I, (size_t)T->cap * sizeof(int));
        T->J = (int *)xrealloc(T->J, (size_t)T->cap * sizeof(int));
        T->V = (double *)xrealloc(T->V, (size_t)T->cap * sizeof(double));
    }
    T->I[T->nz] = i; T->J[T->nz] = j; T->V[T->nz] = v; T->nz++;
}

static void trip_free(Trip *T) { free(T->I); free(T->J); free(T->V); }

typedef struct { int c; double v; } CV;
static int cmp_cv(const void *a, const void *b) { return ((const CV *)a)->c - ((const CV *)b)->c; }

static void csr_build(CSR *A, const Trip *T) {
    int n = T->n;
    int *cnt = (int *)xcalloc((size_t)n, sizeof(int));
    for (long k = 0; k < T->nz; k++) cnt[T->I[k]]++;
    int *rp = (int *)xmalloc((size_t)(n + 1) * sizeof(int));
    rp[0] = 0;
    for (int i = 0; i < n; i++) rp[i + 1] = rp[i] + cnt[i];
    CV *buf = (CV *)xmalloc((size_t)T->nz * sizeof(CV));
    int *pos = (int *)xmalloc((size_t)n * sizeof(int));
    memcpy(pos, rp, (size_t)n * sizeof(int));
    for (long k = 0; k < T->nz; k++) {
        int r = T->I[k], p = pos[r]++;
        buf[p].c = T->J[k]; buf[p].v = T->V[k];
    }
    free(cnt); free(pos);

    A->n = n;
    A->rp = (int *)xmalloc((size_t)(n + 1) * sizeof(int));
    A->ci = (int *)xmalloc((size_t)(T->nz + n) * sizeof(int));
    A->v  = (double *)xmalloc((size_t)(T->nz + n) * sizeof(double));
    int nz = 0; A->rp[0] = 0;
    for (int i = 0; i < n; i++) {
        int s = rp[i], len = rp[i + 1] - s;
        if (len > 1) qsort(buf + s, (size_t)len, sizeof(CV), cmp_cv);
        for (int k = 0; k < len; ) {
            int c = buf[s + k].c; double sum = buf[s + k].v; k++;
            while (k < len && buf[s + k].c == c) { sum += buf[s + k].v; k++; }
            if (fabs(sum) > 1e-16 || c == i) { A->ci[nz] = c; A->v[nz] = sum; nz++; }
        }
        A->rp[i + 1] = nz;
    }
    A->nz = nz;
    free(rp); free(buf);
}

static void csr_free(CSR *A) { free(A->rp); free(A->ci); free(A->v); }

static void spmv(const CSR *A, const double *x, double *y) {
    for (int i = 0; i < A->n; i++) {
        double s = 0.0;
        for (int k = A->rp[i]; k < A->rp[i + 1]; k++) s += A->v[k] * x[A->ci[k]];
        y[i] = s;
    }
}

static double dot_v(int n, const double *a, const double *b) {
    double s = 0.0;
    for (int i = 0; i < n; i++) s += a[i] * b[i];
    return s;
}

static int ilu0_factor(const CSR *A, double *LU, int *dg) {
    int n = A->n;
    memcpy(LU, A->v, (size_t)A->nz * sizeof(double));
    for (int i = 0; i < n; i++) {
        dg[i] = -1;
        for (int k = A->rp[i]; k < A->rp[i + 1]; k++) if (A->ci[k] == i) { dg[i] = k; break; }
        if (dg[i] < 0) return -1;
    }
    int *map = (int *)xmalloc((size_t)n * sizeof(int));
    for (int i = 0; i < n; i++) map[i] = -1;
    for (int i = 0; i < n; i++) {
        for (int k = A->rp[i]; k < A->rp[i + 1]; k++) map[A->ci[k]] = k;
        for (int k = A->rp[i]; k < dg[i]; k++) {
            int j = A->ci[k];
            double piv = LU[dg[j]];
            if (fabs(piv) < 1e-14) { free(map); return -1; }
            double mult = LU[k] / piv;
            LU[k] = mult;
            for (int l = dg[j] + 1; l < A->rp[j + 1]; l++) {
                int p = map[A->ci[l]];
                if (p >= 0) LU[p] -= mult * LU[l];
            }
        }
        for (int k = A->rp[i]; k < A->rp[i + 1]; k++) map[A->ci[k]] = -1;
        if (fabs(LU[dg[i]]) < 1e-14) { free(map); return -1; }
    }
    free(map);
    return 0;
}

static void ilu0_apply(const CSR *A, const double *LU, const int *dg, const double *r, double *z) {
    int n = A->n;
    for (int i = 0; i < n; i++) {
        double s = r[i];
        for (int k = A->rp[i]; k < dg[i]; k++) s -= LU[k] * z[A->ci[k]];
        z[i] = s;
    }
    for (int i = n - 1; i >= 0; i--) {
        double s = z[i];
        for (int k = dg[i] + 1; k < A->rp[i + 1]; k++) s -= LU[k] * z[A->ci[k]];
        z[i] = s / LU[dg[i]];
    }
}

static int bicgstab_ilu0(const CSR *A, const double *b, double *x, double tol, int maxit) {
    int n = A->n;
    double *LU = (double *)xmalloc((size_t)A->nz * sizeof(double));
    int    *dg = (int *)xmalloc((size_t)n * sizeof(int));
    int ok_ilu = (ilu0_factor(A, LU, dg) == 0);
    double *r  = (double *)xmalloc((size_t)n * sizeof(double));
    double *rh = (double *)xmalloc((size_t)n * sizeof(double));
    double *p  = (double *)xcalloc((size_t)n, sizeof(double));
    double *v  = (double *)xcalloc((size_t)n, sizeof(double));
    double *s  = (double *)xmalloc((size_t)n * sizeof(double));
    double *t  = (double *)xmalloc((size_t)n * sizeof(double));
    double *ph = (double *)xmalloc((size_t)n * sizeof(double));
    double *sh = (double *)xmalloc((size_t)n * sizeof(double));

    memset(x, 0, (size_t)n * sizeof(double));
    memcpy(r, b, (size_t)n * sizeof(double));
    memcpy(rh, r, (size_t)n * sizeof(double));
    double nb = sqrt(dot_v(n, b, b)); if (nb < 1e-15) nb = 1.0;

    double rho = 1.0, alpha = 1.0, omega = 1.0;
    int it = 0;
    for (it = 1; it <= maxit; it++) {
        double rn = dot_v(n, rh, r);
        if (fabs(rn) < 1e-30) break;
        if (it == 1) memcpy(p, r, (size_t)n * sizeof(double));
        else {
            double beta = (rn / rho) * (alpha / omega);
            for (int i = 0; i < n; i++) p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }
        rho = rn;
        if (ok_ilu) ilu0_apply(A, LU, dg, p, ph); else memcpy(ph, p, (size_t)n * sizeof(double));
        spmv(A, ph, v);
        double rv = dot_v(n, rh, v);
        if (fabs(rv) < 1e-30) break;
        alpha = rho / rv;
        for (int i = 0; i < n; i++) s[i] = r[i] - alpha * v[i];
        if (sqrt(dot_v(n, s, s)) / nb < tol) {
            for (int i = 0; i < n; i++) x[i] += alpha * ph[i];
            break;
        }
        if (ok_ilu) ilu0_apply(A, LU, dg, s, sh); else memcpy(sh, s, (size_t)n * sizeof(double));
        spmv(A, sh, t);
        double tt = dot_v(n, t, t);
        if (tt < 1e-30) break;
        omega = dot_v(n, t, s) / tt;
        for (int i = 0; i < n; i++) { x[i] += alpha * ph[i] + omega * sh[i]; r[i] = s[i] - omega * t[i]; }
        if (sqrt(dot_v(n, r, r)) / nb < tol) break;
    }
    free(LU); free(dg); free(r); free(rh); free(p); free(v); free(s); free(t); free(ph); free(sh);
    return it;
}

/* ============================================================================
 * 5. Assemblage VF TPFA (Voronoï dual par triangle) + erreurs + adaptation Dörfler
 * ============================================================================ */

static void solve_tpfa(const Mesh *M, double *u, double *nnz_row, int *iters, double *t_out) {
    double t0 = now_sec();
    int nn = M->nn;
    Trip T; trip_init(&T, nn);
    double *b = (double *)xcalloc((size_t)nn, sizeof(double));
    double *g = (double *)xmalloc((size_t)nn * sizeof(double));
    for (int i = 0; i < nn; i++) g[i] = M->bdry[i] ? u_ex(M->x[i], M->y[i]) : 0.0;

    /* Assemblage par triangle de Delaunay :
     * Pour chaque arête (i, j) du triangle T, la demi-face de Voronoï a pour
     * longueur |sigma_{ij,T}| = 0.5 * cot(angle opposé k) * d_{ij},
     * normale unitaire n = (x_j - x_i) / d_{ij},
     * et transmissivité TPFA = |sigma_{ij,T}| * (n^T K n) / d_{ij}
     *                        = 0.5 * cot(theta_k) * (n^T K n) !
     * Exactement la formule Voronoï TPFA, sans jamais sortir du domaine non convexe ! */
    for (int t = 0; t < M->nt; t++) {
        int v[3] = { M->tri[3 * t], M->tri[3 * t + 1], M->tri[3 * t + 2] };
        double x[3] = { M->x[v[0]], M->x[v[1]], M->x[v[2]] };
        double y[3] = { M->y[v[0]], M->y[v[1]], M->y[v[2]] };
        double area = 0.5 * fabs(orient2d(x[0], y[0], x[1], y[1], x[2], y[2]));
        double cx = (x[0] + x[1] + x[2]) / 3.0, cy = (y[0] + y[1] + y[2]) / 3.0;
        double K[2][2]; eval_K(cx, K);

        /* Terme source sur la cellule duale */
        double fv = f_src(cx, cy) * (area / 3.0);
        for (int a = 0; a < 3; a++) if (!M->bdry[v[a]]) b[v[a]] += fv;

        for (int e = 0; e < 3; e++) {
            int i = e, j = (e + 1) % 3, k = (e + 2) % 3;
            double ux = x[i] - x[k], uy = y[i] - y[k];
            double vx = x[j] - x[k], vy = y[j] - y[k];
            double dot_k = ux * vx + uy * vy;
            double cot_k = dot_k / (2.0 * area);
            double ex = x[j] - x[i], ey = y[j] - y[i];
            double d2 = ex * ex + ey * ey;
            double nx = ex / sqrt(d2), ny = ey / sqrt(d2);
            double kn = nx * (K[0][0] * nx + K[0][1] * ny) + ny * (K[1][0] * nx + K[1][1] * ny);
            double tr = 0.5 * cot_k * kn;
            int gi = v[i], gj = v[j];
            if (!M->bdry[gi]) {
                trip_add(&T, gi, gi, tr);
                if (M->bdry[gj]) b[gi] += tr * g[gj]; else trip_add(&T, gi, gj, -tr);
            }
            if (!M->bdry[gj]) {
                trip_add(&T, gj, gj, tr);
                if (M->bdry[gi]) b[gj] += tr * g[gi]; else trip_add(&T, gj, gi, -tr);
            }
        }
    }
    for (int i = 0; i < nn; i++) if (M->bdry[i]) { trip_add(&T, i, i, 1.0); b[i] = g[i]; }

    CSR A; csr_build(&A, &T);
    *nnz_row = (double)A.nz / (double)A.n;
    *iters = bicgstab_ilu0(&A, b, u, 1e-12, 2000);
    csr_free(&A); trip_free(&T); free(b); free(g);
    *t_out = now_sec() - t0;
}

typedef struct { double l2, linf, ener; } ErrStats;

static ErrStats compute_errors(const Mesh *M, const double *u) {
    double e2 = 0.0, r2 = 0.0, linf = 0.0;
    for (int i = 0; i < M->nn; i++) {
        double ue = u_ex(M->x[i], M->y[i]);
        double d = fabs(u[i] - ue);
        if (d > linf) linf = d;
        e2 += M->area_v[i] * d * d;
        r2 += M->area_v[i] * ue * ue;
    }
    double en_e = 0.0, en_r = 0.0;
    for (int t = 0; t < M->nt; t++) {
        int a = M->tri[3 * t], b = M->tri[3 * t + 1], c = M->tri[3 * t + 2];
        double x0 = M->x[a], y0 = M->y[a], x1 = M->x[b], y1 = M->y[b], x2 = M->x[c], y2 = M->y[c];
        double D = orient2d(x0, y0, x1, y1, x2, y2), area = 0.5 * fabs(D);
        double gx = ((y1 - y2) * u[a] + (y2 - y0) * u[b] + (y0 - y1) * u[c]) / D;
        double gy = ((x2 - x1) * u[a] + (x0 - x2) * u[b] + (x1 - x0) * u[c]) / D;
        double cx = (x0 + x1 + x2) / 3.0, cy = (y0 + y1 + y2) / 3.0;
        double d[6], K[2][2]; uex_d(cx, cy, d); eval_K(cx, K);
        double ex = gx - d[1], ey = gy - d[2];
        en_e += area * (K[0][0] * ex * ex + 2.0 * K[0][1] * ex * ey + K[1][1] * ey * ey);
        en_r += area * (K[0][0] * d[1] * d[1] + 2.0 * K[0][1] * d[1] * d[2] + K[1][1] * d[2] * d[2]);
    }
    ErrStats E = { sqrt(e2 / (r2 > 1e-30 ? r2 : 1.0)), linf, sqrt(en_e / (en_r > 1e-30 ? en_r : 1.0)) };
    return E;
}

static double estimate_eta(const Mesh *M, const double *u, double *eta_t) {
    int nt = M->nt;
    double *gux = (double *)xmalloc((size_t)nt * sizeof(double));
    double *guy = (double *)xmalloc((size_t)nt * sizeof(double));
    double *K00 = (double *)xmalloc((size_t)nt * sizeof(double));
    double *K01 = (double *)xmalloc((size_t)nt * sizeof(double));
    double *K11 = (double *)xmalloc((size_t)nt * sizeof(double));

    for (int t = 0; t < nt; t++) {
        int a = M->tri[3 * t], b = M->tri[3 * t + 1], c = M->tri[3 * t + 2];
        double D = orient2d(M->x[a], M->y[a], M->x[b], M->y[b], M->x[c], M->y[c]);
        gux[t] = ((M->y[b] - M->y[c]) * u[a] + (M->y[c] - M->y[a]) * u[b] + (M->y[a] - M->y[b]) * u[c]) / D;
        guy[t] = ((M->x[c] - M->x[b]) * u[a] + (M->x[a] - M->x[c]) * u[b] + (M->x[b] - M->x[a]) * u[c]) / D;
        double cx = (M->x[a] + M->x[b] + M->x[c]) / 3.0;
        double K[2][2]; eval_K(cx, K);
        K00[t] = K[0][0]; K01[t] = K[0][1]; K11[t] = K[1][1];
    }

    double tot = 0.0;
    for (int t = 0; t < nt; t++) {
        int v[3] = { M->tri[3 * t], M->tri[3 * t + 1], M->tri[3 * t + 2] };
        double area = 0.5 * fabs(orient2d(M->x[v[0]], M->y[v[0]], M->x[v[1]], M->y[v[1]], M->x[v[2]], M->y[v[2]]));
        double hT = 0.0;
        for (int e = 0; e < 3; e++) {
            double l = hypot(M->x[v[(e + 1) % 3]] - M->x[v[e]], M->y[v[(e + 1) % 3]] - M->y[v[e]]);
            if (l > hT) hT = l;
        }
        double cx = (M->x[v[0]] + M->x[v[1]] + M->x[v[2]]) / 3.0;
        double cy = (M->y[v[0]] + M->y[v[1]] + M->y[v[2]]) / 3.0;
        double fv = f_src(cx, cy);
        double e2 = hT * hT * area * fv * fv;
        for (int e = 0; e < 3; e++) {
            int t2 = M->et[3 * t + e];
            if (t2 < 0) continue;
            int a = v[e], b = v[(e + 1) % 3];
            double ex = M->x[b] - M->x[a], ey = M->y[b] - M->y[a], le = hypot(ex, ey);
            double nx = ey / le, ny = -ex / le;
            double q1 = (K00[t] * gux[t] + K01[t] * guy[t]) * nx + (K01[t] * gux[t] + K11[t] * guy[t]) * ny;
            double q2 = (K00[t2] * gux[t2] + K01[t2] * guy[t2]) * nx + (K01[t2] * gux[t2] + K11[t2] * guy[t2]) * ny;
            e2 += 0.5 * le * le * (q1 - q2) * (q1 - q2);
        }
        eta_t[t] = e2;
        tot += e2;
    }
    free(gux); free(guy); free(K00); free(K01); free(K11);
    return sqrt(tot);
}

typedef struct { int id; double val; } IdVal;
static int cmp_idval_desc(const void *a, const void *b) {
    double da = ((const IdVal *)a)->val, db = ((const IdVal *)b)->val;
    return (db > da) - (db < da);
}

static int refine_dorfler(Pts *P, const Mesh *M, const double *eta_t, double theta) {
    int nt = M->nt;
    IdVal *ord = (IdVal *)xmalloc((size_t)nt * sizeof(IdVal));
    double tot = 0.0;
    for (int t = 0; t < nt; t++) { ord[t].id = t; ord[t].val = eta_t[t]; tot += eta_t[t]; }
    qsort(ord, (size_t)nt, sizeof(IdVal), cmp_idval_desc);

    int *mark = (int *)xcalloc((size_t)M->nn, sizeof(int));
    double acc = 0.0;
    for (int k = 0; k < nt && acc < theta * theta * tot; k++) {
        int t = ord[k].id;
        acc += ord[k].val;
        for (int a = 0; a < 3; a++) mark[M->tri[3 * t + a]] = 1;
    }
    free(ord);

    int added = 0;
    for (int t = 0; t < nt; t++) {
        for (int e = 0; e < 3; e++) {
            int a = M->tri[3 * t + e], b = M->tri[3 * t + (e + 1) % 3];
            if (!mark[a] && !mark[b]) continue;
            double mx = 0.5 * (M->x[a] + M->x[b]), my = 0.5 * (M->y[a] + M->y[b]);
            if (!in_lshape(mx, my, 1e-10)) continue;
            double L = hypot(M->x[b] - M->x[a], M->y[b] - M->y[a]);
            int ok = 1;
            for (int j = 0; j < P->n; j++) {
                if (hypot(P->x[j] - mx, P->y[j] - my) < 0.35 * L) { ok = 0; break; }
            }
            if (ok) { pts_add(P, mx, my); added++; }
        }
    }
    free(mark);
    return added;
}

/* ============================================================================
 * 6. Export Gnuplot, temps réel et GIF (-o output.gif)
 * ============================================================================ */

static void export_mesh(const char *prefix, const Mesh *M, const double *u) {
    char fn[128];
    snprintf(fn, sizeof(fn), "%s_tri.dat", prefix);
    FILE *ft = fopen(fn, "w");
    if (ft) {
        for (int t = 0; t < M->nt; t++) {
            int a = M->tri[3 * t], b = M->tri[3 * t + 1], c = M->tri[3 * t + 2];
            fprintf(ft, "%.6f %.6f\n%.6f %.6f\n%.6f %.6f\n%.6f %.6f\n\n",
                    M->x[a], M->y[a], M->x[b], M->y[b], M->x[c], M->y[c], M->x[a], M->y[a]);
        }
        fclose(ft);
    }
    snprintf(fn, sizeof(fn), "%s_sol.dat", prefix);
    FILE *fs = fopen(fn, "w");
    if (fs) {
        for (int i = 0; i < M->nn; i++) fprintf(fs, "%.6f %.6f %.6e\n", M->x[i], M->y[i], u[i]);
        fclose(fs);
    }
}

static void plot_realtime(FILE *gp, int ex_id, int iter, const Mesh *M,
                          const double *u, ErrStats E, int delay_ms, int is_gif) {
    if (!gp) return;
    export_mesh("live_tpfa", M, u);
    if (!is_gif) fprintf(gp, "clear\n");
    fprintf(gp, "set multiplot layout 1,2 title 'VF TPFA classique (Ex. %d) | Iter %d | Volumes N=%d | e_L2=%.2e'\n",
            ex_id, iter, M->nn, E.l2);

    fprintf(gp, "set title 'Maillage adaptatif (%d volumes)'\n", M->nn);
    fprintf(gp, "set xrange [-0.05:2.05]; set yrange [-0.05:2.05]; set size ratio -1\n");
    fprintf(gp, "set xlabel 'x'; set ylabel 'y'; unset key\n");
    fprintf(gp, "plot 'live_tpfa_tri.dat' with lines lw 0.45 lc rgb '#444444'\n");

    fprintf(gp, "set title 'Solution VF TPFA u_h(x,y)'\n");
    fprintf(gp, "set xrange [0:2]; set yrange [0:2]; set size nosquare\n");
    fprintf(gp, "set view 58, 35; set xlabel 'x'; set ylabel 'y'; set zlabel 'u_h'\n");
    fprintf(gp, "splot 'live_tpfa_sol.dat' using 1:2:3 with points pt 7 ps 0.55 palette\n");

    fprintf(gp, "unset multiplot\n");
    fflush(gp);
    if (!is_gif && delay_ms > 0) SLEEP_MS(delay_ms);
}

static void run_example(int ex_id, const char *title, int max_iter, int nmax,
                        FILE *gp, int delay_ms, int is_gif) {
    EX_ID = ex_id;
    printf("\n==========================================================================================\n");
    printf(" EXEMPLE %d (VF TPFA) : %s\n", ex_id, title);
    printf("==========================================================================================\n");
    printf("%4s %7s %7s %6s %6s %5s %10s %10s %10s %6s\n",
           "Iter", "Nb_Vol", "Nb_Tri", "angmin", "obtus", "iterK", "e_L2", "e_max", "eta", "t(s)");
    printf("------------------------------------------------------------------------------------------\n");

    Pts P; pts_lshape_grid(&P, 6);
    char fhist[64], prefF[64];
    snprintf(fhist, sizeof(fhist), "ex%d_history.dat", ex_id);
    snprintf(prefF, sizeof(prefF), "ex%d_final", ex_id);
    FILE *fh = fopen(fhist, "w");
    if (fh) fprintf(fh, "# iter Nb_Vol Nb_Tri L2 Linf ener eta t emax_noeuds angle_min n_obtus nnz_ligne\n");

    for (int it = 0; it < max_iter; it++) {
        Mesh M; mesh_build(&M, &P);
        double *u = (double *)xcalloc((size_t)M.nn, sizeof(double));
        double *eta_t = (double *)xcalloc((size_t)M.nt, sizeof(double));
        double nnz_r = 0.0, ts = 0.0;
        int it_s = 0;
        solve_tpfa(&M, u, &nnz_r, &it_s, &ts);
        ErrStats E = compute_errors(&M, u);
        double eta = estimate_eta(&M, u, eta_t);
        int stop0 = (M.nn >= nmax || it == max_iter - 1);
        double amin; int nob; mesh_quality(&M, &amin, &nob);
        double emax = point_errors(&M, u, ex_id, 0);
        printf("%4d %7d %7d %6.1f %6d %5d %10.3e %10.3e %10.3e %6.3f\n",
               it, M.nn, M.nt, amin, nob, it_s, E.l2, emax, eta, ts);
        (void)nnz_r; (void)stop0;
        if (fh) fprintf(fh, "%d %d %d %.6e %.6e %.6e %.6e %.4f %.6e %.2f %d %.2f\n",
                        it, M.nn, M.nt, E.l2, E.linf, E.ener, eta, ts, emax, amin, nob, nnz_r);
        plot_realtime(gp, ex_id, it, &M, u, E, delay_ms, is_gif);

        int stop = (M.nn >= nmax || it == max_iter - 1);
        if (stop) {
            export_mesh(prefF, &M, u);
            if (PRINT_PTS) printf("   --- %s, Ex. %d : erreur ponctuelle aux %d noeuds (fichier ex%d_erreur_points.dat) ---\n", "VF TPFA", ex_id, M.nn, ex_id);
            point_errors(&M, u, ex_id, 1);
            printf("   Bilan : %d noeuds, %d triangles, e_L2 = %.3e, e_max(noeuds) = %.3e\n", M.nn, M.nt, E.l2, emax);
        }
        int added = stop ? 0 : refine_dorfler(&P, &M, eta_t, 0.50);
        free(eta_t); free(u); mesh_free(&M);
        if (stop || added == 0) break;
    }
    if (fh) fclose(fh);
    pts_free(&P);
}

int main(int argc, char **argv) {
    int gui = 1, delay_ms = 650, only_ex = -1;
    const char *gif_out = NULL;
    for (int i = 1; i < argc; i++) {
        if (!strncmp(argv[i], "gui=", 4))        gui = atoi(argv[i] + 4);
        else if (!strncmp(argv[i], "delay=", 6)) delay_ms = atoi(argv[i] + 6);
        else if (!strncmp(argv[i], "ex=", 3))    only_ex = atoi(argv[i] + 3);
        else if (!strncmp(argv[i], "pts=", 4))   PRINT_PTS = atoi(argv[i] + 4);
        else if (!strcmp(argv[i], "-o") && i + 1 < argc) gif_out = argv[++i];
    }

    int is_gif = (gif_out != NULL);
    FILE *gp = NULL;
    if (is_gif) {
        gp = POPEN("gnuplot", "w");
        if (gp) {
            int gd = delay_ms / 10; if (gd < 10) gd = 10;
            fprintf(gp, "set terminal gif animate optimize delay %d size 1100,520\n", gd);
            fprintf(gp, "set output '%s'\n", gif_out);
            fflush(gp);
        }
    } else if (gui) {
        gp = POPEN("gnuplot -persistent", "w");
    }

    if (only_ex < 0 || only_ex == 0)
        run_example(0, "TEMOIN ISOTROPE K = I (TPFA consistant)",                        7, 2000, gp, delay_ms, is_gif);
    if (only_ex < 0 || only_ex == 1)
        run_example(1, "Tenseur plein constant K=[[2,-1],[-1,2]] (TPFA non consistant)", 7, 2000, gp, delay_ms, is_gif);
    if (only_ex < 0 || only_ex == 2)
        run_example(2, "Tenseur plein discontinu x=1 (I vs [[8,-7],[-7,20]])",           7, 2000, gp, delay_ms, is_gif);
    if (only_ex < 0 || only_ex == 3)
        run_example(3, "Forte anisotropie 50:1 (30 deg) + pic au coin rentrant (1,1)",   8, 2000, gp, delay_ms, is_gif);

    if (gp) {
        if (is_gif) { fprintf(gp, "set output\n"); fflush(gp); }
        PCLOSE(gp);
        if (is_gif) printf("\nAnimation GIF enregistree dans : %s\n", gif_out);
    }
    return 0;
}
