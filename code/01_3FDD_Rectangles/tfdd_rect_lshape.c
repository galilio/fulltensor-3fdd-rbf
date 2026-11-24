/*
 * ============================================================================
 *  Méthode 1 : Décomposition de domaine à 3 champs (3FDD-RBF) sur RECTANGLES
 *  pour l'équation de diffusion anisotrope à tenseur plein :
 *
 *      -div( K(x) grad u(x) ) = f(x)   dans Omega (domaine NON CONVEXE en L)
 *      u(x) = g_D(x)                   sur dOmega
 *
 *  Domaine non convexe en L initialisé par EXACTEMENT 2 rectangles :
 *      R_1 = [0, 1] x [0, 2]   (barre verticale gauche)
 *      R_2 = [1, 2] x [0, 1]   (carré inférieur droit)
 *      Coin rentrant non convexe en (1, 1) (angle 270 deg)
 *      Interface initiale Gamma_{12} = {1} x [0, 1]
 *
 *  Trois exemples de complexité croissante :
 *    Ex 1 : K = [[2,-1],[-1,2]] constant, solution lisse sin(pi x/2)sin(pi y/2) + 0.25 xy
 *    Ex 2 : K discontinu en x=1 (I vs [[8,-7],[-7,20]]), u = exp(0.5(x+y)) (flux continu)
 *    Ex 3 : K discontinu anisotrope (I vs [[12,-8],[-8,18]]), gradient raide au coin
 *           rentrant (1,1) : u = exp(-35((x-0.92)^2+(y-0.92)^2)) + 0.2(x+y)
 *
 *  Double adaptativité :
 *    - Niveau 1 (nœuds dans chaque rectangle R_k) : ARS par résidu fort EDP
 *    - Niveau 2 (décomposition en rectangles)     : bissection de R_k en 2 rectangles
 *      fils selon son plus grand côté dès que N_{u,k} > nsplit (gestion exacte
 *      des arêtes non-conformes / hanging interfaces sur Sigma)
 *
 *  Compilation :
 *      gcc -O3 -std=c99 -Wall -Wextra -Wpedantic -o tfdd_rect_lshape tfdd_rect_lshape.c -lm
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

#define MAX_SEGS 64

/* ============================================================================
 * 1. Utilitaires mémoire, temps et LU dense (pivot partiel)
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

static double now_sec(void) {
    return (double)clock() / (double)CLOCKS_PER_SEC;
}

static int lu_factor(double *A, int n, int *piv) {
    for (int i = 0; i < n; i++) piv[i] = i;
    for (int k = 0; k < n; k++) {
        int p = k;
        double mx = fabs(A[k * n + k]);
        for (int i = k + 1; i < n; i++) {
            double v = fabs(A[i * n + k]);
            if (v > mx) { mx = v; p = i; }
        }
        if (mx < 1e-14) return -1;
        if (p != k) {
            int ti = piv[k]; piv[k] = piv[p]; piv[p] = ti;
            for (int j = 0; j < n; j++) {
                double t = A[k * n + j];
                A[k * n + j] = A[p * n + j];
                A[p * n + j] = t;
            }
        }
        double inv = 1.0 / A[k * n + k];
        for (int i = k + 1; i < n; i++) {
            double m = A[i * n + k] * inv;
            A[i * n + k] = m;
            for (int j = k + 1; j < n; j++) A[i * n + j] -= m * A[k * n + j];
        }
    }
    return 0;
}

static void lu_solve(const double *LU, int n, const int *piv, double *b) {
    double *x = (double *)xmalloc((size_t)n * sizeof(double));
    for (int i = 0; i < n; i++) x[i] = b[piv[i]];
    for (int i = 0; i < n; i++) {
        double s = x[i];
        for (int j = 0; j < i; j++) s -= LU[i * n + j] * x[j];
        x[i] = s;
    }
    for (int i = n - 1; i >= 0; i--) {
        double s = x[i];
        for (int j = i + 1; j < n; j++) s -= LU[i * n + j] * x[j];
        x[i] = s / LU[i * n + i];
    }
    memcpy(b, x, (size_t)n * sizeof(double));
    free(x);
}

/* ============================================================================
 * 2. Trois exemples de complexité croissante sur le domaine non convexe en L
 * ============================================================================ */

static int EX_ID = 1;

static void eval_K(double x, double K[2][2]) {
    if (EX_ID == 1) {
        K[0][0] =  2.0; K[0][1] = -1.0;
        K[1][0] = -1.0; K[1][1] =  2.0;
    } else if (EX_ID == 2) {
        if (x < 1.0) {
            K[0][0] = 1.0; K[0][1] = 0.0;
            K[1][0] = 0.0; K[1][1] = 1.0;
        } else {
            K[0][0] =  8.0; K[0][1] = -7.0;
            K[1][0] = -7.0; K[1][1] = 20.0;
        }
    } else {
        /* Exemple 3 : forte anisotropie tournée (ratio 50:1, angle 30 deg)
         * K = R(30 deg) diag(50, 1) R(30 deg)^T = [[37.75, 21.2176], [21.2176, 13.25]] */
        double th = 30.0 * (M_PI / 180.0);
        double c = cos(th), s = sin(th);
        double l1 = 50.0, l2 = 1.0;
        K[0][0] = l1 * c * c + l2 * s * s;
        K[0][1] = (l1 - l2) * c * s;
        K[1][0] = K[0][1];
        K[1][1] = l1 * s * s + l2 * c * c;
    }
}

/* d[0]=u, d[1]=ux, d[2]=uy, d[3]=uxx, d[4]=uxy, d[5]=uyy */
static void uex_d(double x, double y, double d[6]) {
    if (EX_ID == 1) {
        double k = 0.5 * M_PI;
        double sx = sin(k * x), cx = cos(k * x);
        double sy = sin(k * y), cy = cos(k * y);
        d[0] = sx * sy + 0.25 * x * y;
        d[1] = k * cx * sy + 0.25 * y;
        d[2] = k * sx * cy + 0.25 * x;
        d[3] = -k * k * sx * sy;
        d[4] =  k * k * cx * cy + 0.25;
        d[5] = -k * k * sx * sy;
    } else if (EX_ID == 2) {
        double a = 0.5;
        double e = exp(a * (x + y));
        d[0] = e;
        d[1] = a * e;
        d[2] = a * e;
        d[3] = a * a * e;
        d[4] = a * a * e;
        d[5] = a * a * e;
    } else {
        /* Pic raide centré juste près du coin rentrant non convexe (1, 1) */
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
 * 3. Quadratures (Gauss-Legendre 5 pts 1D et 5x5 2D sur rectangle)
 * ============================================================================ */

static const double GL5_X[5] = {
    -0.9061798459386640, -0.5384693101056831, 0.0,
     0.5384693101056831,  0.9061798459386640
};
static const double GL5_W[5] = {
    0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
    0.4786286704993665, 0.2369268850561891
};

static double leg1d(int m, double s) {
    double t = 2.0 * s - 1.0;
    if (m == 0) return 1.0;
    if (m == 1) return t;
    if (m == 2) return 0.5 * (3.0 * t * t - 1.0);
    if (m == 3) return 0.5 * (5.0 * t * t * t - 3.0 * t);
    if (m == 4) return 0.125 * (35.0 * t * t * t * t - 30.0 * t * t + 3.0);
    return 0.125 * (63.0 * pow(t, 5) - 70.0 * t * t * t + 15.0 * t);
}

/* ============================================================================
 * 4. Structures des sous-domaines rectangulaires et interfaces non-conformes
 * ============================================================================ */

typedef struct {
    double x0, y0, x1, y1;    /* extrémités du sous-segment de bord */
    int    iface;             /* indice global sur Sigma, ou -1 si bord externe dOmega */
    int    iorient;           /* +1 ou -1 selon l'orientation globale de l'interface */
    int    nb_e;              /* nombre de modes Legendre pour lambda_k sur ce segment */
} BndSeg;

typedef struct {
    double a, b, c, d;        /* rectangle [a, b] x [c, d] */
    double cx, cy, lx, ly;    /* centre et largeurs */
    int    lev;               /* niveau de bissection */

    /* Sous-segments de bord (gère les arêtes non-conformes / hanging nodes !) */
    int    nseg;
    BndSeg seg[MAX_SEGS];

    /* Nuage local de centres RBF dans [a, b] x [c, d] */
    int     nu, cap_u;
    double *ux, *uy;

    /* Base cardinale locale PHS r^5 + P2 (dim = nu + 6) */
    double *G_lu;
    int    *G_piv;

    /* Solution 3FDD locale */
    int     nl;
    double *u_nodal;
    double *lambda;
    double *c_phs;
    double  res_max;
    double  alpha;            /* parametre alpha_W adaptatif local (Wendland) */
} RectSub;

typedef struct {
    int    k1, k2;
    double x0, y0, x1, y1;
    double len;
    int    nc_e;
    int    dof0;
} IfaceSeg;

typedef struct {
    int       ns, cap_s;
    RectSub  *sub;
    int       nif;
    IfaceSeg *ifs;
    int       nb, nc;         /* nb > nc : condition inf-sup garantie */
    int       n_sigma;
    double   *psi;
} DomainRect3FDD;

/* Découpe chaque côté des M rectangles aux sommets de tous les rectangles adjacents
 * et apparie exactement toutes les interfaces partagées Gamma_{ij} sur Sigma */
static void build_rect_interfaces(DomainRect3FDD *D) {
    const double eps = 1e-9;

    for (int k = 0; k < D->ns; k++) {
        RectSub *R = &D->sub[k];
        R->cx = 0.5 * (R->a + R->b); R->cy = 0.5 * (R->c + R->d);
        R->lx = R->b - R->a;         R->ly = R->d - R->c;
        R->nseg = 0;

        /* 4 côtés : 0=bas(y=c), 1=droite(x=b), 2=haut(y=d), 3=gauche(x=a) */
        for (int side = 0; side < 4; side++) {
            double pts[64];
            int np = 2;
            int horiz = (side == 0 || side == 2);
            double tmin = horiz ? R->a : R->c;
            double tmax = horiz ? R->b : R->d;
            double fixed = (side == 0) ? R->c : (side == 1 ? R->b : (side == 2 ? R->d : R->a));
            pts[0] = tmin; pts[1] = tmax;

            /* Collecter tous les points de coupure induits par des rectangles voisins */
            for (int j = 0; j < D->ns; j++) {
                if (j == k) continue;
                const RectSub *Q = &D->sub[j];
                if (horiz) {
                    if (fabs(Q->c - fixed) < eps || fabs(Q->d - fixed) < eps) {
                        double cands[2] = { Q->a, Q->b };
                        for (int m = 0; m < 2; m++) {
                            if (cands[m] > tmin + eps && cands[m] < tmax - eps) {
                                int dup = 0;
                                for (int q = 0; q < np; q++) if (fabs(pts[q] - cands[m]) < eps) dup = 1;
                                if (!dup && np < 64) pts[np++] = cands[m];
                            }
                        }
                    }
                } else {
                    if (fabs(Q->a - fixed) < eps || fabs(Q->b - fixed) < eps) {
                        double cands[2] = { Q->c, Q->d };
                        for (int m = 0; m < 2; m++) {
                            if (cands[m] > tmin + eps && cands[m] < tmax - eps) {
                                int dup = 0;
                                for (int q = 0; q < np; q++) if (fabs(pts[q] - cands[m]) < eps) dup = 1;
                                if (!dup && np < 64) pts[np++] = cands[m];
                            }
                        }
                    }
                }
            }
            /* Tri croissant de pts */
            for (int i = 0; i < np; i++) {
                for (int j = i + 1; j < np; j++) {
                    if (pts[j] < pts[i]) { double tmp = pts[i]; pts[i] = pts[j]; pts[j] = tmp; }
                }
            }
            for (int i = 0; i < np - 1 && R->nseg < MAX_SEGS; i++) {
                BndSeg *sg = &R->seg[R->nseg++];
                if (horiz) {
                    sg->x0 = pts[i]; sg->y0 = fixed; sg->x1 = pts[i + 1]; sg->y1 = fixed;
                } else {
                    sg->x0 = fixed; sg->y0 = pts[i]; sg->x1 = fixed; sg->y1 = pts[i + 1];
                }
                sg->iface = -1;
                sg->iorient = 1;
                double slen = hypot(sg->x1 - sg->x0, sg->y1 - sg->y0);
                double ref_l = horiz ? R->lx : R->ly;
                sg->nb_e = (slen < 0.25 * ref_l) ? 2 : D->nb;
            }
        }
    }

    /* Appariement 1-à-1 des sous-segments internes partagés entre deux rectangles k1 < k2 */
    free(D->ifs);
    int cap_if = 8 * D->ns + 16;
    D->ifs = (IfaceSeg *)xmalloc((size_t)cap_if * sizeof(IfaceSeg));
    D->nif = 0;
    int dof_cur = 0;

    for (int k1 = 0; k1 < D->ns; k1++) {
        RectSub *R1 = &D->sub[k1];
        for (int e1 = 0; e1 < R1->nseg; e1++) {
            BndSeg *s1 = &R1->seg[e1];
            if (s1->iface >= 0) continue;
            for (int k2 = k1 + 1; k2 < D->ns; k2++) {
                RectSub *R2 = &D->sub[k2];
                for (int e2 = 0; e2 < R2->nseg; e2++) {
                    BndSeg *s2 = &R2->seg[e2];
                    if (s2->iface >= 0) continue;
                    int same = (hypot(s1->x0 - s2->x0, s1->y0 - s2->y0) < 1e-7 &&
                                hypot(s1->x1 - s2->x1, s1->y1 - s2->y1) < 1e-7);
                    int rev  = (hypot(s1->x0 - s2->x1, s1->y0 - s2->y1) < 1e-7 &&
                                hypot(s1->x1 - s2->x0, s1->y1 - s2->y0) < 1e-7);
                    if (same || rev) {
                        if (D->nif >= cap_if) {
                            cap_if *= 2;
                            D->ifs = (IfaceSeg *)xrealloc(D->ifs, (size_t)cap_if * sizeof(IfaceSeg));
                        }
                        int id = D->nif++;
                        int min_nb = s1->nb_e < s2->nb_e ? s1->nb_e : s2->nb_e;
                        int nc_e = D->nc < (min_nb - 1) ? D->nc : (min_nb - 1);
                        if (nc_e < 1) nc_e = 1;
                        D->ifs[id].k1 = k1; D->ifs[id].k2 = k2;
                        D->ifs[id].x0 = s1->x0; D->ifs[id].y0 = s1->y0;
                        D->ifs[id].x1 = s1->x1; D->ifs[id].y1 = s1->y1;
                        D->ifs[id].len = hypot(s1->x1 - s1->x0, s1->y1 - s1->y0);
                        D->ifs[id].nc_e = nc_e;
                        D->ifs[id].dof0 = dof_cur;
                        dof_cur += nc_e;
                        s1->iface = id; s1->iorient = 1;
                        s2->iface = id; s2->iorient = same ? 1 : -1;
                        break;
                    }
                }
            }
        }
    }
    D->n_sigma = dof_cur;
}

/* ============================================================================
 * 5. Points locaux dans un rectangle R_k et base cardinale PHS r^5 + P2
 * ============================================================================ */

static void rect_free_arrays(RectSub *R) {
    free(R->ux); free(R->uy); free(R->G_lu); free(R->G_piv);
    free(R->u_nodal); free(R->lambda); free(R->c_phs);
    R->ux = R->uy = R->G_lu = R->u_nodal = R->lambda = R->c_phs = NULL;
    R->G_piv = NULL;
    R->nu = R->cap_u = 0;
}

static void rect_add_pt(RectSub *R, double x, double y, double min_dist) {
    for (int i = 0; i < R->nu; i++) {
        if (hypot(R->ux[i] - x, R->uy[i] - y) < min_dist) return;
    }
    if (R->nu >= R->cap_u) {
        R->cap_u = R->cap_u ? 2 * R->cap_u : 64;
        R->ux = (double *)xrealloc(R->ux, (size_t)R->cap_u * sizeof(double));
        R->uy = (double *)xrealloc(R->uy, (size_t)R->cap_u * sizeof(double));
    }
    R->ux[R->nu] = x; R->uy[R->nu] = y; R->nu++;
}

/* Ensemencement structuré (nx+1) x (ny+1) sur le rectangle R_k, plus les extrémités
 * des sous-segments de bord pour assurer la stabilité inf-sup même sur arêtes non-conformes */
static void rect_seed_points(RectSub *R, int nx, int ny) {
    rect_free_arrays(R);
    double dmin = 0.008 * fmin(R->lx, R->ly);
    for (int i = 0; i <= nx; i++) {
        for (int j = 0; j <= ny; j++) {
            double x = R->a + R->lx * (double)i / (double)nx;
            double y = R->c + R->ly * (double)j / (double)ny;
            rect_add_pt(R, x, y, dmin);
        }
    }
    /* Garantir des points de bord sur chaque sous-segment non-conforme (inf-sup : nb_e+1 intervalles),
     * sauf s'ils tombent pres d'un point de la grille n x n (distance < h/2) */
    double hg = fmin(R->lx / (double)(nx > 0 ? nx : 1), R->ly / (double)(ny > 0 ? ny : 1));
    for (int e = 0; e < R->nseg; e++) {
        const BndSeg *sg = &R->seg[e];
        int m = sg->nb_e + 1;
        /* la grille porte deja >= m intervalles sur ce segment : on n'ajoute que loin de la grille */
        double seglen = hypot(sg->x1 - sg->x0, sg->y1 - sg->y0);
        double dseg = (seglen / hg >= (double)m - 1e-9) ? 0.5 * hg : dmin;
        for (int i = 0; i <= m; i++) {
            double t = (double)i / (double)m;
            rect_add_pt(R, sg->x0 + t * (sg->x1 - sg->x0),
                           sg->y0 + t * (sg->y1 - sg->y0), dseg);
        }
    }
}

/* Base RBF locale (toutes augmentees de P2) :
 *   0 = PHS r^5 (defaut)   1 = Wendland C4 support volume   2 = PHS r^3   3 = TPS r^4 log r
 *   4 = gaussienne exp(-r^2/s^2)   5 = multiquadrique sqrt(r^2+s^2)   6 = MQ inverse 1/sqrt(r^2+s^2)
 * Pour 4-6 : s = KAPPA * h_k, h_k = 1/(sqrt(N_{u,k})-1) en coordonnees normalisees (cf. 3FDD_ARS_tenseur). */
static int RBF_TYPE = 0;
static double KAPPA = 1.0;
static int N_INIT = 6;
static int NB_LEG = 4;          /* modes de Legendre de lambda_k par arete (nc = NB_LEG - 1) */          /* n x n points sur chaque sous-domaine (initial et apres division) */
static const char *RBF_NAME[7] = { "PHS r^5", "Wendland C4", "PHS r^3", "TPS r^4 log r", "Gaussienne", "MQ", "MQ inverse" };

/* alpha_W courant (positionne par sous-domaine avant chaque evaluation de base) */
static double ALPHA_W0 = 18.0;          /* valeur initiale alpha_W^(0) */
static double CUR_ALPHA = 18.0;
static int    ALPHA_ADAPT = 0;          /* 1 : alpha_k adaptatif par residu local (option alpha=adapt) */
static const double ALPHA_MIN = 8.0, ALPHA_MAX = 40.0, ALPHA_GAMMA = 1.3;
static void set_alpha(const RectSub *R) { CUR_ALPHA = (R->alpha > 0.0) ? R->alpha : ALPHA_W0; }

static double rect_wendland_delta(int nu) {
    double d = CUR_ALPHA / sqrt((double)(nu > 1 ? nu : 1));
    if (d < 2.2) d = 2.2;
    if (d > 5.5) d = 5.5;
    return d;
}

static void wendland_rect_eval(double dx, double dy, int nu,
                               double *f, double *fx, double *fy,
                               double *fxx, double *fxy, double *fyy) {
    /* Jauge de Minkowski lissée p=4 symétrisée sur le rectangle normalisé [-0.5, 0.5]^2 :
     * 4 faces : (±1, 0) à distance h=0.5, (0, ±1) à distance h=0.5 => a_e = 2, b_e = 2 */
    double delta = rect_wendland_delta(nu);
    double tx = 2.0 * dx, ty = 2.0 * dy;
    double tx2 = tx * tx, ty2 = ty * ty;
    double Q = tx2 * tx2 + ty2 * ty2; /* (1/2) * [2*tx^4 + 2*ty^4] */
    if (Q < 1e-30) {
        *f = 3.0;
        if (fx) { *fx = 0.0; *fy = 0.0; }
        if (fxx) { *fxx = -168.0 * 4.0 / (delta * delta); *fxy = 0.0; *fyy = -168.0 * 4.0 / (delta * delta); }
        return;
    }
    double rho = pow(Q, 0.25) / delta;
    if (rho >= 1.0) {
        *f = 0.0;
        if (fx) { *fx = 0.0; *fy = 0.0; }
        if (fxx) { *fxx = 0.0; *fxy = 0.0; *fyy = 0.0; }
        return;
    }
    double u = 1.0 - rho, u2 = u * u, u4 = u2 * u2;
    *f = u4 * u2 * (35.0 * rho * rho + 18.0 * rho + 3.0);
    if (!fx && !fxx) return;
    double Wp  = -56.0 * rho * u4 * u * (5.0 * rho + 1.0);
    double Qx  = 8.0 * tx2 * tx;
    double Qy  = 8.0 * ty2 * ty;
    double c1  = rho / (4.0 * Q);
    double rx  = c1 * Qx;
    double ry  = c1 * Qy;
    if (fx) { *fx = Wp * rx; *fy = Wp * ry; }
    if (fxx) {
        double Wpp = 56.0 * u4 * (35.0 * rho * rho - 4.0 * rho - 1.0);
        double Qxx = 48.0 * tx2;
        double Qyy = 48.0 * ty2;
        double rxx = c1 * (Qxx - 0.75 * Qx * Qx / Q);
        double rxy = c1 * (-0.75 * Qx * Qy / Q);
        double ryy = c1 * (Qyy - 0.75 * Qy * Qy / Q);
        *fxx = Wpp * rx * rx + Wp * rxx;
        *fxy = Wpp * rx * ry + Wp * rxy;
        *fyy = Wpp * ry * ry + Wp * ryy;
    }
}

static void phs5_eval(double dx, double dy, double *f, double *fx, double *fy,
                      double *fxx, double *fxy, double *fyy) {
    double r2 = dx * dx + dy * dy;
    double r  = sqrt(r2);
    double r3 = r2 * r;
    *f = r3 * r2;
    if (fx) { *fx = 5.0 * dx * r3; *fy = 5.0 * dy * r3; }
    if (fxx) {
        *fxx = 5.0 * r3 + 15.0 * dx * dx * r;
        *fxy = 15.0 * dx * dy * r;
        *fyy = 5.0 * r3 + 15.0 * dy * dy * r;
    }
}

static void other_rbf_eval(double dx, double dy, int nu,
                           double *f, double *fx, double *fy,
                           double *fxx, double *fxy, double *fyy) {
    double r2 = dx * dx + dy * dy;
    double hk = 1.0 / fmax(sqrt((double)nu) - 1.0, 1.0), sg = KAPPA * hk, s2 = sg * sg;
    double F, Fx = 0, Fy = 0, Fxx = 0, Fxy = 0, Fyy = 0;
    switch (RBF_TYPE) {
    case 2: {                                   /* r^3 */
        double r = sqrt(r2);
        F = r2 * r; Fx = 3.0 * dx * r; Fy = 3.0 * dy * r;
        if (r > 1e-14) { Fxx = 3.0 * r + 3.0 * dx * dx / r; Fxy = 3.0 * dx * dy / r; Fyy = 3.0 * r + 3.0 * dy * dy / r; }
        break; }
    case 3: {                                   /* r^4 log r = 1/2 rho^2 log rho, rho = r^2 */
        if (r2 < 1e-30) { F = 0; break; }
        double lg = 2.0 * log(r2) + 1.0;
        F = 0.5 * r2 * r2 * log(r2);
        Fx = dx * r2 * lg; Fy = dy * r2 * lg;
        Fxx = r2 * lg + 2.0 * dx * dx * lg + 4.0 * dx * dx;
        Fyy = r2 * lg + 2.0 * dy * dy * lg + 4.0 * dy * dy;
        Fxy = 2.0 * dx * dy * lg + 4.0 * dx * dy;
        break; }
    case 4: {                                   /* gaussienne */
        double e = exp(-r2 / s2), a = 2.0 / s2;
        F = e; Fx = -a * dx * e; Fy = -a * dy * e;
        Fxx = (a * a * dx * dx - a) * e; Fyy = (a * a * dy * dy - a) * e; Fxy = a * a * dx * dy * e;
        break; }
    case 5: {                                   /* MQ */
        double p = sqrt(r2 + s2), p3 = p * p * p;
        F = p; Fx = dx / p; Fy = dy / p;
        Fxx = 1.0 / p - dx * dx / p3; Fyy = 1.0 / p - dy * dy / p3; Fxy = -dx * dy / p3;
        break; }
    default: {                                  /* IMQ */
        double q = r2 + s2, p = 1.0 / sqrt(q), p3 = p / q, p5 = p3 / q;
        F = p; Fx = -dx * p3; Fy = -dy * p3;
        Fxx = -p3 + 3.0 * dx * dx * p5; Fyy = -p3 + 3.0 * dy * dy * p5; Fxy = 3.0 * dx * dy * p5;
        break; }
    }
    *f = F;
    if (fx) { *fx = Fx; *fy = Fy; }
    if (fxx) { *fxx = Fxx; *fxy = Fxy; *fyy = Fyy; }
}

static void rbf_rect_eval(double dx, double dy, int nu,
                          double *f, double *fx, double *fy,
                          double *fxx, double *fxy, double *fyy) {
    if (RBF_TYPE == 1)      wendland_rect_eval(dx, dy, nu, f, fx, fy, fxx, fxy, fyy);
    else if (RBF_TYPE == 0) phs5_eval(dx, dy, f, fx, fy, fxx, fxy, fyy);
    else                    other_rbf_eval(dx, dy, nu, f, fx, fy, fxx, fxy, fyy);
}

static void rect_build_basis(RectSub *R) {
    set_alpha(R);
    int n = R->nu, dim = n + 6;
    free(R->G_lu); free(R->G_piv);
    R->G_lu  = (double *)xcalloc((size_t)(dim * dim), sizeof(double));
    R->G_piv = (int *)xmalloc((size_t)dim * sizeof(int));
    double ix = 1.0 / R->lx, iy = 1.0 / R->ly;

    for (int i = 0; i < n; i++) {
        double xi = (R->ux[i] - R->cx) * ix, yi = (R->uy[i] - R->cy) * iy;
        for (int j = 0; j < n; j++) {
            double xj = (R->ux[j] - R->cx) * ix, yj = (R->uy[j] - R->cy) * iy;
            double f;
            rbf_rect_eval(xi - xj, yi - yj, n, &f, NULL, NULL, NULL, NULL, NULL);
            R->G_lu[i * dim + j] = f;
        }
        double p[6] = { 1.0, xi, yi, xi * xi, xi * yi, yi * yi };
        for (int m = 0; m < 6; m++) {
            R->G_lu[i * dim + (n + m)] = p[m];
            R->G_lu[(n + m) * dim + i] = p[m];
        }
    }
    if (lu_factor(R->G_lu, dim, R->G_piv) != 0) {
        fprintf(stderr, "Avertissement : base locale PHS+P2 mal conditionnee (nu=%d)\n", n);
    }
}

static void rect_eval_cardinal(const RectSub *R, double x, double y,
                               double *phi, double *gx, double *gy) {
    set_alpha(R);
    int n = R->nu, dim = n + 6;
    double ix = 1.0 / R->lx, iy = 1.0 / R->ly;
    double xi = (x - R->cx) * ix, yi = (y - R->cy) * iy;
    double *r0 = (double *)xmalloc((size_t)dim * sizeof(double));
    double *rx = gx ? (double *)xmalloc((size_t)dim * sizeof(double)) : NULL;
    double *ry = gy ? (double *)xmalloc((size_t)dim * sizeof(double)) : NULL;

    for (int j = 0; j < n; j++) {
        double xj = (R->ux[j] - R->cx) * ix, yj = (R->uy[j] - R->cy) * iy;
        if (gx) rbf_rect_eval(xi - xj, yi - yj, n, &r0[j], &rx[j], &ry[j], NULL, NULL, NULL);
        else    rbf_rect_eval(xi - xj, yi - yj, n, &r0[j], NULL, NULL, NULL, NULL, NULL);
    }
    double p[6]  = { 1.0, xi, yi, xi * xi, xi * yi, yi * yi };
    double px[6] = { 0.0, 1.0, 0.0, 2.0 * xi, yi, 0.0 };
    double py[6] = { 0.0, 0.0, 1.0, 0.0, xi, 2.0 * yi };
    for (int m = 0; m < 6; m++) {
        r0[n + m] = p[m];
        if (gx) { rx[n + m] = px[m]; ry[n + m] = py[m]; }
    }
    lu_solve(R->G_lu, dim, R->G_piv, r0);
    if (phi) memcpy(phi, r0, (size_t)n * sizeof(double));
    if (gx) {
        lu_solve(R->G_lu, dim, R->G_piv, rx);
        lu_solve(R->G_lu, dim, R->G_piv, ry);
        for (int j = 0; j < n; j++) { gx[j] = rx[j] * ix; gy[j] = ry[j] * iy; }
    }
    free(r0); free(rx); free(ry);
}

static double rect_eval_u_and_pde(const RectSub *R, double x, double y,
                                  double *ux, double *uy, double *L_u) {
    set_alpha(R);
    int n = R->nu;
    double ix = 1.0 / R->lx, iy = 1.0 / R->ly;
    double xi = (x - R->cx) * ix, yi = (y - R->cy) * iy;
    double u = 0.0, gx = 0.0, gy = 0.0, hxx = 0.0, hxy = 0.0, hyy = 0.0;

    for (int j = 0; j < n; j++) {
        double c = R->c_phs[j];
        double xj = (R->ux[j] - R->cx) * ix, yj = (R->uy[j] - R->cy) * iy;
        double f, fx, fy, fxx, fxy, fyy;
        rbf_rect_eval(xi - xj, yi - yj, n, &f, &fx, &fy, &fxx, &fxy, &fyy);
        u   += c * f;
        gx  += c * fx;  gy  += c * fy;
        hxx += c * fxx; hxy += c * fxy; hyy += c * fyy;
    }
    const double *cp = R->c_phs + n;
    u   += cp[0] + cp[1] * xi + cp[2] * yi + cp[3] * xi * xi + cp[4] * xi * yi + cp[5] * yi * yi;
    gx  += cp[1] + 2.0 * cp[3] * xi + cp[4] * yi;
    gy  += cp[2] + cp[4] * xi + 2.0 * cp[5] * yi;
    hxx += 2.0 * cp[3];
    hxy += cp[4];
    hyy += 2.0 * cp[5];

    gx *= ix; gy *= iy;
    hxx *= ix * ix; hxy *= ix * iy; hyy *= iy * iy;
    if (ux) { *ux = gx; *uy = gy; }
    if (L_u) {
        double K[2][2]; eval_K(R->cx, K);
        *L_u = -(K[0][0] * hxx + 2.0 * K[0][1] * hxy + K[1][1] * hyy);
    }
    return u;
}

/* ============================================================================
 * 6. Assemblage 3FDD symétrique et complément de Schur global SPD sur Sigma
 * ============================================================================ */

typedef struct {
    double *K_lu;
    int    *K_piv;
    int     n_loc;
    double *f_rhs;
    double *C_mat;
    int     n_sig_k;
    int    *sig_map;
} LocalBlock;

static void assemble_local(DomainRect3FDD *D, int k, LocalBlock *L) {
    RectSub *R = &D->sub[k];
    rect_build_basis(R);
    int nu = R->nu;
    int nl = 0;
    for (int e = 0; e < R->nseg; e++) nl += R->seg[e].nb_e;
    R->nl = nl;
    int nloc = nu + nl;

    double *Kmat = (double *)xcalloc((size_t)(nloc * nloc), sizeof(double));
    double *rhs  = (double *)xcalloc((size_t)nloc, sizeof(double));
    double *phi  = (double *)xmalloc((size_t)nu * sizeof(double));
    double *gx   = (double *)xmalloc((size_t)nu * sizeof(double));
    double *gy   = (double *)xmalloc((size_t)nu * sizeof(double));

    /* 1. Quadrature 2D sur le rectangle R_k : 2x2 sous-rectangles x (5x5 Gauss-Legendre) */
    double Ktens[2][2];
    eval_K(R->cx, Ktens);
    int nsub_q = (int)ceil(sqrt((double)nu) / 3.0);   /* quadrature 2D : ~3 points de Gauss par centre et par direction */
    if (nsub_q < 2) nsub_q = 2;
    double sub_lx = R->lx / (double)nsub_q, sub_ly = R->ly / (double)nsub_q;

    for (int sx = 0; sx < nsub_q; sx++) {
        for (int sy = 0; sy < nsub_q; sy++) {
            double x0 = R->a + sx * sub_lx, y0 = R->c + sy * sub_ly;
            for (int qx_i = 0; qx_i < 5; qx_i++) {
                double qx = x0 + 0.5 * sub_lx * (1.0 + GL5_X[qx_i]);
                double wx = 0.5 * sub_lx * GL5_W[qx_i];
                for (int qy_i = 0; qy_i < 5; qy_i++) {
                    double qy = y0 + 0.5 * sub_ly * (1.0 + GL5_X[qy_i]);
                    double w  = wx * 0.5 * sub_ly * GL5_W[qy_i];
                    rect_eval_cardinal(R, qx, qy, phi, gx, gy);
                    double fv = f_src(qx, qy);
                    for (int i = 0; i < nu; i++) {
                        rhs[i] += w * fv * phi[i];
                        double kgi_x = Ktens[0][0] * gx[i] + Ktens[0][1] * gy[i];
                        double kgi_y = Ktens[1][0] * gx[i] + Ktens[1][1] * gy[i];
                        for (int j = i; j < nu; j++) {
                            double val = w * (kgi_x * gx[j] + kgi_y * gy[j]);
                            Kmat[i * nloc + j] += val;
                            if (j != i) Kmat[j * nloc + i] += val;
                        }
                    }
                }
            }
        }
    }

    /* 2. Blocs B_k, g_k et C_k sur les sous-segments de bord */
    int n_sig_k = 0;
    for (int e = 0; e < R->nseg; e++) {
        if (R->seg[e].iface >= 0) n_sig_k += D->ifs[R->seg[e].iface].nc_e;
    }
    int *sig_map = (int *)xmalloc((size_t)(n_sig_k > 0 ? n_sig_k : 1) * sizeof(int));
    double *Cmat = (double *)xcalloc((size_t)(nl * (n_sig_k > 0 ? n_sig_k : 1)), sizeof(double));
    int row_off = 0, col_off = 0;

    for (int e = 0; e < R->nseg; e++) {
        const BndSeg *sg = &R->seg[e];
        int nb_e = sg->nb_e;
        int if_id = sg->iface;
        int nc_e = (if_id >= 0) ? D->ifs[if_id].nc_e : 0;
        if (if_id >= 0) {
            for (int m = 0; m < nc_e; m++) sig_map[col_off + m] = D->ifs[if_id].dof0 + m;
        }
        double elen = hypot(sg->x1 - sg->x0, sg->y1 - sg->y0);
        for (int sub = 0; sub < 2; sub++) {
            double s0 = 0.5 * sub, s1 = 0.5 * (sub + 1);
            for (int q = 0; q < 5; q++) {
                double s = 0.5 * (s0 + s1) + 0.5 * (s1 - s0) * GL5_X[q];
                double w = 0.5 * (0.5 * elen) * GL5_W[q];
                double qx = sg->x0 + s * (sg->x1 - sg->x0);
                double qy = sg->y0 + s * (sg->y1 - sg->y0);
                rect_eval_cardinal(R, qx, qy, phi, NULL, NULL);
                for (int l = 0; l < nb_e; l++) {
                    int row_l = row_off + l;
                    double mu = leg1d(l, s);
                    for (int i = 0; i < nu; i++) {
                        double v = -w * phi[i] * mu;
                        Kmat[i * nloc + (nu + row_l)] += v;
                        Kmat[(nu + row_l) * nloc + i] += v;
                    }
                    if (if_id < 0) {
                        rhs[nu + row_l] -= w * u_ex(qx, qy) * mu;
                    } else {
                        double s_if = (sg->iorient > 0) ? s : (1.0 - s);
                        for (int m = 0; m < nc_e; m++) {
                            Cmat[row_l * n_sig_k + (col_off + m)] += w * mu * leg1d(m, s_if);
                        }
                    }
                }
            }
        }
        row_off += nb_e;
        if (if_id >= 0) col_off += nc_e;
    }

    /* 3. Factorisation locale et assemblage du complément de Schur global SPD */
    int *piv = (int *)xmalloc((size_t)nloc * sizeof(int));
    if (lu_factor(Kmat, nloc, piv) != 0) {
        fprintf(stderr, "Erreur : bloc local 3FDD singulier sur rectangle %d (nu=%d, nl=%d)\n", k, nu, nl);
        exit(1);
    }

    free(phi); free(gx); free(gy);
    L->K_lu = Kmat; L->K_piv = piv; L->n_loc = nloc;
    L->f_rhs = rhs; L->C_mat = Cmat; L->n_sig_k = n_sig_k; L->sig_map = sig_map;
}

static void local_block_free(LocalBlock *L) {
    free(L->K_lu); free(L->K_piv); free(L->f_rhs); free(L->C_mat); free(L->sig_map);
    memset(L, 0, sizeof(LocalBlock));
}

/* Retro-substitution locale (psi connu) : u_k, lambda_k puis coefficients RBF (c_k, d_k) */
static void backsub_local(DomainRect3FDD *D, int k, const LocalBlock *L) {
    RectSub *R = &D->sub[k];
    int nu = R->nu, nl = R->nl, nloc = L->n_loc;
    double *rhs_tot = (double *)xmalloc((size_t)nloc * sizeof(double));
    memcpy(rhs_tot, L->f_rhs, (size_t)nloc * sizeof(double));
    for (int l = 0; l < nl; l++) {
        double cpsi = 0.0;
        for (int c = 0; c < L->n_sig_k; c++) {
            cpsi += L->C_mat[l * L->n_sig_k + c] * D->psi[L->sig_map[c]];
        }
        rhs_tot[nu + l] -= cpsi;
    }
    lu_solve(L->K_lu, nloc, L->K_piv, rhs_tot);

    free(R->u_nodal); free(R->lambda); free(R->c_phs);
    R->u_nodal = (double *)xmalloc((size_t)nu * sizeof(double));
    R->lambda  = (double *)xmalloc((size_t)nl * sizeof(double));
    memcpy(R->u_nodal, rhs_tot, (size_t)nu * sizeof(double));
    memcpy(R->lambda,  rhs_tot + nu, (size_t)nl * sizeof(double));

    int dim = nu + 6;
    R->c_phs = (double *)xcalloc((size_t)dim, sizeof(double));
    memcpy(R->c_phs, R->u_nodal, (size_t)nu * sizeof(double));
    lu_solve(R->G_lu, dim, R->G_piv, R->c_phs);

    free(rhs_tot);
}

static void solve_rect_3fdd(DomainRect3FDD *D, double *t_out) {
    double t0 = now_sec();
    int M = D->ns, nS = D->n_sigma;
    double *S_glob = (double *)xcalloc((size_t)(nS * nS), sizeof(double));
    double *r_glob = (double *)xcalloc((size_t)nS, sizeof(double));
    LocalBlock *loc = (LocalBlock *)xcalloc((size_t)M, sizeof(LocalBlock));

    for (int k = 0; k < M; k++) {
        assemble_local(D, k, &loc[k]);
        int nu = D->sub[k].nu, nl = D->sub[k].nl, nloc = loc[k].n_loc, n_sig_k = loc[k].n_sig_k;
        double *Kmat = loc[k].K_lu, *Cmat = loc[k].C_mat, *rhs = loc[k].f_rhs;
        int *piv = loc[k].K_piv, *sig_map = loc[k].sig_map;
        double *sol0 = (double *)xmalloc((size_t)nloc * sizeof(double));
        memcpy(sol0, rhs, (size_t)nloc * sizeof(double));
        lu_solve(Kmat, nloc, piv, sol0);
        for (int c = 0; c < n_sig_k; c++) {
            double sv = 0.0;
            for (int l = 0; l < nl; l++) sv += Cmat[l * n_sig_k + c] * sol0[nu + l];
            r_glob[sig_map[c]] -= sv;
        }

        double *col_v = (double *)xmalloc((size_t)nloc * sizeof(double));
        for (int c = 0; c < n_sig_k; c++) {
            memset(col_v, 0, (size_t)nloc * sizeof(double));
            for (int l = 0; l < nl; l++) col_v[nu + l] = Cmat[l * n_sig_k + c];
            lu_solve(Kmat, nloc, piv, col_v);
            int gc = sig_map[c];
            for (int r = 0; r < n_sig_k; r++) {
                double sv = 0.0;
                for (int l = 0; l < nl; l++) sv += Cmat[l * n_sig_k + r] * col_v[nu + l];
                S_glob[sig_map[r] * nS + gc] -= sv;
            }
        }

        free(sol0); free(col_v);
    }
    /* 4. Résolution du complément de Schur global SPD sur Sigma */
    free(D->psi);
    D->psi = (double *)xcalloc((size_t)(nS > 0 ? nS : 1), sizeof(double));
    if (nS > 0) {
        int *pivS = (int *)xmalloc((size_t)nS * sizeof(int));
        memcpy(D->psi, r_glob, (size_t)nS * sizeof(double));
        if (lu_factor(S_glob, nS, pivS) != 0) {
            fprintf(stderr, "Erreur : complement de Schur singulier (nS=%d)\n", nS);
            exit(1);
        }
        lu_solve(S_glob, nS, pivS, D->psi);
        free(pivS);
    }

    /* 5. Back-substitution locale par rectangle */
    for (int k = 0; k < M; k++) { backsub_local(D, k, &loc[k]); local_block_free(&loc[k]); }

    free(loc); free(S_glob); free(r_glob);
    if (t_out) *t_out = now_sec() - t0;
}

/* ============================================================================
 * 7. Calcul des erreurs (L2, Linf, énergie, saut d'interface)
 * ============================================================================ */

typedef struct { double l2, linf, ener, jump; } ErrStats;

static ErrStats compute_errors(const DomainRect3FDD *D) {
    double e2 = 0.0, r2 = 0.0, en_e = 0.0, en_r = 0.0, linf = 0.0;
    for (int k = 0; k < D->ns; k++) {
        const RectSub *R = &D->sub[k];
        double K[2][2]; eval_K(R->cx, K);
        for (int qx_i = 0; qx_i < 5; qx_i++) {
            double qx = R->a + 0.5 * R->lx * (1.0 + GL5_X[qx_i]);
            double wx = 0.5 * R->lx * GL5_W[qx_i];
            for (int qy_i = 0; qy_i < 5; qy_i++) {
                double qy = R->c + 0.5 * R->ly * (1.0 + GL5_X[qy_i]);
                double w  = wx * 0.5 * R->ly * GL5_W[qy_i];
                double gx, gy, d[6];
                double uh = rect_eval_u_and_pde(R, qx, qy, &gx, &gy, NULL);
                uex_d(qx, qy, d);
                double diff = fabs(uh - d[0]);
                if (diff > linf) linf = diff;
                e2 += w * diff * diff;
                r2 += w * d[0] * d[0];
                double ex = gx - d[1], ey = gy - d[2];
                en_e += w * (K[0][0] * ex * ex + 2.0 * K[0][1] * ex * ey + K[1][1] * ey * ey);
                en_r += w * (K[0][0] * d[1] * d[1] + 2.0 * K[0][1] * d[1] * d[2] + K[1][1] * d[2] * d[2]);
            }
        }
    }
    double max_jump = 0.0;
    for (int f = 0; f < D->nif; f++) {
        const RectSub *R1 = &D->sub[D->ifs[f].k1], *R2 = &D->sub[D->ifs[f].k2];
        for (int q = 0; q < 5; q++) {
            double s = 0.5 * (1.0 + GL5_X[q]);
            double qx = D->ifs[f].x0 + s * (D->ifs[f].x1 - D->ifs[f].x0);
            double qy = D->ifs[f].y0 + s * (D->ifs[f].y1 - D->ifs[f].y0);
            double u1 = rect_eval_u_and_pde(R1, qx, qy, NULL, NULL, NULL);
            double u2 = rect_eval_u_and_pde(R2, qx, qy, NULL, NULL, NULL);
            double jp = fabs(u1 - u2);
            if (jp > max_jump) max_jump = jp;
        }
    }
    ErrStats E = {
        sqrt(e2 / (r2 > 1e-30 ? r2 : 1.0)),
        linf,
        sqrt(en_e / (en_r > 1e-30 ? en_r : 1.0)),
        max_jump
    };
    return E;
}

/* ============================================================================
 * 8. Double adaptativité :
 *    Niveau 1 : insertion de points ARS dans chaque rectangle (résidu EDP)
 *    Niveau 2 : bissection d'un rectangle en 2 fils si N_{u,k} > nsplit
 * ============================================================================ */

/* ============================================================================
 * 8bis. Adaptation du parametre alpha_W par residu local (option alpha=adapt)
 *   Pour chaque sous-domaine Omega_k, a trace psi fixee (probleme local bien pose),
 *   on teste alpha in {alpha_k/gamma, alpha_k, gamma*alpha_k} et on retient celui
 *   qui minimise le residu fort RMS  r_k(alpha) = D_k * RMS_m |f + div(K grad u_k)|
 *   sur la grille de test 8x8 (recherche lineaire geometrique, un pas par iteration).
 * ============================================================================ */

static double rect_local_residual_rms(const RectSub *R, int NT) {
    double diam = hypot(R->lx, R->ly), s2 = 0.0; int cnt = 0;
    for (int i = 1; i <= NT; i++)
        for (int j = 1; j <= NT; j++) {
            double qx = R->a + R->lx * (double)i / (double)(NT + 1);
            double qy = R->c + R->ly * (double)j / (double)(NT + 1);
            double Lu;
            rect_eval_u_and_pde(R, qx, qy, NULL, NULL, &Lu);
            double r = (Lu - f_src(qx, qy)) * diam;
            s2 += r * r; cnt++;
        }
    return sqrt(s2 / (double)(cnt > 0 ? cnt : 1));
}

/* Ecart RMS de trace  || u_k - psi ||_{L2(dOmega_k)} / sqrt|dOmega_k|  (psi = g_D sur dOmega) */
static double rect_trace_mismatch(const DomainRect3FDD *D, const RectSub *R) {
    double s2 = 0.0, len = 0.0;
    for (int e = 0; e < R->nseg; e++) {
        const BndSeg *sg = &R->seg[e];
        double elen = hypot(sg->x1 - sg->x0, sg->y1 - sg->y0);
        for (int q = 0; q < 5; q++) {
            double s = 0.5 * (1.0 + GL5_X[q]), w = 0.5 * elen * GL5_W[q];
            double qx = sg->x0 + s * (sg->x1 - sg->x0), qy = sg->y0 + s * (sg->y1 - sg->y0);
            double ps;
            if (sg->iface < 0) ps = u_ex(qx, qy);
            else {
                const IfaceSeg *I = &D->ifs[sg->iface];
                double s_if = (sg->iorient > 0) ? s : (1.0 - s);
                ps = 0.0;
                for (int m = 0; m < I->nc_e; m++) ps += D->psi[I->dof0 + m] * leg1d(m, s_if);
            }
            double d = rect_eval_u_and_pde(R, qx, qy, NULL, NULL, NULL) - ps;
            s2 += w * d * d;
        }
        len += elen;
    }
    return sqrt(s2 / (len > 0.0 ? len : 1.0));
}

/* Retourne (residu fort RMS, ecart de trace RMS) pour alpha donne, psi fixe */
static void rect_try_alpha(DomainRect3FDD *D, int k, double alpha, double *res, double *mis) {
    RectSub *R = &D->sub[k];
    R->alpha = alpha;
    LocalBlock L; memset(&L, 0, sizeof(L));
    assemble_local(D, k, &L);
    backsub_local(D, k, &L);
    local_block_free(&L);
    *res = rect_local_residual_rms(R, 8);
    *mis = rect_trace_mismatch(D, R);
}

/* Critere adimensionnel J(alpha) = r(alpha)/r(alpha_k) + m(alpha)/m(alpha_k) */
static void adapt_alpha_residual(DomainRect3FDD *D) {
    for (int k = 0; k < D->ns; k++) {
        RectSub *R = &D->sub[k];
        double a0 = (R->alpha > 0.0) ? R->alpha : ALPHA_W0;
        double cand[3] = { a0, a0 / ALPHA_GAMMA, a0 * ALPHA_GAMMA };
        double r[3], m[3], best_a = a0, best_J = 2.0;
        for (int c = 0; c < 3; c++) {
            double a = fmin(fmax(cand[c], ALPHA_MIN), ALPHA_MAX);
            cand[c] = a;
            rect_try_alpha(D, k, a, &r[c], &m[c]);
            double J = r[c] / fmax(r[0], 1e-300) + m[c] / fmax(m[0], 1e-300);
            if (c > 0 && J < best_J * 0.98) { best_J = J; best_a = a; }
        }
        double rr, mm;
        rect_try_alpha(D, k, best_a, &rr, &mm);   /* solution locale finale avec alpha_k retenu */
    }
}

static int adapt_rect_and_nodes(DomainRect3FDD *D, double thr_rel, int nsplit, int max_subs) {
    const int NT = 8;
    double global_rmax = 0.0;
    if (RBF_TYPE == 1 && ALPHA_ADAPT) adapt_alpha_residual(D);

    for (int k = 0; k < D->ns; k++) {
        RectSub *R = &D->sub[k];
        double diam = hypot(R->lx, R->ly);
        double rmax = 0.0;
        for (int i = 1; i <= NT; i++) {
            for (int j = 1; j <= NT; j++) {
                double qx = R->a + R->lx * (double)i / (double)(NT + 1);
                double qy = R->c + R->ly * (double)j / (double)(NT + 1);
                double Lu;
                rect_eval_u_and_pde(R, qx, qy, NULL, NULL, &Lu);
                double res = fabs(Lu - f_src(qx, qy)) * diam;
                if (res > rmax) rmax = res;
            }
        }
        R->res_max = rmax;
        if (rmax > global_rmax) global_rmax = rmax;
    }

    double thresh = thr_rel * global_rmax;
    int pts_added = 0;

    /* Niveau 1 : raffinement ARS des nœuds dans les rectangles à fort résidu */
    for (int k = 0; k < D->ns; k++) {
        RectSub *R = &D->sub[k];
        if (R->res_max < thresh) continue;
        double diam = hypot(R->lx, R->ly);
        double dmin = 0.07 * fmin(R->lx, R->ly);
        int added_k = 0;
        for (int i = 1; i <= NT && added_k < 12; i++) {
            for (int j = 1; j <= NT && added_k < 12; j++) {
                double qx = R->a + R->lx * (double)i / (double)(NT + 1);
                double qy = R->c + R->ly * (double)j / (double)(NT + 1);
                double Lu;
                rect_eval_u_and_pde(R, qx, qy, NULL, NULL, &Lu);
                double res = fabs(Lu - f_src(qx, qy)) * diam;
                if (res >= thresh) {
                    int old_nu = R->nu;
                    rect_add_pt(R, qx, qy, dmin);
                    if (R->nu > old_nu) { added_k++; pts_added++; }
                }
            }
        }
    }

    /* Niveau 2 : division d'un rectangle en 2 fils selon son plus grand côté si nu > nsplit */
    int splits = 0;
    int old_ns = D->ns;
    int *was_split = (int *)xcalloc((size_t)(max_subs + 16), sizeof(int));

    for (int k = 0; k < old_ns && D->ns < max_subs; k++) {
        if (D->sub[k].nu <= nsplit) continue;
        if (D->ns >= D->cap_s) {
            D->cap_s *= 2;
            D->sub = (RectSub *)xrealloc(D->sub, (size_t)D->cap_s * sizeof(RectSub));
        }
        RectSub *R1 = &D->sub[k];
        RectSub *R2 = &D->sub[D->ns];
        memset(R2, 0, sizeof(RectSub));
        int new_id = D->ns++;

        R2->a = R1->a; R2->b = R1->b; R2->c = R1->c; R2->d = R1->d;
        R2->alpha = R1->alpha;   /* le fils herite du alpha_k adapte du pere */
        R1->lev++; R2->lev = R1->lev;

        if (R1->lx >= R1->ly) {
            double xm = 0.5 * (R1->a + R1->b);
            R1->b = xm;
            R2->a = xm;
        } else {
            double ym = 0.5 * (R1->c + R1->d);
            R1->d = ym;
            R2->c = ym;
        }
        was_split[k] = 1;
        was_split[new_id] = 1;
        splits++;
    }

    /* Equilibrage 2:1 : un rectangle est divise s'il touche un voisin plus de 2.5 fois plus petit
     * (limite le nombre de noeuds pendants par cote et donc de modes lambda face a peu de centres) */
    for (int changed = 1; changed && D->ns < max_subs; ) {
        changed = 0;
        int nsc = D->ns;
        for (int k = 0; k < nsc && D->ns < max_subs; k++) {
            RectSub *A = &D->sub[k];
            double sa = fmax(A->b - A->a, A->d - A->c), need = 0;
            for (int j = 0; j < nsc && !need; j++) {
                if (j == k) continue;
                const RectSub *B = &D->sub[j];
                double ox = fmin(A->b, B->b) - fmax(A->a, B->a), oy = fmin(A->d, B->d) - fmax(A->c, B->c);
                int adj = (fabs(ox) < 1e-12 && oy > 1e-12) || (fabs(oy) < 1e-12 && ox > 1e-12);
                if (adj && sa > 2.5 * fmax(B->b - B->a, B->d - B->c)) need = 1;
            }
            if (!need) continue;
            if (D->ns >= D->cap_s) { D->cap_s *= 2; D->sub = (RectSub *)xrealloc(D->sub, (size_t)D->cap_s * sizeof(RectSub)); }
            RectSub *R1 = &D->sub[k], *R2 = &D->sub[D->ns];
            memset(R2, 0, sizeof(RectSub));
            int new_id = D->ns++;
            R2->a = R1->a; R2->b = R1->b; R2->c = R1->c; R2->d = R1->d; R2->alpha = R1->alpha;
            R1->lev++; R2->lev = R1->lev;
            if (R1->b - R1->a >= R1->d - R1->c) { double xm = 0.5 * (R1->a + R1->b); R1->b = xm; R2->a = xm; }
            else { double ym = 0.5 * (R1->c + R1->d); R1->d = ym; R2->c = ym; }
            was_split[k] = 1; was_split[new_id] = 1; splits++; changed = 1;
        }
    }

    if (splits > 0) {
        build_rect_interfaces(D);
        for (int k = 0; k < D->ns; k++) {
            if (was_split[k]) {
                rect_seed_points(&D->sub[k], N_INIT - 1, N_INIT - 1);
            } else {
                /* S'assurer que les nouveaux sommets d'interface non-conformes sont présents */
                RectSub *R = &D->sub[k];
                double dmin = 0.008 * fmin(R->lx, R->ly);
                for (int e = 0; e < R->nseg; e++) {
                    const BndSeg *sg = &R->seg[e];
                    int m = sg->nb_e + 1;
                    for (int i = 0; i <= m; i++) {
                        double t = (double)i / (double)m;
                        rect_add_pt(R, sg->x0 + t * (sg->x1 - sg->x0),
                                       sg->y0 + t * (sg->y1 - sg->y0), dmin);
                    }
                }
            }
        }
    }
    free(was_split);
    return pts_added + splits;
}

/* ============================================================================
 * 9. Exports Gnuplot (rectangles, nœuds, interfaces Sigma, historique)
 * ============================================================================ */

static void export_snapshot(const char *prefix, const DomainRect3FDD *D) {
    char fn[256];
    snprintf(fn, sizeof(fn), "%s_rects.dat", prefix);
    FILE *fr = fopen(fn, "w");
    if (fr) {
        for (int k = 0; k < D->ns; k++) {
            const RectSub *R = &D->sub[k];
            fprintf(fr, "%.6f %.6f\n%.6f %.6f\n%.6f %.6f\n%.6f %.6f\n%.6f %.6f\n\n",
                    R->a, R->c, R->b, R->c, R->b, R->d, R->a, R->d, R->a, R->c);
        }
        fclose(fr);
    }
    snprintf(fn, sizeof(fn), "%s_nodes.dat", prefix);
    FILE *fn_p = fopen(fn, "w");
    if (fn_p) {
        for (int k = 0; k < D->ns; k++) {
            const RectSub *R = &D->sub[k];
            for (int i = 0; i < R->nu; i++) {
                fprintf(fn_p, "%.6f %.6f %.6e %d\n", R->ux[i], R->uy[i], R->u_nodal[i], k);
            }
        }
        fclose(fn_p);
    }
}

/* ============================================================================
 * 9b. Affichage temps réel via pipe Gnuplot (multiplot 1x2 : domaine + surface 3D)
 * ============================================================================ */

static void export_surface_lshape(const char *fname, const DomainRect3FDD *D, int ng) {
    FILE *fs = fopen(fname, "w");
    if (!fs) return;
    for (int i = 0; i <= ng; i++) {
        double x = 2.0 * (double)i / (double)ng;
        for (int j = 0; j <= ng; j++) {
            double y = 2.0 * (double)j / (double)ng;
            if (x > 1.0 + 1e-9 && y > 1.0 + 1e-9) {
                fprintf(fs, "%.6f %.6f NaN\n", x, y);
                continue;
            }
            /* Trouver le rectangle contenant (x, y) */
            int hit = -1;
            for (int k = 0; k < D->ns; k++) {
                const RectSub *R = &D->sub[k];
                if (x >= R->a - 1e-9 && x <= R->b + 1e-9 &&
                    y >= R->c - 1e-9 && y <= R->d + 1e-9) { hit = k; break; }
            }
            if (hit >= 0) {
                double u = rect_eval_u_and_pde(&D->sub[hit], x, y, NULL, NULL, NULL);
                fprintf(fs, "%.6f %.6f %.6e\n", x, y, u);
            } else {
                fprintf(fs, "%.6f %.6f NaN\n", x, y);
            }
        }
        fprintf(fs, "\n");
    }
    fclose(fs);
}

static void plot_realtime(
    FILE *gp, int ex_id, int iter, const DomainRect3FDD *D,
    ErrStats E, int delay_ms, int is_gif
) {
    if (!gp) return;
    export_snapshot("live_frame", D);
    export_surface_lshape("live_surf.dat", D, 44);

    FILE *fsig = fopen("live_sigma.dat", "w");
    if (fsig) {
        for (int f = 0; f < D->nif; f++) {
            fprintf(fsig, "%.6f %.6f\n%.6f %.6f\n\n",
                    D->ifs[f].x0, D->ifs[f].y0, D->ifs[f].x1, D->ifs[f].y1);
        }
        fclose(fsig);
    }

    int nu_tot = 0;
    for (int k = 0; k < D->ns; k++) nu_tot += D->sub[k].nu;

    if (!is_gif) fprintf(gp, "clear\n");
    fprintf(gp, "set multiplot layout 1,2 title '3FDD-RBF Rectangles (Ex. %d) | Iter %d | Domaines M=%d | N_u=%d, N_Sig=%d | e_L2=%.2e'\n",
            ex_id, iter, D->ns, nu_tot, D->n_sigma, E.l2);

    /* Panneau gauche : domaine non convexe en L, rectangles, interfaces Sigma et nœuds ARS */
    fprintf(gp, "set title 'Decomposition (%d rectangles) + noeuds ARS'\n", D->ns);
    fprintf(gp, "set xrange [-0.05:2.05]; set yrange [-0.05:2.05]; set size ratio -1\n");
    fprintf(gp, "set xlabel 'x'; set ylabel 'y'; unset key\n");
    fprintf(gp, "plot 'live_frame_rects.dat' with lines lw 1.2 lc rgb '#222222', "
                "'live_sigma.dat' with lines lw 2.0 lc rgb '#0055CC', "
                "'live_frame_nodes.dat' using 1:2 with points pt 7 ps 0.55 lc rgb '#CC3311'\n");

    /* Panneau droit : solution 3D u_h(x, y) sur le domaine en L */
    fprintf(gp, "set title 'Solution approchee u_h(x,y) sur domaine en L'\n");
    fprintf(gp, "set xrange [0:2]; set yrange [0:2]; set size nosquare\n");
    fprintf(gp, "set pm3d; set hidden3d; set view 58, 35\n");
    fprintf(gp, "set xlabel 'x'; set ylabel 'y'; set zlabel 'u_h'\n");
    fprintf(gp, "splot 'live_surf.dat' using 1:2:3 with pm3d\n");

    fprintf(gp, "unset multiplot\n");
    fflush(gp);
    if (!is_gif && delay_ms > 0) SLEEP_MS(delay_ms);
}

/* ============================================================================
 * 10. Initialisation du domaine non convexe en L (EXACTEMENT 2 rectangles)
 * ============================================================================ */

static void init_lshape_2rects(DomainRect3FDD *D, int nx0, int ny0, int nb, int nc) {
    memset(D, 0, sizeof(DomainRect3FDD));
    D->nb = nb; D->nc = nc;
    D->cap_s = 64;
    D->sub = (RectSub *)xcalloc((size_t)D->cap_s, sizeof(RectSub));
    D->ns = 2;

    /* Rectangle 1 : [0, 1] x [0, 2] ; Rectangle 2 : [1, 2] x [0, 1] */
    D->sub[0].a = 0.0; D->sub[0].b = 1.0; D->sub[0].c = 0.0; D->sub[0].d = 2.0; D->sub[0].lev = 0;
    D->sub[1].a = 1.0; D->sub[1].b = 2.0; D->sub[1].c = 0.0; D->sub[1].d = 1.0; D->sub[1].lev = 0;

    build_rect_interfaces(D);
    (void)nx0; (void)ny0;
    /* decomposition minimale du L (2 rectangles) et n x n points reguliers par sous-domaine */
    rect_seed_points(&D->sub[0], N_INIT - 1, N_INIT - 1);
    rect_seed_points(&D->sub[1], N_INIT - 1, N_INIT - 1);
}

static void domain_free(DomainRect3FDD *D) {
    for (int k = 0; k < D->ns; k++) rect_free_arrays(&D->sub[k]);
    free(D->sub); free(D->ifs); free(D->psi);
}

static int total_nodes(const DomainRect3FDD *D) {
    int s = 0;
    for (int k = 0; k < D->ns; k++) s += D->sub[k].nu;
    return s;
}

static int MMAX = 0;        /* mmax= : nombre max de sous-domaines (0 : valeur par exemple) */
static int PRINT_PTS = 0;   /* pts=1 : tableau u_h - u_exact en chaque centre */
static double point_errors(const DomainRect3FDD *D, int ex_id, int final) {
    double emax = 0.0; FILE *f = NULL;
    if (final) { char fn[64]; snprintf(fn, sizeof(fn), "ex%d_erreur_points.dat", ex_id); f = fopen(fn, "w"); }
    if (f) fprintf(f, "# k i x y u_h u_exact u_h-u_exact   (k = sous-domaine)\n");
    if (final && PRINT_PTS) printf("   %4s %4s %10s %10s %14s %14s %12s\n", "k", "i", "x", "y", "u_h", "u_exact", "u_h-u_ex");
    for (int k = 0; k < D->ns; k++) {
        const RectSub *R = &D->sub[k];
        for (int i = 0; i < R->nu; i++) {
            double ue = u_ex(R->ux[i], R->uy[i]), d = R->u_nodal[i] - ue;
            if (fabs(d) > emax) emax = fabs(d);
            if (f) fprintf(f, "%d %d %.8f %.8f %.10e %.10e %.4e\n", k, i, R->ux[i], R->uy[i], R->u_nodal[i], ue, d);
            if (final && PRINT_PTS) printf("   %4d %4d %10.5f %10.5f %14.8f %14.8f %12.3e\n", k, i, R->ux[i], R->uy[i], R->u_nodal[i], ue, d);
        }
    }
    if (f) fclose(f);
    return emax;
}

static void run_example(
    int ex_id, const char *title, int max_iter, double thr_rel,
    int nsplit, int max_subs, FILE *gp, int delay_ms, int is_gif
) {
    EX_ID = ex_id;
    printf("\n==========================================================================================\n");
    printf(" EXEMPLE %d : %s\n", ex_id, title);
    printf("   Base %s + P2%s, %d x %d points par sous-domaine\n", RBF_NAME[RBF_TYPE],
           (RBF_TYPE >= 4) ? " (s = kappa h_k)" : "", N_INIT, N_INIT);
    printf("   Domaine non convexe en L = [0,1]x[0,2] U [1,2]x[0,1] (depart M = 2 sous-domaines)\n");
    printf("==========================================================================================\n");
    printf("%4s %7s %9s %6s %6s %10s %10s %10s %10s %6s\n",
           "Iter", "Nb_Dom", "nu(m..M)", "N_u", "N_Sig", "e_L2", "e_max", "e_Linf_q", "saut_Sig", "t(s)");
    printf("------------------------------------------------------------------------------------------\n");

    DomainRect3FDD D;
    /* nb modes de Legendre pour lambda_k, nc = nb - 1 pour psi (these : nb = 4, nc = 3).
     * Si la grille n x n porte moins de nb+1 intervalles par cote, des points de bord sont ajoutes (inf-sup). */
    if (N_INIT < 3) N_INIT = 3;
    int nb_l = NB_LEG, nc_l = nb_l - 1;
    init_lshape_2rects(&D, 5, 5, nb_l, nc_l);
    printf("   Multiplicateurs : nb = %d modes (lambda_k), nc = %d modes (psi sur Sigma)\n", nb_l, nc_l);
    /* seuil de division relatif a la densite initiale : un sous-domaine neuf (n x n points)
     * ne doit pas etre redivise immediatement */
    if (nsplit < N_INIT * N_INIT + 8) nsplit = N_INIT * N_INIT + 8;
    if (MMAX > 0) max_subs = MMAX;

    char fhist[64], pref0[64], prefF[64];
    snprintf(fhist, sizeof(fhist), "ex%d_history.dat", ex_id);
    snprintf(pref0, sizeof(pref0), "ex%d_init", ex_id);
    snprintf(prefF, sizeof(prefF), "ex%d_final", ex_id);
    FILE *fh = fopen(fhist, "w");
    if (fh) fprintf(fh, "# iter Nb_Dom N_u N_Sig L2 Linf ener jump t alpha_min alpha_max emax_noeuds\n");

    for (int it = 0; it < max_iter; it++) {
        double ts = 0.0;
        solve_rect_3fdd(&D, &ts);
        ErrStats E = compute_errors(&D);
        int nu_tot = total_nodes(&D);
        int nu_min = D.sub[0].nu, nu_max = D.sub[0].nu;
        for (int k = 1; k < D.ns; k++) {
            if (D.sub[k].nu < nu_min) nu_min = D.sub[k].nu;
            if (D.sub[k].nu > nu_max) nu_max = D.sub[k].nu;
        }
        char range_s[16];
        snprintf(range_s, sizeof(range_s), "%2d..%-2d", nu_min, nu_max);

        double emax = point_errors(&D, ex_id, 0);
        printf("%4d %7d %9s %6d %6d %10.3e %10.3e %10.3e %10.3e %6.3f\n",
               it, D.ns, range_s, nu_tot, D.n_sigma, E.l2, emax, E.linf, E.jump, ts);
        double amin = 1e300, amax = 0.0;
        for (int k = 0; k < D.ns; k++) {
            double a = (D.sub[k].alpha > 0.0) ? D.sub[k].alpha : ALPHA_W0;
            if (a < amin) amin = a;
            if (a > amax) amax = a;
        }
        if (RBF_TYPE == 1) printf("       alpha_k in [%.2f, %.2f]\n", amin, amax);
        if (fh) fprintf(fh, "%d %d %d %d %.6e %.6e %.6e %.6e %.4f %.3f %.3f %.6e\n",
                        it, D.ns, nu_tot, D.n_sigma, E.l2, E.linf, E.ener, E.jump, ts, amin, amax, emax);
        plot_realtime(gp, ex_id, it, &D, E, delay_ms, is_gif);
        if (it == 0) export_snapshot(pref0, &D);
        int fin = (it == max_iter - 1);
        if (fin) {
            export_snapshot(prefF, &D);
        } else {
            int chg = adapt_rect_and_nodes(&D, thr_rel, nsplit, max_subs);
            if (chg == 0) {
                fin = 1;
                export_snapshot(prefF, &D);
                break;
            }
        }
    }
    if (fh) fclose(fh);
    {
        if (PRINT_PTS) printf("   --- 3FDD-Rectangles, Ex. %d : erreur ponctuelle aux centres (fichier ex%d_erreur_points.dat) ---\n", ex_id, ex_id);
        double emax = point_errors(&D, ex_id, 1);
        ErrStats E = compute_errors(&D);
        printf("   Bilan : %d sous-domaines, %d centres au total, e_L2 = %.3e, e_max(centres) = %.3e\n",
               D.ns, total_nodes(&D), E.l2, emax);
    }
    domain_free(&D);
}

int main(int argc, char **argv) {
    int gui = 1, delay_ms = 650, only_ex = 0;
    const char *gif_out = NULL;

    for (int i = 1; i < argc; i++) {
        if (!strncmp(argv[i], "gui=", 4))        gui = atoi(argv[i] + 4);
        else if (!strncmp(argv[i], "delay=", 6)) delay_ms = atoi(argv[i] + 6);
        else if (!strncmp(argv[i], "ex=", 3))    only_ex = atoi(argv[i] + 3);
        else if (!strcmp(argv[i], "-o") && i + 1 < argc) gif_out = argv[++i];
        else if (!strncmp(argv[i], "out=", 4))   gif_out = argv[i] + 4;
        else if (!strcmp(argv[i], "rbf=wendland") || !strcmp(argv[i], "rbf=w4")) RBF_TYPE = 1;
        else if (!strcmp(argv[i], "rbf=phs") || !strcmp(argv[i], "rbf=phs5")) RBF_TYPE = 0;
        else if (!strcmp(argv[i], "rbf=phs3"))   RBF_TYPE = 2;
        else if (!strcmp(argv[i], "rbf=tps"))    RBF_TYPE = 3;
        else if (!strcmp(argv[i], "rbf=gauss"))  RBF_TYPE = 4;
        else if (!strcmp(argv[i], "rbf=mq"))     RBF_TYPE = 5;
        else if (!strcmp(argv[i], "rbf=imq"))    RBF_TYPE = 6;
        else if (!strncmp(argv[i], "kappa=", 6)) KAPPA = atof(argv[i] + 6);
        else if (!strncmp(argv[i], "n=", 2))     N_INIT = atoi(argv[i] + 2);
        else if (!strncmp(argv[i], "nb=", 3))    NB_LEG = atoi(argv[i] + 3);
        else if (!strncmp(argv[i], "pts=", 4))   PRINT_PTS = atoi(argv[i] + 4);
        else if (!strncmp(argv[i], "mmax=", 5))  MMAX = atoi(argv[i] + 5);
        else if (!strcmp(argv[i], "alpha=adapt")) ALPHA_ADAPT = 1;
        else if (!strncmp(argv[i], "alpha=", 6))  ALPHA_W0 = atof(argv[i] + 6);
    }

    int is_gif = (gif_out != NULL);
    FILE *gp = NULL;
    if (is_gif) {
        gp = POPEN("gnuplot", "w");
        if (gp) {
            int gif_delay = delay_ms / 10; /* en centièmes de seconde */
            if (gif_delay < 10) gif_delay = 10;
            fprintf(gp, "set terminal gif animate optimize delay %d size 1100,520\n", gif_delay);
            fprintf(gp, "set output '%s'\n", gif_out);
            fflush(gp);
        } else {
            fprintf(stderr, "Erreur : impossible de lancer gnuplot pour generer %s\n", gif_out);
        }
    } else if (gui) {
        gp = POPEN("gnuplot -persistent", "w");
        if (!gp) fprintf(stderr, "Avertissement : impossible de lancer gnuplot (mode console seul).\n");
    }

    if (only_ex == 0 || only_ex == 1)
        run_example(1, "Lisse, tenseur plein constant K=[[2,-1],[-1,2]]",                6, 0.25, 44, 128, gp, delay_ms, is_gif);
    if (only_ex == 0 || only_ex == 2)
        run_example(2, "Tenseur plein discontinu x=1 (I vs [[8,-7],[-7,20]])",          6, 0.25, 44, 128, gp, delay_ms, is_gif);
    if (only_ex == 0 || only_ex == 3)
        run_example(3, "Forte anisotropie 50:1 (30 deg) + pic raide au coin rentrant", 10, 0.15, 42, 128, gp, delay_ms, is_gif);

    if (gp) {
        if (is_gif) { fprintf(gp, "set output\n"); fflush(gp); }
        PCLOSE(gp);
        if (is_gif) printf("\nAnimation GIF enregistree dans : %s\n", gif_out);
    }
    return 0;
}
