/* sge.c -- structured sparse Gaussian elimination mod p (p < 2^31) with Markowitz-style pivoting.
 *
 * Rows are sorted sparse vectors; column occupancy lists are kept lazily (validated on use) while the
 * exact per-column counts of active rows are maintained and bucketed, so the cheapest pivot column is
 * found in amortized O(1).  The caller stops the elimination when the active submatrix is small/dense
 * enough (status 1/2/3) and finishes it densely; sge_backsub then recovers the pivot unknowns.
 *
 * Column rhs_col (if >= 0) holds the right-hand side and is never pivoted; a row whose only entry is in
 * rhs_col proves inconsistency (status -1).
 *
 * Loaded from python via ctypes (hodge/sge.py).  Compile: gcc -O3 -shared -fPIC -o sge.so sge.c
 */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef uint32_t u32; typedef uint64_t u64; typedef int32_t i32; typedef int64_t i64;

typedef struct { i32 len, cap; i32 *c; u32 *v; } Row;
typedef struct { i32 len, cap; i32 *r; } CList;

typedef struct {
    i32 nrows, ncols; u32 p;
    Row *rows;
    unsigned char *rstate;      /* 0 removed, 1 active, 2 pivot */
    i32 *rpiv;                  /* pivot column of pivot rows */
    CList *cols;
    i32 *ccnt;                  /* exact number of active rows containing the column */
    unsigned char *cstate;      /* 0 active, 1 pivot, 2 retired (empty) */
    i32 *bnext, *bprev, *bhead; /* buckets by count */
    i32 curmin;
    i32 *piv_order; i32 npiv;
    i32 rhs_col;
    int inconsistent;
    i64 active_nnz; i32 active_rows, active_cols;
    i32 *tmpc; u32 *tmpv; i32 tmpcap;
    i64 fill_events;
} SGE;

static u64 inv_mod(u64 a, u64 p) {
    i64 t = 0, newt = 1; i64 r = (i64)p, newr = (i64)(a % p);
    while (newr) { i64 q = r / newr; i64 tmp = t - q * newt; t = newt; newt = tmp; tmp = r - q * newr; r = newr; newr = tmp; }
    if (t < 0) t += (i64)p;
    return (u64)t;
}

static void bucket_remove(SGE *h, i32 c) {
    i32 n = h->bnext[c], pv = h->bprev[c];
    if (pv >= 0) h->bnext[pv] = n; else h->bhead[h->ccnt[c]] = n;
    if (n >= 0) h->bprev[n] = pv;
    h->bnext[c] = h->bprev[c] = -1;
}

static void bucket_insert(SGE *h, i32 c) {
    i32 k = h->ccnt[c];
    h->bprev[c] = -1; h->bnext[c] = h->bhead[k];
    if (h->bhead[k] >= 0) h->bprev[h->bhead[k]] = c;
    h->bhead[k] = c;
    if (k < h->curmin) h->curmin = k;
}

/* count change for an active column */
static void col_count_add(SGE *h, i32 c, i32 delta) {
    if (h->cstate[c] != 0) return;
    bucket_remove(h, c);
    h->ccnt[c] += delta;
    bucket_insert(h, c);
}

static void clist_push(CList *L, i32 r) {
    if (L->len == L->cap) { L->cap = L->cap ? 2 * L->cap : 4; L->r = (i32 *)realloc(L->r, sizeof(i32) * L->cap); }
    L->r[L->len++] = r;
}

static int row_find(const Row *R, i32 c) {
    i32 lo = 0, hi = R->len - 1;
    while (lo <= hi) { i32 m = (lo + hi) >> 1; if (R->c[m] == c) return m; if (R->c[m] < c) lo = m + 1; else hi = m - 1; }
    return -1;
}

SGE *sge_create(i32 nrows, i32 ncols, const i64 *indptr, const i32 *indices, const u32 *vals, u32 p, i32 rhs_col) {
    SGE *h = (SGE *)calloc(1, sizeof(SGE));
    h->nrows = nrows; h->ncols = ncols; h->p = p; h->rhs_col = rhs_col; h->curmin = 0;
    h->rows = (Row *)calloc(nrows, sizeof(Row));
    h->rstate = (unsigned char *)calloc(nrows, 1);
    h->rpiv = (i32 *)malloc(sizeof(i32) * nrows);
    h->cols = (CList *)calloc(ncols, sizeof(CList));
    h->ccnt = (i32 *)calloc(ncols, sizeof(i32));
    h->cstate = (unsigned char *)calloc(ncols, 1);
    h->bnext = (i32 *)malloc(sizeof(i32) * ncols); h->bprev = (i32 *)malloc(sizeof(i32) * ncols);
    h->bhead = (i32 *)malloc(sizeof(i32) * (nrows + 2));
    for (i32 k = 0; k <= nrows + 1; k++) h->bhead[k] = -1;
    h->piv_order = (i32 *)malloc(sizeof(i32) * (ncols + 1)); h->npiv = 0;
    h->tmpcap = 1024; h->tmpc = (i32 *)malloc(sizeof(i32) * h->tmpcap); h->tmpv = (u32 *)malloc(sizeof(u32) * h->tmpcap);
    for (i32 r = 0; r < nrows; r++) {
        i64 a = indptr[r], b = indptr[r + 1]; i32 len = (i32)(b - a);
        Row *R = &h->rows[r];
        R->len = len; R->cap = len > 0 ? len : 1;
        R->c = (i32 *)malloc(sizeof(i32) * R->cap); R->v = (u32 *)malloc(sizeof(u32) * R->cap);
        /* entries are assumed sorted by column with no duplicates and nonzero values (python side guarantees) */
        for (i32 k = 0; k < len; k++) { R->c[k] = indices[a + k]; R->v[k] = vals[a + k] % p; }
        if (len == 0) { h->rstate[r] = 0; continue; }
        h->rstate[r] = 1; h->active_rows++; h->active_nnz += len;
        for (i32 k = 0; k < len; k++) { h->ccnt[R->c[k]]++; clist_push(&h->cols[R->c[k]], r); }
    }
    for (i32 c = 0; c < ncols; c++) {
        h->bnext[c] = h->bprev[c] = -1;
        if (c == rhs_col) { h->cstate[c] = 3; continue; }       /* 3 = rhs: never pivoted, never counted */
        if (h->ccnt[c] == 0) { h->cstate[c] = 2; continue; }
        h->cstate[c] = 0; h->active_cols++; bucket_insert(h, c);
    }
    h->curmin = 1;
    /* a row whose only entry is the rhs => inconsistent */
    for (i32 r = 0; r < nrows; r++)
        if (h->rstate[r] == 1 && h->rows[r].len == 1 && h->rows[r].c[0] == rhs_col) h->inconsistent = 1;
    return h;
}

void sge_free(SGE *h) {
    if (!h) return;
    for (i32 r = 0; r < h->nrows; r++) { free(h->rows[r].c); free(h->rows[r].v); }
    for (i32 c = 0; c < h->ncols; c++) free(h->cols[c].r);
    free(h->rows); free(h->rstate); free(h->rpiv); free(h->cols); free(h->ccnt); free(h->cstate);
    free(h->bnext); free(h->bprev); free(h->bhead); free(h->piv_order); free(h->tmpc); free(h->tmpv); free(h);
}

/* remove row r from the active set (count updates); state becomes `newstate` */
static void row_deactivate(SGE *h, i32 r, unsigned char newstate) {
    Row *R = &h->rows[r];
    for (i32 k = 0; k < R->len; k++) col_count_add(h, R->c[k], -1);
    h->rstate[r] = newstate; h->active_rows--; h->active_nnz -= R->len;
}

/* row_k <- row_k - f * row_p  (mod p); maintains counts / lists / nnz.  Column pivcol is being retired. */
static void row_axpy(SGE *h, i32 k, const Row *P, u64 f, i32 pivcol) {
    Row *K = &h->rows[k]; u64 p = h->p;
    i32 need = K->len + P->len;
    if (need > h->tmpcap) { h->tmpcap = 2 * need; h->tmpc = (i32 *)realloc(h->tmpc, sizeof(i32) * h->tmpcap); h->tmpv = (u32 *)realloc(h->tmpv, sizeof(u32) * h->tmpcap); }
    i32 i = 0, j = 0, n = 0; i32 oldlen = K->len;
    while (i < K->len || j < P->len) {
        i32 ck = i < K->len ? K->c[i] : 0x7fffffff, cp = j < P->len ? P->c[j] : 0x7fffffff;
        if (ck < cp) { h->tmpc[n] = ck; h->tmpv[n] = K->v[i]; n++; i++; }
        else if (cp < ck) {
            u64 val = (p - (f * (u64)P->v[j]) % p) % p;
            if (val) { h->tmpc[n] = cp; h->tmpv[n] = (u32)val; n++; if (cp != pivcol) { col_count_add(h, cp, +1); clist_push(&h->cols[cp], k); h->fill_events++; } }
            j++;
        } else {
            u64 val = ((u64)K->v[i] + p - (f * (u64)P->v[j]) % p) % p;
            if (val) { h->tmpc[n] = ck; h->tmpv[n] = (u32)val; n++; }
            else if (ck != pivcol) col_count_add(h, ck, -1);
            i++; j++;
        }
    }
    if (n > K->cap) { K->cap = n + n / 4 + 4; K->c = (i32 *)realloc(K->c, sizeof(i32) * K->cap); K->v = (u32 *)realloc(K->v, sizeof(u32) * K->cap); }
    memcpy(K->c, h->tmpc, sizeof(i32) * n); memcpy(K->v, h->tmpv, sizeof(u32) * n); K->len = n;
    h->active_nnz += (i64)n - oldlen;
    if (n == 0) { h->rstate[k] = 0; h->active_rows--; }
    else if (n == 1 && K->c[0] == h->rhs_col) { h->inconsistent = 1; }
}

/* choose pivot: scan the lowest buckets, up to `cands` candidate columns, minimize (len-1)*(cnt-1).
   returns 1 and sets *pr,*pc,*cost; 0 if no active column remains */
static int choose_pivot(SGE *h, int cands, i32 *pr, i32 *pc, i64 *cost) {
    while (h->curmin <= h->nrows && h->bhead[h->curmin] < 0) h->curmin++;
    if (h->curmin > h->nrows) return 0;
    i64 best = -1; i32 br = -1, bc = -1; int seen = 0;
    for (i32 b = h->curmin; b <= h->nrows && seen < cands; b++) {
        for (i32 c = h->bhead[b]; c >= 0 && seen < cands; c = h->bnext[c]) {
            seen++;
            /* find the shortest active row containing c (validate lazily, compact the list) */
            CList *L = &h->cols[c]; i32 w = 0; i32 bestr = -1; i32 bestlen = 0x7fffffff;
            for (i32 t = 0; t < L->len; t++) {
                i32 r = L->r[t];
                if (h->rstate[r] != 1) continue;
                if (row_find(&h->rows[r], c) < 0) continue;
                L->r[w++] = r;
                if (h->rows[r].len < bestlen) { bestlen = h->rows[r].len; bestr = r; }
            }
            L->len = w;
            if (bestr < 0) continue;      /* should not happen (counts exact) */
            i64 cst = (i64)(bestlen - 1) * (i64)(b - 1);
            if (best < 0 || cst < best) { best = cst; br = bestr; bc = c; }
            if (best == 0) break;
        }
        if (best == 0) break;
    }
    if (br < 0) return 0;
    *pr = br; *pc = bc; *cost = best; return 1;
}

/* Eliminate up to max_pivots pivots.  Stops early with:
 *  -1 inconsistent (only when stop_on_inconsistent; otherwise the sticky flag is set and elimination goes on so the
 *     rank of the coefficient part is still available), 0 finished (no active columns / rows),
 *   1 active_cols <= stop_cols, 2 best pivot cost > max_cost, 3 active_nnz > max_nnz, 4 max_pivots reached. */
int sge_eliminate(SGE *h, i32 max_pivots, i32 stop_cols, i64 max_cost, i64 max_nnz, int cands, int stop_on_inconsistent) {
    for (i32 it = 0; it < max_pivots; it++) {
        if (h->inconsistent && stop_on_inconsistent) return -1;
        if (h->active_cols <= stop_cols || h->active_rows == 0) return h->active_cols == 0 || h->active_rows == 0 ? 0 : 1;
        if (h->active_nnz > max_nnz) return 3;
        i32 r, c; i64 cost;
        if (!choose_pivot(h, cands, &r, &c, &cost)) return 0;
        if (cost > max_cost) return 2;
        /* pivot row r, column c */
        row_deactivate(h, r, 2);
        h->rpiv[r] = c; h->piv_order[h->npiv++] = r;
        Row *P = &h->rows[r];
        i32 kp = row_find(P, c);
        u64 inv = inv_mod(P->v[kp], h->p);
        CList *L = &h->cols[c];
        /* copy the list: row_axpy may append to other lists, not to this one (c is being retired) */
        for (i32 t = 0; t < L->len; t++) {
            i32 k = L->r[t];
            if (h->rstate[k] != 1) continue;
            i32 kk = row_find(&h->rows[k], c);
            if (kk < 0) continue;
            u64 f = ((u64)h->rows[k].v[kk] * inv) % h->p;
            row_axpy(h, k, P, f, c);
            if (h->inconsistent && stop_on_inconsistent) return -1;
        }
        L->len = 0;
        /* retire column c */
        bucket_remove(h, c); h->ccnt[c] = 0; h->cstate[c] = 1; h->active_cols--;
        /* columns whose count dropped to 0 are free: retire them */
        /* (they stay in bucket 0; sweep bucket 0 lazily) */
        while (h->bhead[0] >= 0) { i32 z = h->bhead[0]; bucket_remove(h, z); h->cstate[z] = 2; h->active_cols--; }
    }
    return 4;
}

i32 sge_npiv(SGE *h) { return h->npiv; }
i32 sge_active_rows(SGE *h) { return h->active_rows; }
i32 sge_active_cols(SGE *h) { return h->active_cols; }
i64 sge_active_nnz(SGE *h) { return h->active_nnz; }
i64 sge_fill_events(SGE *h) { return h->fill_events; }
int sge_inconsistent(SGE *h) { return h->inconsistent; }

/* column states into out[ncols] (0 active, 1 pivot, 2 free/retired, 3 rhs); pivot columns in order into piv[npiv] */
void sge_col_states(SGE *h, unsigned char *out, i32 *piv) {
    memcpy(out, h->cstate, h->ncols);
    for (i32 i = 0; i < h->npiv; i++) piv[i] = h->rpiv[h->piv_order[i]];
}

/* export the active rows as CSR with ORIGINAL column ids.  First call with indptr only (nr+1 entries, needs nr = active_rows) */
void sge_export_active(SGE *h, i64 *indptr, i32 *indices, u32 *vals, i32 *rowids) {
    i64 pos = 0; i32 n = 0;
    for (i32 r = 0; r < h->nrows; r++) {
        if (h->rstate[r] != 1) continue;
        Row *R = &h->rows[r];
        if (rowids) rowids[n] = r;
        indptr[n] = pos;
        if (indices) { memcpy(indices + pos, R->c, sizeof(i32) * R->len); memcpy(vals + pos, R->v, sizeof(u32) * R->len); }
        pos += R->len; n++;
    }
    indptr[n] = pos;
}

/* Back-substitution: x[ncols] holds values of all non-pivot columns (and anything for rhs/pivot columns);
 * pivot unknowns are filled in reverse elimination order.  If use_rhs == 0 the right-hand side is taken as 0
 * (homogeneous solve, e.g. nullspace vectors). */
void sge_backsub(SGE *h, u32 *x, int use_rhs) {
    u64 p = h->p;
    for (i32 i = h->npiv - 1; i >= 0; i--) {
        i32 r = h->piv_order[i]; i32 j = h->rpiv[r]; Row *R = &h->rows[r];
        u64 s = 0, aj = 1, b = 0;
        for (i32 k = 0; k < R->len; k++) {
            i32 c = R->c[k];
            if (c == j) { aj = R->v[k]; continue; }
            if (c == h->rhs_col) { b = R->v[k]; continue; }
            s = (s + (u64)R->v[k] * (u64)x[c]) % p;
        }
        if (!use_rhs) b = 0;
        u64 val = ((b + p - s) % p) * inv_mod(aj, p) % p;
        x[j] = (u32)val;
    }
}

/* Batched back-substitution for nvec vectors stored row-major X[nvec][ncols]. */
void sge_backsub_many(SGE *h, u32 *X, i32 nvec, int use_rhs) {
    for (i32 v = 0; v < nvec; v++) sge_backsub(h, X + (i64)v * h->ncols, use_rhs);
}
