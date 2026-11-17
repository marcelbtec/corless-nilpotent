// Enumerate all sign sequences satisfying the mod-4 constraints (affine space t0 + span(basis)) and print
// those whose path matching polynomial D(x) = sum_M prod s_e x^|M| equals 1 exactly (= nilpotent matrix).
// Usage: ./enum params_K.txt slice nslices
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned __int128 u128;
static const uint64_t P = (1ULL << 61) - 1;
static inline uint64_t mulmod(uint64_t a, uint64_t b) { u128 z = (u128)a * b; uint64_t lo = (uint64_t)(z & P), hi = (uint64_t)(z >> 61); uint64_t r = lo + hi; return r >= P ? r - P : r; }
static inline uint64_t addmod(uint64_t a, uint64_t b) { uint64_t r = a + b; return r >= P ? r - P : r; }
static int m, n; static uint64_t t0, basis[64];
// D(P_k) = D(P_{k-1}) + s_{k-1} x D(P_{k-2}), D(P_0) = D(P_1) = 1; evaluated at x mod P
static inline uint64_t evalD(uint64_t t, uint64_t x) {
  uint64_t a = 1, b = 1; // a = D(P_{k-2}), b = D(P_{k-1})
  for (int k = 2; k <= m; k++) {
    uint64_t term = mulmod(x, a);
    uint64_t c = ((t >> (k - 2)) & 1) ? addmod(b, P - term) : addmod(b, term);
    a = b; b = c;
  }
  return b;
}
static int exact(uint64_t t) { // full integer check of all coefficients
  static long long A[64], B[64], C[64];
  for (int j = 0; j < 64; j++) A[j] = B[j] = 0;
  A[0] = 1; B[0] = 1;
  for (int k = 2; k <= m; k++) {
    int s = ((t >> (k - 2)) & 1) ? -1 : 1;
    C[0] = B[0];
    for (int j = 1; j <= k / 2; j++) C[j] = B[j] + s * A[j - 1];
    for (int j = 0; j <= k / 2; j++) { A[j] = B[j]; B[j] = C[j]; }
    for (int j = k / 2 + 1; j < 64; j++) { A[j] = B[j]; B[j] = 0; }
  }
  if (B[0] != 1) return 0;
  for (int j = 1; j <= m / 2; j++) if (B[j] != 0) return 0;
  return 1;
}
int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "r"); int slice = atoi(argv[2]), ns = atoi(argv[3]);
  if (fscanf(f, "%d %d", &m, &n) != 2) return 1;
  if (fscanf(f, "%llx", (unsigned long long *)&t0) != 1) return 1;
  for (int i = 0; i < n; i++) if (fscanf(f, "%llx", (unsigned long long *)&basis[i]) != 1) return 1;
  uint64_t total = 1ULL << n, per = (total + ns - 1) / ns, lo = per * slice < total ? per * slice : total, hi = lo + per < total ? lo + per : total;
  uint64_t x1 = 1234567891011ULL % P, x2 = 987654321987ULL % P;
  // start state: gray(lo)
  uint64_t g = lo ^ (lo >> 1), t = t0;
  for (int i = 0; i < n; i++) if ((g >> i) & 1) t ^= basis[i];
  uint64_t found = 0, filt = 0;
  for (uint64_t i = lo; i < hi; i++) {
    if (i > lo) { int b = __builtin_ctzll(i); t ^= basis[b]; }
    if (evalD(t, x1) == 1 && evalD(t, x2) == 1) {
      filt++;
      if (exact(t)) { found++; printf("%llx\n", (unsigned long long)t); }
    }
  }
  fprintf(stderr, "slice %d: candidates %llu, passed filter %llu, nilpotent %llu\n", slice, (unsigned long long)(hi - lo), (unsigned long long)filt, (unsigned long long)found);
  return 0;
}
