#ifndef HE_UTILS_H
#define HE_UTILS_H

#include <vector>
#include <assert.h>
#include "CtxtBit.h"

float clock_diff(clock_t &t1, clock_t &t2);

/*************************************************************************************/
/* these are used to evaluate the "comparison" */
CtxtBit compute_z (int i, int j, vector<CtxtBit>& ct_x, vector<CtxtBit>& ct_y);
CtxtBit compute_t (int i, int j, vector<CtxtBit>& ct_x, vector<CtxtBit>& ct_y);
CtxtBit compute_s (int i, int j, vector<CtxtBit>& ct_x, vector<CtxtBit>& ct_y);


/*************************************************************************************/
/* comparisons evaluations */
void evaluate_X_gt_Y (vector<CtxtBit>& ct_x, vector<CtxtBit>& ct_y, int t_bits);
void evaluate_X_ge_Y (vector<CtxtBit>& ct_x, vector<CtxtBit>& ct_y, int t_bits);
void evaluate_X_eq_Y (vector<CtxtBit>& ct_x, vector<CtxtBit>& ct_y, int t_bits);

/*************************************************************************************/
/* these are used to evaluate the "maximum" */
vector<CtxtBit> select (CtxtBit& c, vector<CtxtBit>& a, vector<CtxtBit>& b);
vector<CtxtBit> getmax (vector<vector<CtxtBit> >& vvct);                   // consumes a large number of levels 
vector<CtxtBit> getmax (vector<vector<CtxtBit> >& vvct, int start, int n); // A tree approach with a small nr. of levels consumed

/* ciphertext maintenance */
/*void batchRecrypt (vector<CtxtBit>& vct);
void ctxtRecrypt (CtxtBit& ct);*/

// test scope
int gmax (int v[], int start, int size);
int gmax (vector<int>& v, int start, int n);

#endif // HE_UTILS_H
