#include "curve_extras.h"
#include "tedwards.h"
#include <ec_params.h>
#include <assert.h>


bool ec_is_zero(ec_point_t const* P)
{
    return fp2_is_zero(&P->z);
}

void ec_init(ec_point_t* P)
{ // Initialize point as identity element (1:0)
    fp_t one = {0};
    
    memset((digit_t*)P, 0, NWORDS_FIELD*RADIX*4/8);
    one[0] = 1;
    fp_tomont(P->x.re, one);
}

void xDBL(ec_point_t* Q, ec_point_t const* P, ec_point_t const* AC)
{
    // This version computes the coefficient values A+2C and 4C on-the-fly 
    // The curve coefficients are passed via AC = (A:C)
    fp2_t t0, t1, t2, t3;

    fp2_add(&t0, &P->x, &P->z);
    fp2_sqr(&t0, &t0);
    fp2_sub(&t1, &P->x, &P->z);
    fp2_sqr(&t1, &t1);
    fp2_sub(&t2, &t0, &t1);
    fp2_add(&t3, &AC->z, &AC->z);  
    fp2_mul(&t1, &t1, &t3);
    fp2_add(&t1, &t1, &t1);
    fp2_mul(&Q->x, &t0, &t1);
    fp2_add(&t0, &t3, &AC->x);
    fp2_mul(&t0, &t0, &t2);
    fp2_add(&t0, &t0, &t1);
    fp2_mul(&Q->z, &t0, &t2);
}

void xDBLv2(ec_point_t* Q, ec_point_t const* P, ec_point_t const* A24)
{
    // This version receives the coefficient value A24 = (A+2C:4C) 
    fp2_t t0, t1, t2;

    fp2_add(&t0, &P->x, &P->z);
    fp2_sqr(&t0, &t0);
    fp2_sub(&t1, &P->x, &P->z);
    fp2_sqr(&t1, &t1);
    fp2_sub(&t2, &t0, &t1);
    fp2_mul(&t1, &t1, &A24->z);
    fp2_mul(&Q->x, &t0, &t1);
    fp2_mul(&t0, &t2, &A24->x);
    fp2_add(&t0, &t0, &t1);
    fp2_mul(&Q->z, &t0, &t2);
}

void xADD(ec_point_t* R, ec_point_t const* P, ec_point_t const* Q, ec_point_t const* PQ)
{
    fp2_t t0, t1, t2, t3;

    fp2_add(&t0, &P->x, &P->z);
    fp2_sub(&t1, &P->x, &P->z);
    fp2_add(&t2, &Q->x, &Q->z);
    fp2_sub(&t3, &Q->x, &Q->z);
    fp2_mul(&t0, &t0, &t3);
    fp2_mul(&t1, &t1, &t2);
    fp2_add(&t2, &t0, &t1);
    fp2_sub(&t3, &t0, &t1);
    fp2_sqr(&t2, &t2);
    fp2_sqr(&t3, &t3);
    fp2_mul(&t2, &PQ->z, &t2);
    fp2_mul(&R->z, &PQ->x, &t3);
    fp2_copy(&R->x, &t2);
}

void xDBLADD(ec_point_t* R, ec_point_t* S, ec_point_t const* P, ec_point_t const* Q, ec_point_t const* PQ, ec_point_t const* A24)
{
    // Requires precomputation of A24 = (A+2C:4C)
    fp2_t t0, t1, t2;

    fp2_add(&t0, &P->x, &P->z);
    fp2_sub(&t1, &P->x, &P->z);
    fp2_sqr(&R->x, &t0);
    fp2_sub(&t2, &Q->x, &Q->z);
    fp2_add(&S->x, &Q->x, &Q->z);
    fp2_mul(&t0, &t0, &t2);
    fp2_sqr(&R->z, &t1);
    fp2_mul(&t1, &t1, &S->x);
    fp2_sub(&t2, &R->x, &R->z);
    fp2_mul(&R->z, &R->z, &A24->z);
    fp2_mul(&R->x, &R->x, &R->z);
    fp2_mul(&S->x, &A24->x, &t2);
    fp2_sub(&S->z, &t0, &t1);
    fp2_add(&R->z, &R->z, &S->x);
    fp2_add(&S->x, &t0, &t1);
    fp2_mul(&R->z, &R->z, &t2);
    fp2_sqr(&S->z, &S->z);
    fp2_sqr(&S->x, &S->x);
    fp2_mul(&S->z, &S->z, &PQ->x);
    fp2_mul(&S->x, &S->x, &PQ->z);
}

bool is_point_equal(const ec_point_t* P, const ec_point_t* Q)
{ // Evaluate if two points in Montgomery coordinates (X:Z) are equal
  // Returns 1 (true) if P=Q, 0 (false) otherwise
    fp2_t t0, t1;

    if ((fp2_is_zero(&P->x) && fp2_is_zero(&P->z)) || (fp2_is_zero(&Q->x) && fp2_is_zero(&Q->z))) {
        return fp2_is_zero(&P->x) && fp2_is_zero(&P->z) && fp2_is_zero(&Q->x) && fp2_is_zero(&Q->z);
    }

    fp2_mul(&t0, &P->x, &Q->z);
    fp2_mul(&t1, &Q->x, &P->z);
    fp2_sub(&t0, &t0, &t1);

    return fp2_is_zero(&t0);
}

void swap_points(ec_point_t* P, ec_point_t* Q, const digit_t option)
{ // Swap points
  // If option = 0 then P <- P and Q <- Q, else if option = 0xFF...FF then P <- Q and Q <- P
    digit_t temp;

    for (int i = 0; i < NWORDS_FIELD; i++) {
        temp = option & (P->x.re[i] ^ Q->x.re[i]);
        P->x.re[i] = temp ^ P->x.re[i];
        Q->x.re[i] = temp ^ Q->x.re[i];
        temp = option & (P->x.im[i] ^ Q->x.im[i]);
        P->x.im[i] = temp ^ P->x.im[i];
        Q->x.im[i] = temp ^ Q->x.im[i];
        temp = option & (P->z.re[i] ^ Q->z.re[i]);
        P->z.re[i] = temp ^ P->z.re[i];
        Q->z.re[i] = temp ^ Q->z.re[i];
        temp = option & (P->z.im[i] ^ Q->z.im[i]);
        P->z.im[i] = temp ^ P->z.im[i];
        Q->z.im[i] = temp ^ Q->z.im[i];
    }
}

void copy_point(ec_point_t* P, ec_point_t const* Q)
{
    fp2_copy(&(P->x), &(Q->x));
    fp2_copy(&(P->z), &(Q->z));
}

void ec_normalize(ec_point_t* P){
    fp2_inv(&P->z);
    fp2_mul(&P->x, &P->x, &P->z);
    fp_mont_setone(P->z.re);
    fp_set(P->z.im, 0);
}

void ec_neg(ec_point_t* res, const ec_point_t* P){
    // DOES NOTHING
    copy_point(res, P);
}

void xMUL(ec_point_t* Q, ec_point_t const* P, digit_t const* k, ec_curve_t const* curve)
{
    ec_point_t R0, R1, A24;
    digit_t mask;
    unsigned int bit = 0, prevbit = 0, swap;

    fp2_add(&A24.x, &curve->C, &curve->C);    // Precomputation of A24=(A+2C:4C)
    fp2_add(&A24.z, &A24.x, &A24.x);
    fp2_add(&A24.x, &A24.x, &curve->A);

    // R0 <- (1:0), R1 <- P
    ec_init(&R0);
    fp2_copy(&R1.x, &P->x);
    fp2_copy(&R1.z, &P->z);

    // Main loop
    for (int i = BITS-1; i >= 0; i--) {                                          
        bit = (k[i >> LOG2RADIX] >> (i & (RADIX-1))) & 1;                         
        swap = bit ^ prevbit;
        prevbit = bit;
        mask = 0 - (digit_t)swap;

        swap_points(&R0, &R1, mask);
        xDBLADD(&R0, &R1, &R0, &R1, P, &A24);
    }
    swap = 0 ^ prevbit;
    mask = 0 - (digit_t)swap;
    swap_points(&R0, &R1, mask);

    fp2_copy(&Q->x, &R0.x);
    fp2_copy(&Q->z, &R0.z);
}

void xMULv2(ec_point_t* Q, ec_point_t const* P, digit_t const* k, const int kbits, ec_point_t const* A24)
{
    // This version receives the coefficient value A24 = (A+2C:4C) 
    ec_point_t R0, R1;
    digit_t mask;
    unsigned int bit = 0, prevbit = 0, swap;

    // R0 <- (1:0), R1 <- P
    ec_init(&R0);
    fp2_copy(&R1.x, &P->x);
    fp2_copy(&R1.z, &P->z);

    // Main loop
    for (int i = kbits-1; i >= 0; i--) {                              
        bit = (k[i >> LOG2RADIX] >> (i & (RADIX-1))) & 1;                      
        swap = bit ^ prevbit;
        prevbit = bit;
        mask = 0 - (digit_t)swap;

        swap_points(&R0, &R1, mask);
        xDBLADD(&R0, &R1, &R0, &R1, P, A24);
    }
    swap = 0 ^ prevbit;
    mask = 0 - (digit_t)swap;
    swap_points(&R0, &R1, mask);

    fp2_copy(&Q->x, &R0.x);
    fp2_copy(&Q->z, &R0.z);
}

static void mp_add(digit_t* c, const digit_t* a, const digit_t* b, const unsigned int nwords)
{ // Multiprecision addition
    unsigned int i, carry = 0;

    for (i = 0; i < nwords; i++) {
        ADDC(c[i], carry, a[i], b[i], carry);
    }
}

static void mp_sub(digit_t* c, digit_t const* a, digit_t const* b, const unsigned int nwords)
{ // Multiprecision subtraction, assuming a > b
    unsigned int i, borrow = 0;

    for (i = 0; i < nwords; i++) {
        SUBC(c[i], borrow, a[i], b[i], borrow);
    }
}

void select_ct(digit_t* c, const digit_t* a, const digit_t* b, const digit_t mask, const int nwords)
{ // Select c <- a if mask = 0, select c <- b if mask = 1...1

    for (int i = 0; i < nwords; i++) {
        c[i] = ((a[i] ^ b[i]) & mask) ^ a[i];
    }
}

void swap_ct(digit_t* a, digit_t* b, const digit_t option, const int nwords)
{ // Swap entries
  // If option = 0 then P <- P and Q <- Q, else if option = 0xFF...FF then a <- b and b <- a
    digit_t temp;

    for (int i = 0; i < nwords; i++) {
        temp = option & (a[i] ^ b[i]);
        a[i] = temp ^ a[i];
        b[i] = temp ^ b[i];
    }
}

// Compute S = k*P + l*Q, with PQ = P+Q
void xDBLMUL(ec_point_t* S, ec_point_t const* P, digit_t const* k, ec_point_t const* Q, digit_t const* l, ec_point_t const* PQ, ec_curve_t const* curve)
{
    int i;
    digit_t evens, mevens, bitk0, bitl0, maskk, maskl, temp, bs1_ip1, bs2_ip1, bs1_i, bs2_i, h;
    digit_t sigma[2] = {0}, pre_sigma = 0;
    digit_t k_t[NWORDS_ORDER], l_t[NWORDS_ORDER], one[NWORDS_ORDER] = {0}, r[2*BITS] = {0};            
    ec_point_t A24, DIFF1a, DIFF1b, DIFF2a, DIFF2b, R[3] = {0}, T[3];

    // Derive sigma according to parity
    bitk0 = (k[0] & 1);
    bitl0 = (l[0] & 1);
    maskk = 0 - bitk0;               // Parity masks: 0 if even, otherwise 1...1
    maskl = 0 - bitl0;
    sigma[0] = (bitk0 ^ 1);
    sigma[1] = (bitl0 ^ 1);
    evens = sigma[0] + sigma[1];     // Count number of even scalars
    mevens = 0 - (evens & 1);        // Mask mevens <- 0 if # even scalars = 0 or 2, otherwise mevens = 1...1

    // If k and l are both even or both odd, pick sigma = (0,1)
    sigma[0] = (sigma[0] & mevens);
    sigma[1] = (sigma[1] & mevens) | (1 & ~mevens);

    // Convert even scalars to odd
    one[0] = 1;
    mp_sub(k_t, k, one, NWORDS_ORDER);
    mp_sub(l_t, l, one, NWORDS_ORDER);
    select_ct(k_t, k_t, k, maskk, NWORDS_ORDER);                                        
    select_ct(l_t, l_t, l, maskl, NWORDS_ORDER);

    // Scalar recoding
    for (i = 0; i < BITS; i++) {
        // If sigma[0] = 1 swap k_t and l_t
        maskk = 0 - (sigma[0] ^ pre_sigma);
        swap_ct(k_t, l_t, maskk, NWORDS_ORDER);
        
        if (i == BITS-1) {
            bs1_ip1 = 0;
            bs2_ip1 = 0;
        } else {
            bs1_ip1 = mp_shiftr(k_t, 1, NWORDS_ORDER);
            bs2_ip1 = mp_shiftr(l_t, 1, NWORDS_ORDER);
        }
        bs1_i = k_t[0] & 1;
        bs2_i = l_t[0] & 1;

        r[2*i]   = bs1_i ^ bs1_ip1;
        r[2*i+1] = bs2_i ^ bs2_ip1;

        // Revert sigma if second bit, r_(2i+1), is 1 
        pre_sigma = sigma[0];
        maskk = 0 - r[2*i+1];
        select_ct(&temp, &sigma[0], &sigma[1], maskk, 1);
        select_ct(&sigma[1], &sigma[1], &sigma[0], maskk, 1);
        sigma[0] = temp;
    }

    // Point initialization
    ec_init(&R[0]);
    maskk = 0 - sigma[0];
    select_ct((digit_t*)&R[1], (digit_t*)P, (digit_t*)Q, maskk, 4*NWORDS_FIELD);
    select_ct((digit_t*)&R[2], (digit_t*)Q, (digit_t*)P, maskk, 4*NWORDS_FIELD);
    fp2_copy(&DIFF1a.x, &R[1].x);
    fp2_copy(&DIFF1a.z, &R[1].z);
    fp2_copy(&DIFF1b.x, &R[2].x);
    fp2_copy(&DIFF1b.z, &R[2].z);

    // Initialize DIFF2a <- P+Q, DIFF2b <- P-Q
    xADD(&R[2], &R[1], &R[2], PQ);
    fp2_copy(&DIFF2a.x, &R[2].x);
    fp2_copy(&DIFF2a.z, &R[2].z);
    fp2_copy(&DIFF2b.x, &PQ->x);
    fp2_copy(&DIFF2b.z, &PQ->z);

    fp2_add(&A24.x, &curve->C, &curve->C);    // Precomputation of A24=(A+2C:4C)
    fp2_add(&A24.z, &A24.x, &A24.x);
    fp2_add(&A24.x, &A24.x, &curve->A);

    // Main loop
    for (i = BITS-1; i>=0; i--) {
        h = r[2*i] + r[2*i+1];    // in {0, 1, 2}
        maskk = 0 - (h & 1);
        select_ct((digit_t*)&T[0], (digit_t*)&R[0], (digit_t*)&R[1], maskk, 4*NWORDS_FIELD);
        maskk = 0 - (h >> 1);
        select_ct((digit_t*)&T[0], (digit_t*)&T[0], (digit_t*)&R[2], maskk, 4*NWORDS_FIELD);
        xDBLv2(&T[0], &T[0], &A24);

        maskk = 0 - r[2*i+1];     // in {0, 1}
        select_ct((digit_t*)&T[1], (digit_t*)&R[0], (digit_t*)&R[1], maskk, 4*NWORDS_FIELD);
        select_ct((digit_t*)&T[2], (digit_t*)&R[1], (digit_t*)&R[2], maskk, 4*NWORDS_FIELD);
        swap_points(&DIFF1a, &DIFF1b, maskk);
        xADD(&T[1], &T[1], &T[2], &DIFF1a);
        xADD(&T[2], &R[0], &R[2], &DIFF2a);
        
        // If hw (mod 2) = 1 then swap DIFF2a and DIFF2b
        maskk = 0 - (h & 1);
        swap_points(&DIFF2a, &DIFF2b, maskk);

        // R <- T
        memcpy((digit_t*)&R[0], (digit_t*)&T[0], NWORDS_FIELD*RADIX*4/8);
        memcpy((digit_t*)&R[1], (digit_t*)&T[1], NWORDS_FIELD*RADIX*4/8);
        memcpy((digit_t*)&R[2], (digit_t*)&T[2], NWORDS_FIELD*RADIX*4/8);
    }

    // Output R[evens]
    select_ct((digit_t*)S, (digit_t*)&R[0], (digit_t*)&R[1], mevens, 4*NWORDS_FIELD);
    maskk = 0 - (bitk0 & bitl0);
    select_ct((digit_t*)S, (digit_t*)S, (digit_t*)&R[2], maskk, 4*NWORDS_FIELD);
}

void ec_ladder3pt(ec_point_t *R, fp_t const m, ec_point_t const *P, ec_point_t const *Q, ec_point_t const *PQ, ec_curve_t const *A)
{
    // Curve constant in the form A24=(A+2C:4C)
    ec_point_t A24;
    fp2_add(&A24.z, &A->C, &A->C);
    fp2_add(&A24.x, &A->A, &A24.z);
    fp2_add(&A24.z, &A24.z, &A24.z);

	ec_point_t X0, X1, X2;
	copy_point(&X0, Q);
	copy_point(&X1, P);
	copy_point(&X2, PQ);

	int i,j;
	uint64_t t;
	for (i = 0; i < NWORDS_FIELD; i++)
	{
		t = 1;
		for (j = 0 ; j < 64; j++)
		{
			swap_points(&X1, &X2, -((t & m[i]) == 0));
			xDBLADD(&X0, &X1, &X0, &X1, &X2, &A24);
			swap_points(&X1, &X2, -((t & m[i]) == 0));
			t <<= 1;
		};
	};
	copy_point(R, &X1);
}

void ec_j_inv(fp2_t* j_inv, const ec_curve_t* curve){
    /* j-invariant computation for montgommery coefficient A2=(A+2C:4C) */
    fp2_t t0, t1;
    
    fp2_sqr(&t1, &curve->C);
    fp2_sqr(j_inv, &curve->A);
    fp2_add(&t0, &t1, &t1);
    fp2_sub(&t0, j_inv, &t0);
    fp2_sub(&t0, &t0, &t1);
    fp2_sub(j_inv, &t0, &t1);
    fp2_sqr(&t1, &t1);
    fp2_mul(j_inv, j_inv, &t1);
    fp2_add(&t0, &t0, &t0);
    fp2_add(&t0, &t0, &t0);
    fp2_sqr(&t1, &t0);
    fp2_mul(&t0, &t0, &t1);
    fp2_add(&t0, &t0, &t0);
    fp2_add(&t0, &t0, &t0);
    fp2_inv(j_inv);
    fp2_mul(j_inv, &t0, j_inv);    
}

static void jac_init(jac_point_t* P)
{ // Initialize Montgomery in Jacobian coordinates as identity element (0:1:0)
    fp_t one = {0};

    memset((digit_t*)P, 0, NWORDS_FIELD*RADIX*6/8);
    one[0] = 1;
    fp_tomont(P->y.re, one);
}

static bool is_jac_equal(const jac_point_t* P, const jac_point_t* Q)
{ // Evaluate if two points in Jacobian coordinates (X:Y:Z) are equal
  // Returns 1 (true) if P=Q, 0 (false) otherwise
    fp2_t t0, t1, t2, t3;

    fp2_sqr(&t0, &Q->z);
    fp2_mul(&t2, &P->x, &t0);       // x1*z2^2
    fp2_sqr(&t1, &P->z);
    fp2_mul(&t3, &Q->x, &t1);       // x2*z1^2
    fp2_sub(&t2, &t2, &t3);

    fp2_mul(&t0, &t0, &Q->z);
    fp2_mul(&t0, &P->y, &t0);       // y1*z2^3
    fp2_mul(&t1, &t1, &P->z);
    fp2_mul(&t1, &Q->y, &t1);       // y2*z1^3
    fp2_sub(&t0, &t0, &t1);

    return fp2_is_zero(&t0) && fp2_is_zero(&t2);
}

static bool is_jac_xz_equal(const jac_point_t* P, const ec_point_t* Q)
{ // Evaluate if point P in Jacobian coordinates is equal to Q in homogeneous projective coordinates (X:Z) 
  // Comparison is up to sign (only compares X and Z coordinates)
  // Returns 1 (true) if P=Q, 0 (false) otherwise
    fp2_t t0, t1;

    fp2_mul(&t0, &P->x, &Q->z);     // x1*z2
    fp2_sqr(&t1, &P->z);
    fp2_mul(&t1, &Q->x, &t1);       // x2*z1^2
    fp2_sub(&t0, &t0, &t1);

    return fp2_is_zero(&t0);
}

static void copy_jac_point(jac_point_t* P, jac_point_t const* Q)
{
    fp2_copy(&(P->x), &(Q->x));
    fp2_copy(&(P->y), &(Q->y));
    fp2_copy(&(P->z), &(Q->z));
}

static void jac_neg(jac_point_t* Q, jac_point_t const* P)
{
    fp2_copy(&Q->x, &P->x);
    fp2_neg(&Q->y, &P->y);
    fp2_copy(&Q->z, &P->z);
}

void DBL(jac_point_t* Q, jac_point_t const* P, ec_curve_t const* AC)
{ // Doubling on a Montgomery curve, representation in Jacobian coordinates (X:Y:Z) corresponding to (X/Z^2,Y/Z^3) 
  // This version receives the coefficient value A 
    fp2_t t0, t1, t2, t3;

    if (fp2_is_zero(&P->x) && fp2_is_zero(&P->z)) {
        jac_init(Q);
        return;
    }

    fp2_sqr(&t0, &P->x);            // t0 = x1^2
    fp2_add(&t1, &t0, &t0);
    fp2_add(&t0, &t0, &t1);         // t0 = 3x1^2
    fp2_sqr(&t1, &P->z);            // t1 = z1^2
    fp2_mul(&t2, &P->x, &AC->A);
    fp2_add(&t2, &t2, &t2);         // t2 = 2Ax1  
    fp2_add(&t2, &t1, &t2);         // t2 = 2Ax1+z1^2
    fp2_mul(&t2, &t1, &t2);         // t2 = z1^2(2Ax1+z1^2)
    fp2_add(&t2, &t0, &t2);         // t2 = alpha = 3x1^2 + z1^2(2Ax1+z1^2)
    fp2_mul(&Q->z, &P->y, &P->z);
    fp2_add(&Q->z, &Q->z, &Q->z);   // z2 = 2y1z1
    fp2_sqr(&t0, &Q->z);
    fp2_mul(&t0, &t0, &AC->A);      // t0 = 4Ay1^2z1^2
    fp2_sqr(&t1, &P->y);
    fp2_add(&t1, &t1, &t1);         // t1 = 2y1^2
    fp2_add(&t3, &P->x, &P->x);     // t3 = 2x1
    fp2_mul(&t3, &t1, &t3);         // t3 = 4x1y1^2
    fp2_sqr(&Q->x, &t2);            // x2 = alpha^2
    fp2_sub(&Q->x, &Q->x, &t0);     // x2 = alpha^2 - 4Ay1^2z1^2
    fp2_sub(&Q->x, &Q->x, &t3);
    fp2_sub(&Q->x, &Q->x, &t3);     // x2 = alpha^2 - 4Ay1^2z1^2 - 8x1y1^2
    fp2_sub(&Q->y, &t3, &Q->x);     // y2 = 4x1y1^2 - x2
    fp2_mul(&Q->y, &Q->y, &t2);     // y2 = alpha(4x1y1^2 - x2)
    fp2_sqr(&t1, &t1);              // t1 = 4y1^4
    fp2_sub(&Q->y, &Q->y, &t1);
    fp2_sub(&Q->y, &Q->y, &t1);     // y2 = alpha(4x1y1^2 - x2) - 8y1^4
}

void ADD(jac_point_t* R, jac_point_t const* P, jac_point_t const* Q, ec_curve_t const* AC)
{ // Addition on a Montgomery curve, representation in Jacobian coordinates (X:Y:Z) corresponding to (X/Z^2,Y/Z^3)
  // This version receives the coefficient value A 
    fp2_t t0, t1, t2, t3, t4, t5, t6;
    jac_point_t T;

    if (is_jac_equal(P, Q)) {
        DBL(R, P, AC);
        return;
    }
    jac_neg(&T, P);
    if (is_jac_equal(&T, Q)) {
        jac_init(R);
        return;
    }
    if (fp2_is_zero(&P->x) && fp2_is_zero(&P->z)) {
        copy_jac_point(R, Q);
        return;
    } else if (fp2_is_zero(&Q->x) && fp2_is_zero(&Q->z)) {
        copy_jac_point(R, P);
        return;
    }

    fp2_sqr(&t0, &P->z);            // t0 = z1^2
    fp2_mul(&t1, &t0, &P->z);       // t1 = z1^3
    fp2_sqr(&t2, &Q->z);            // t2 = z2^2
    fp2_mul(&t3, &t2, &Q->z);       // t3 = z2^3
    fp2_mul(&t1, &t1, &Q->y);       // t1 = y2z1^3
    fp2_mul(&t3, &t3, &P->y);       // t3 = y1z2^3
    fp2_sub(&t1, &t1, &t3);         // t1 = lambda1 = y2z1^3 - y1z2^3
    fp2_mul(&t0, &t0, &Q->x);       // t0 = x2z1^2
    fp2_mul(&t2, &t2, &P->x);       // t2 = x1z2^2
    fp2_sub(&t4, &t0, &t2);         // t4 = lambda3 = x2z1^2 - x1z2^2
    fp2_add(&t0, &t0, &t2);         // t0 = lambda2 = x2z1^2 + x1z2^2
    fp2_mul(&t5, &P->z, &Q->z);     // t5 = z1z2
    fp2_mul(&R->z, &t4, &t5);       // z3 = z1z2*lambda3
    fp2_sqr(&t5, &t5);              // t5 = z1^2z2^2
    fp2_mul(&t5, &AC->A, &t5);      // t5 = Az1^2z2^2
    fp2_add(&t0, &t0, &t5);         // t0 = Az1^2z2^2 + lambda2
    fp2_sqr(&t6, &t4);              // t6 = lambda3^2
    fp2_mul(&t5, &t0, &t6);         // t5 = lambda3^2(Az1^2z2^2 + lambda2)
    fp2_sqr(&R->x, &t1);            // x3 = lambda1^2
    fp2_sub(&R->x, &R->x, &t5);     // x3 = lambda1^2 - lambda3^2(Az1^2z2^2 + lambda2)
    fp2_mul(&t3, &t3, &t4);         // t3 = y1z2^3*lambda3
    fp2_mul(&t3, &t3, &t6);         // t3 = y1z2^3*lambda3^3
    fp2_mul(&t2, &t2, &t6);         // t2 = x1z2^2*lambda3^2
    fp2_sub(&R->y, &t2, &R->x);     // y3 = x1z2^2*lambda3^2 - x3
    fp2_mul(&R->y, &R->y, &t1);     // y3 = lambda1(x1z2^2*lambda3^2 - x3)
    fp2_sub(&R->y, &R->y, &t3);     // y3 = lambda1(x1z2^2*lambda3^2 - x3) - y1z2^3*lambda3^3
}

void TPL(jac_point_t* Q, jac_point_t const* P, ec_curve_t const* AC)
{ // Naive tripling on a Montgomery curve, representation in Jacobian coordinates (X:Y:Z) corresponding to (X/Z^2,Y/Z^3) 
  // This version receives the coefficient value A 
    jac_point_t R;

    DBL(&R, P, AC);
    ADD(Q, &R, P, AC);
}

void recover_y(fp2_t* y, fp2_t const* Px, ec_curve_t const* curve)
{ // Recover y-coordinate of a point on the Montgomery curve y^2 = x^3 + Ax^2 + x
    fp2_t t0;

    fp2_sqr(&t0, Px);
    fp2_mul(y, &t0, &curve->A);    // Ax^2
    fp2_add(y, y, Px);             // Ax^2 + x
    fp2_mul(&t0, &t0, Px);
    fp2_add(y, y, &t0);            // x^3 + Ax^2 + x
    fp2_sqrt(y);

    //fp2_t t0, t1;

    //fp2_sqr(&t0, &P->x);
    //fp2_sqr(&t1, &P->z);
    //fp2_mul(y, &t0, &t1);
    //fp2_mul(y, y, &curve->A);    // AX^2Z^2
    //fp2_mul(&t0, &t0, &P->x);
    //fp2_add(y, y, &t0);          // X^3 + AX^2Z^2
    //fp2_mul(&t0, &t1, &P->z);
    //fp2_sqr(&t1, &t0);
    //fp2_mul(&t1, &t1, &P->x);
    //fp2_add(y, y, &t1);          // X^3 + AX^2Z^2 + XZ^6
    //fp2_sqrt(y);
}


static int mp_compare(digit_t* a, digit_t* b, unsigned int nwords)
{ // Multiprecision comparison, a=b? : (1) a>b, (0) a=b, (-1) a<b

    for (int i = nwords-1; i >= 0; i--) {
        if (a[i] > b[i]) return 1;
        else if (a[i] < b[i]) return -1;
    }
    return 0;
}

static bool mp_is_zero(const digit_t* a, unsigned int nwords)
{ // Is a multiprecision element zero?
  // Returns 1 (true) if a=0, 0 (false) otherwise
    digit_t r = 0;

    for (unsigned int i = 0; i < nwords; i++)
        r |= a[i] ^ 0;

    return (bool)is_digit_zero_ct(r);
}

void DBLMUL(jac_point_t* R, const jac_point_t* P, const digit_t k, const jac_point_t* Q, const digit_t l, const ec_curve_t* curve)
{  // Double-scalar multiplication R <- k*P + l*Q, fixed for 64-bit scalars
    digit_t k_t, l_t;
    jac_point_t PQ;

    ADD(&PQ, P, Q, curve);
    jac_init(R);

    for (int i = 0; i < 64; i++) {
        k_t = k >> (63-i);
        k_t &= 0x01;
        l_t = l >> (63-i);
        l_t &= 0x01;
        DBL(R, R, curve);
        if (k_t == 1 && l_t == 1) {
            ADD(R, R, &PQ, curve);
        } else if (k_t == 1) {
            ADD(R, R, P, curve);
        } else if (l_t == 1) {
            ADD(R, R, Q, curve);
        }
    }
}

#define SCALAR_BITS 128  //// TODO: this could be defined by level

void DBLMUL2(jac_point_t* R, const jac_point_t* P, const digit_t* k, const jac_point_t* Q, const digit_t* l, const ec_curve_t* curve)
{  // Double-scalar multiplication R <- k*P + l*Q, fixed for 128-bit scalars
    digit_t k_t[2], l_t[2];
    jac_point_t PQ;

    ADD(&PQ, P, Q, curve);
    jac_init(R);

    for (int i = 0; i < SCALAR_BITS-64; i++) {
        k_t[1] = k[1] >> (SCALAR_BITS-64-1-i);
        k_t[1] &= 0x01;
        l_t[1] = l[1] >> (SCALAR_BITS-64-1-i);
        l_t[1] &= 0x01;
        DBL(R, R, curve);
        if (k_t[1] == 1 && l_t[1] == 1) {
            ADD(R, R, &PQ, curve);
        } else if (k_t[1] == 1) {
            ADD(R, R, P, curve);
        } else if (l_t[1] == 1) {
            ADD(R, R, Q, curve);
        }
    }

    for (int i = 0; i < 64; i++) {
        k_t[0] = k[0] >> (64-1-i);
        k_t[0] &= 0x01;
        l_t[0] = l[0] >> (64-1-i);
        l_t[0] &= 0x01;
        DBL(R, R, curve);
        if (k_t[0] == 1 && l_t[0] == 1) {
            ADD(R, R, &PQ, curve);
        } else if (k_t[0] == 1) {
            ADD(R, R, P, curve);
        } else if (l_t[0] == 1) {
            ADD(R, R, Q, curve);
        }
    }
}

static void mp_mul2(digit_t* c, const digit_t* a, const digit_t* b)
{ // Multiprecision multiplication fixed to two-digit operands
    unsigned int carry = 0;
    digit_t t0[2], t1[2], t2[2];

    MUL(t0, a[0], b[0]);
    MUL(t1, a[0], b[1]);
    ADDC(t0[1], carry, t0[1], t1[0], carry);
    ADDC(t1[1], carry, 0, t1[1], carry);
    MUL(t2, a[1], b[1]);
    ADDC(t2[0], carry, t2[0], t1[1], carry);
    ADDC(t2[1], carry, 0, t2[1], carry);
    c[0] = t0[0];
    c[1] = t0[1];
    c[2] = t2[0];
    c[3] = t2[1];
}

static void ec_dlog_2_step(digit_t* x, digit_t* y, const jac_point_t* R, const int f, const int B, const jac_point_t* Pe2, const jac_point_t* Qe2, const jac_point_t* PQe2, const ec_curve_t* curve)
{ // Based on Montgomery formulas using Jacobian coordinates
    int i, j;
    digit_t value = 1;
    jac_point_t P, Q, TT, SS, PQp, T[(POWER_OF_2-1)/2], Re[(POWER_OF_2-1)/2];    // Storage could be reduced to e1 points

    *x = 0;        
    *y = 0;

    copy_jac_point(&P, &Pe2[f-1]);
    copy_jac_point(&Q, &Qe2[f-1]);
    copy_jac_point(&Re[f-1], R);

    for (i = 0; i < (f-1); i++) {
        DBL(&Re[f-i-2], &Re[f-i-1], curve);
    }
    ADD(&PQp, &P, &Q, curve);

    // Unrolling the first two iterations
    if (is_jac_equal(&Pe2[0], &Re[0])) {
        *x += 1;
        copy_jac_point(&T[0], &P);
        for (j = 3; j <= B; j++) {
            copy_jac_point(&T[j-1], &Pe2[j-1]);
        }
        jac_neg(&TT, &Pe2[1]);
        ADD(&Re[0], &Re[1], &TT, curve);
    } else if (is_jac_equal(&Qe2[0], &Re[0])) {
        *y += 1;
        copy_jac_point(&T[0], &Q);
        for (j = 3; j <= B; j++) {
            copy_jac_point(&T[j-1], &Qe2[j-1]);
        }
        jac_neg(&TT, &Qe2[1]);
        ADD(&Re[0], &Re[1], &TT, curve);
    } else if (is_jac_equal(&PQe2[0], &Re[0])) {
        *x += 1;
        *y += 1;
        copy_jac_point(&T[0], &PQp);
        for (j = 3; j <= B; j++) {
            copy_jac_point(&T[j-1], &PQe2[j-1]);
        }
        jac_neg(&TT, &PQe2[1]);
        ADD(&Re[0], &Re[1], &TT, curve);
    } else {
        jac_init(&T[0]);
        for (j = 3; j <= B; j++) {
            jac_init(&T[j-1]);
        }
        copy_jac_point(&Re[0], &Re[1]);
    }

    // Unrolling iterations 3-B
    for (i = 3; i <= B; i++) {
        value <<= 1;
        jac_neg(&TT, &T[i-1]);
        ADD(&TT, &Re[i-1], &TT, curve);       // TT <- Re[i-1]-T[i-1]
        if (is_jac_equal(&Pe2[0], &Re[0])) {
            *x += value;
            ADD(&T[0], &T[0], &Pe2[f-i+1], curve);
            for (j = i+1; j <= B; j++) {
                ADD(&T[j-1], &T[j-1], &Pe2[j-i+1], curve);
            }
            jac_neg(&SS, &Pe2[1]);
            ADD(&Re[0], &TT, &SS, curve);
        } else if (is_jac_equal(&Qe2[0], &Re[0])) {
            *y += value;
            ADD(&T[0], &T[0], &Qe2[f-i+1], curve);
            for (j = i+1; j <= B; j++) {
                ADD(&T[j-1], &T[j-1], &Qe2[j-i+1], curve);
            }
            jac_neg(&SS, &Qe2[1]);
            ADD(&Re[0], &TT, &SS, curve);
        } else if (is_jac_equal(&PQe2[0], &Re[0])) {
            *x += value;
            *y += value;
            ADD(&T[0], &T[0], &PQe2[f-i+1], curve);
            for (j = i+1; j <= B; j++) {
                ADD(&T[j-1], &T[j-1], &PQe2[j-i+1], curve);
            }
            jac_neg(&SS, &PQe2[1]);
            ADD(&Re[0], &TT, &SS, curve);
        } else {
            copy_jac_point(&Re[0], &TT);
        }
    }

    // Main Loop
    for (i = B; i < f; i++) {
        value <<= 1;
        if (is_jac_equal(&Pe2[0], &Re[0])) {
            *x += value;
            ADD(&T[0], &T[0], &Pe2[f-i], curve);
        } else if (is_jac_equal(&Qe2[0], &Re[0])) {
            *y += value;
            ADD(&T[0], &T[0], &Qe2[f-i], curve);
        } else if (is_jac_equal(&PQe2[0], &Re[0])) {
            *x += value;
            *y += value;
            ADD(&T[0], &T[0], &PQe2[f-i], curve);
        }
        jac_neg(&TT, &T[0]);
        ADD(&Re[0], R, &TT, curve);       // TT <- R-T[0]
        for (j = 0; j < (f-i-1); j++) {
            DBL(&Re[0], &Re[0], curve);
        }
    }
    value <<= 1;
    if (is_jac_equal(&Pe2[0], &Re[0])) {
        *x += value;
    } else if (is_jac_equal(&Qe2[0], &Re[0])) {
        *y += value;
    } else if (is_jac_equal(&PQe2[0], &Re[0])) {
        *x += value;
        *y += value;
    }
}

void ec_dlog_2(digit_t* scalarP, digit_t* scalarQ, const ec_basis_t* PQ2, const ec_point_t* R, const ec_curve_t* curve)
{ // Optimized implementation based on Montgomery formulas using Jacobian coordinates
    int i;
    digit_t w0 = 0, z0 = 0, x0 = 0, y0 = 0, x1 = 0, y1 = 0, w1 = 0, z1 = 0, e, e1, f, f1, f2, f2div2, f22, w, z;
    digit_t fp2[NWORDS_ORDER] = {0}, xx[NWORDS_ORDER] = {0}, yy[NWORDS_ORDER] = {0}, ww[NWORDS_ORDER] = {0}, zz[NWORDS_ORDER] = {0};
    digit_t f22_p[NWORDS_ORDER] = {0}, w_p[NWORDS_ORDER] = {0}, z_p[NWORDS_ORDER] = {0}, x0_p[NWORDS_ORDER] = {0}, y0_p[NWORDS_ORDER] = {0}, w0_p[NWORDS_ORDER] = {0}, z0_p[NWORDS_ORDER] = {0};
    jac_point_t P, Q, RR, TT, R2r0, R2r, R2r1;
    jac_point_t Pe2[POWER_OF_2], Qe2[POWER_OF_2], PQe2[POWER_OF_2];
    ec_point_t Rnorm;    
    ec_curve_t curvenorm;
    ec_basis_t PQ2norm;

    f = POWER_OF_2;
    memset(scalarP, 0, NWORDS_ORDER*RADIX/8);
    memset(scalarQ, 0, NWORDS_ORDER*RADIX/8);

    // Normalize R,PQ2,curve
    fp2_t D;
    fp2_mul(&D, &PQ2->P.z, &PQ2->Q.z);
    fp2_mul(&D, &D, &PQ2->PmQ.z);
    fp2_mul(&D, &D, &R->z);
    fp2_mul(&D, &D, &curve->C);
    fp2_inv(&D);
    fp_mont_setone(Rnorm.z.re);
    fp_set(Rnorm.z.im, 0);
    fp2_copy(&PQ2norm.P.z, &Rnorm.z);
    fp2_copy(&PQ2norm.Q.z, &Rnorm.z);
    fp2_copy(&PQ2norm.PmQ.z, &Rnorm.z);
    fp2_copy(&curvenorm.C, &Rnorm.z);
    fp2_mul(&Rnorm.x, &R->x, &D);
    fp2_mul(&Rnorm.x, &Rnorm.x, &PQ2->P.z);
    fp2_mul(&Rnorm.x, &Rnorm.x, &PQ2->Q.z);
    fp2_mul(&Rnorm.x, &Rnorm.x, &PQ2->PmQ.z);
    fp2_mul(&Rnorm.x, &Rnorm.x, &curve->C);
    fp2_mul(&PQ2norm.P.x, &PQ2->P.x, &D);
    fp2_mul(&PQ2norm.P.x, &PQ2norm.P.x, &R->z);
    fp2_mul(&PQ2norm.P.x, &PQ2norm.P.x, &PQ2->Q.z);
    fp2_mul(&PQ2norm.P.x, &PQ2norm.P.x, &PQ2->PmQ.z);
    fp2_mul(&PQ2norm.P.x, &PQ2norm.P.x, &curve->C);
    fp2_mul(&PQ2norm.Q.x, &PQ2->Q.x, &D);
    fp2_mul(&PQ2norm.Q.x, &PQ2norm.Q.x, &R->z);
    fp2_mul(&PQ2norm.Q.x, &PQ2norm.Q.x, &PQ2->P.z);
    fp2_mul(&PQ2norm.Q.x, &PQ2norm.Q.x, &PQ2->PmQ.z);
    fp2_mul(&PQ2norm.Q.x, &PQ2norm.Q.x, &curve->C);
    fp2_mul(&PQ2norm.PmQ.x, &PQ2->PmQ.x, &D);
    fp2_mul(&PQ2norm.PmQ.x, &PQ2norm.PmQ.x, &R->z);
    fp2_mul(&PQ2norm.PmQ.x, &PQ2norm.PmQ.x, &PQ2->P.z);
    fp2_mul(&PQ2norm.PmQ.x, &PQ2norm.PmQ.x, &PQ2->Q.z);
    fp2_mul(&PQ2norm.PmQ.x, &PQ2norm.PmQ.x, &curve->C);
    fp2_mul(&curvenorm.A, &curve->A, &D);
    fp2_mul(&curvenorm.A, &curvenorm.A, &R->z);
    fp2_mul(&curvenorm.A, &curvenorm.A, &PQ2->P.z);
    fp2_mul(&curvenorm.A, &curvenorm.A, &PQ2->Q.z);
    fp2_mul(&curvenorm.A, &curvenorm.A, &PQ2->PmQ.z);

    recover_y(&P.y, &PQ2norm.P.x, &curvenorm);
    fp2_copy(&P.x, &PQ2norm.P.x);
    fp2_copy(&P.z, &PQ2norm.P.z);
    recover_y(&Q.y, &PQ2norm.Q.x, &curvenorm);   // TODO: THIS SECOND SQRT CAN BE ELIMINATED
    fp2_copy(&Q.x, &PQ2norm.Q.x);
    fp2_copy(&Q.z, &PQ2norm.Q.z);
    recover_y(&RR.y, &Rnorm.x, &curvenorm);
    fp2_copy(&RR.x, &Rnorm.x);
    fp2_copy(&RR.z, &Rnorm.z);
    
    jac_neg(&TT, &Q);
    ADD(&TT, &P, &TT, &curvenorm);
    if (!is_jac_xz_equal(&TT, &PQ2norm.PmQ))
        jac_neg(&Q, &Q);    
    
    // Computing torsion-2^f points, multiples of P, Q and P+Q 
    copy_jac_point(&Pe2[POWER_OF_2-1], &P);
    copy_jac_point(&Qe2[POWER_OF_2-1], &Q);
    ADD(&PQe2[POWER_OF_2-1], &P, &Q, &curvenorm);       // P+Q

    for (i = 0; i < (POWER_OF_2-1); i++) {
        DBL(&Pe2[POWER_OF_2-i-2], &Pe2[POWER_OF_2-i-1], &curvenorm);
        DBL(&Qe2[POWER_OF_2-i-2], &Qe2[POWER_OF_2-i-1], &curvenorm);
        DBL(&PQe2[POWER_OF_2-i-2], &PQe2[POWER_OF_2-i-1], &curvenorm);
    }

    e = f;
    mp_shiftr(&e, 1, 1);
    f1 = f-e;
    e1 = f1;
    mp_shiftr(&e1, 1, 1);
    f2 = f1-e1;

    copy_jac_point(&TT, &RR);
    for (i = 0; i < (f-f2); i++) {
        DBL(&TT, &TT, &curvenorm);
    }
    // w0, z0 <- dlog2(2^(f-f2)*R, f2, f2 div 2)
    f2div2 = f2;
    mp_shiftr(&f2div2, 1, 1);
    ec_dlog_2_step(&w0, &z0, &TT, (int)f2, (int)f2div2, &Pe2[0], &Qe2[0], &PQe2[0], &curvenorm);

    // R2r0 <- 2^e*R2 - (w0*Pe2[f1-1] + z0*Qe2[f1-1])
    copy_jac_point(&TT, &RR);
    for (i = 0; i < e; i++) {
        DBL(&TT, &TT, &curvenorm);
    }
    DBLMUL(&R2r0, &Pe2[f1-1], w0, &Qe2[f1-1], z0, &curvenorm);
    jac_neg(&R2r0, &R2r0);
    ADD(&R2r0, &TT, &R2r0, &curvenorm);
    ec_dlog_2_step(&x0, &y0, &R2r0, (int)e1, (int)f2div2, &Pe2[0], &Qe2[0], &PQe2[0], &curvenorm);

    // w <- w0 + 2^f2 * x0, z <- z0 + 2^f2 * y0
    //f22 = 1;
    //mp_shiftl(&f22, (unsigned int)f2, 1);
    //w = w0 + f22 * x0;
    //z = z0 + f22 * y0;
    f22_p[0] = 1;
    mp_shiftl(&f22_p[0], (unsigned int)f2, NWORDS_ORDER);
    x0_p[0] = x0;
    y0_p[0] = y0;
    w0_p[0] = w0;
    z0_p[0] = z0;
    mp_mul2(&x0_p[0], &f22_p[0], &x0_p[0]);
    mp_add(&w_p[0], &w0_p[0], &x0_p[0], NWORDS_ORDER);
    mp_mul2(&y0_p[0], &f22_p[0], &y0_p[0]);
    mp_add(&z_p[0], &z0_p[0], &y0_p[0], NWORDS_ORDER);

    // R2r <- R2 - (w*Pe2[f-1] + z*Qe2[f-1]), R2r has order 2^(f-e)
    //DBLMUL(&R2r, &Pe2[f-1], w, &Qe2[f-1], z, &curvenorm);
    DBLMUL2(&R2r, &Pe2[f-1], &w_p[0], &Qe2[f-1], &z_p[0], &curvenorm);
    jac_neg(&R2r, &R2r);
    ADD(&R2r, &RR, &R2r, &curvenorm);
    copy_jac_point(&TT, &R2r);
    for (i = 0; i < e1; i++) {
        DBL(&TT, &TT, &curvenorm);
    }
    ec_dlog_2_step(&w1, &z1, &TT, (int)f2, (int)f2div2, &Pe2[0], &Qe2[0], &PQe2[0], &curvenorm);

    // R2r1 <- R2r - (w1*Pe2[f1-1] + z1*Qe2[f1-1])
    DBLMUL(&R2r1, &Pe2[f1-1], w1, &Qe2[f1-1], z1, &curvenorm);
    jac_neg(&R2r1, &R2r1);
    ADD(&R2r1, &R2r, &R2r1, &curvenorm);
    ec_dlog_2_step(&x1, &y1, &R2r1, (int)e1, (int)f2div2, &Pe2[0], &Qe2[0], &PQe2[0], &curvenorm);

    xx[0] = x1;
    yy[0] = y1;
    ww[0] = w1;
    zz[0] = z1;
    mp_shiftl(&xx[0], (unsigned int)f2, NWORDS_ORDER);
    mp_add(&xx[0], &xx[0], &ww[0], NWORDS_ORDER);
    mp_shiftl(&yy[0], (unsigned int)f2, NWORDS_ORDER);
    mp_add(&yy[0], &yy[0], &zz[0], NWORDS_ORDER);
    //ww[0] = w;
    //zz[0] = z;
    ww[0] = w_p[0];
    zz[0] = z_p[0];
    ww[1] = w_p[1];
    zz[1] = z_p[1];
    if (e < 64) {
        mp_shiftl(&xx[0], (unsigned int)e, NWORDS_ORDER);
        mp_shiftl(&yy[0], (unsigned int)e, NWORDS_ORDER);
    } else {
        x0_p[0] = 0;
        y0_p[0] = 0;
        x0_p[1] = 0x100;
        y0_p[1] = 0x100;
        mp_mul2(&xx[0], &xx[0], &x0_p[0]);
        mp_mul2(&yy[0], &yy[0], &y0_p[0]);
    }
    mp_add(scalarP, &xx[0], &ww[0], NWORDS_ORDER);
    mp_add(scalarQ, &yy[0], &zz[0], NWORDS_ORDER);

    fp_copy(fp2, TWOpFm1);   // 2^(f-1)

    // If scalarP > 2^(f-1) or (scalarQ > 2^(f-1) and (scalarP = 0 or scalarP = 2^(f-1))) then output -scalarP mod 2^f, -scalarQ mod 2^f
    if (mp_compare(scalarP, fp2, NWORDS_ORDER) == 1 ||
        (mp_compare(scalarQ, fp2, NWORDS_ORDER) == 1 && (mp_is_zero(scalarP, NWORDS_ORDER) == 1 || mp_compare(scalarP, fp2, NWORDS_ORDER) == 0))) {
        mp_shiftl(fp2, 1, NWORDS_ORDER);       // Get 2^f
        if (mp_is_zero(scalarP, NWORDS_ORDER) != 1)
            mp_sub(scalarP, fp2, scalarP, NWORDS_ORDER);
        if (mp_is_zero(scalarQ, NWORDS_ORDER) != 1)
            mp_sub(scalarQ, fp2, scalarQ, NWORDS_ORDER);
    }
}

// point doublings using Modified Jacobian coordinates
void DBL_MJ(jac_point_t *Q, fp2_t *T, fp2_t *lam, const jac_point_t *P)
{
    fp2_t t1, t2, t3, t4, t5;
    fp2_sqr(&t1, &P->x);      // t1 = xx^2;
    fp2_sqr(&t2, &P->y);
    fp2_add(&t2, &t2, &t2);   // t2 = 2*yy^2;
    fp2_sqr(&t3, &t2);        // t3 = t2^2;
    fp2_add(&t4, &t3, &t3);   // t4 = 2 * t3;
    fp2_add(&t5, &t2, &P->x);
    fp2_sqr(&t5, &t5);
    fp2_sub(&t5, &t5, &t1);
    fp2_sub(&t5, &t5, &t3);   // t5 = (xx+t2)^2-t1-t3;
    fp2_add(lam, &t1, &t1);
    fp2_add(lam, lam, &t1);
    fp2_add(lam, lam, T);     // tmp = 3*t1+tt;
    fp2_sqr(&Q->x, lam);
    fp2_sub(&Q->x,&Q->x, &t5);
    fp2_sub(&Q->x, &Q->x, &t5);// x3 = tmp^2-2*t5;
    fp2_sub(&Q->y, &t5, &Q->x);
    fp2_mul(&Q->y, &Q->y, lam);
    fp2_sub(&Q->y, &Q->y, &t4);// y3 = tmp*(t5-x3)-t4;
    fp2_mul(&Q->z, &P->z, &P->y);
    fp2_add(&Q->z, &Q->z, &Q->z);// z3 = 2*yy*zz;
    fp2_mul(&t3, T, &t4);
    fp2_add(T, &t3, &t3);     // t3 = 2*t4*tt;
}

// Miller loop
void line_fun(fp2_t *f, const jac_point_t *P, const jac_point_t *Q, fp2_t *lam1, fp2_t *lam2)
{
    fp2_t t1, t2, t3, t4, t5, t6, t7;
    fp2_sqr(&t1, &P->z);
    fp2_mul(&t2, &t1, &P->z);
    fp2_mul(&t3, &t2, lam2);
    fp2_mul(&t4, &t2, &P->y);
    fp2_sqr(f, f);
    fp2_mul(&t5, &t1, &Q->x);
    fp2_sub(&t5, &t5, &P->x);
    fp2_mul(&t6, &t5, lam1);
    fp2_mul(&t7, &t2, &Q->y);
    fp2_add(&t7, &t7, &P->y);
    fp2_sub(&t6, &t6, &t7);
    fp2_mul(f, f, &t6);
    fp2_sqr(f, f);
    fp2_mul(f, f, &P->y);
    fp2_add(f, f, f);
    fp2_mul(&t5, &t5, &t3);
    fp2_mul(&t6, &t4, &t7);
    fp2_add(&t6, &t6, &t6);
    fp2_add(&t5, &t5, &t6);
    fp_neg(t5.im, t5.im);
    fp2_mul(f, f, &t5);
}

// Miller loop (shared the first argument)
void line_fun_share(fp2_t *f1, fp2_t *f2, const jac_point_t *P, const jac_point_t *Q1, const jac_point_t *Q2, fp2_t *lam1, fp2_t *lam2)
{
    fp2_t t1, t2, t3, t4, t5, t6, t7;
    fp2_sqr(&t1, &P->z);
    fp2_mul(&t2, &t1, &P->z);
    fp2_mul(&t3, &t2, lam2);
    fp2_mul(&t4, &t2, &P->y);
    fp2_sqr(f1, f1);
    fp2_mul(&t5, &t1, &Q1->x);
    fp2_sub(&t5, &t5, &P->x);
    fp2_mul(&t6, &t5, lam1);
    fp2_mul(&t7, &t2, &Q1->y);
    fp2_add(&t7, &t7, &P->y);
    fp2_sub(&t6, &t6, &t7);
    fp2_mul(f1, f1, &t6);
    fp2_sqr(f1, f1);
    fp2_mul(f1, f1, &P->y);
    fp2_add(f1, f1, f1);
    fp2_mul(&t5, &t5, &t3);
    fp2_mul(&t6, &t4, &t7);
    fp2_add(&t6, &t6, &t6);
    fp2_add(&t5, &t5, &t6);
    fp_neg(t5.im, t5.im);
    fp2_mul(f1, f1, &t5);
    fp2_sqr(f2, f2);
    fp2_mul(&t5, &t1, &Q2->x);
    fp2_sub(&t5, &t5, &P->x);
    fp2_mul(&t6, &t5, lam1);
    fp2_mul(&t7, &t2, &Q2->y);
    fp2_add(&t7, &t7, &P->y);
    fp2_sub(&t6, &t6, &t7);
    fp2_mul(f2, f2, &t6);
    fp2_sqr(f2, f2);
    fp2_mul(f2, f2, &P->y);
    fp2_add(f2, f2, f2);
    fp2_mul(&t5, &t5, &t3);
    fp2_mul(&t6, &t4, &t7);
    fp2_add(&t6, &t6, &t6);
    fp2_add(&t5, &t5, &t6);
    fp_neg(t5.im, t5.im);
    fp2_mul(f2, f2, &t5);
}

static inline void Lucassequences(fp2_t *f, fp_t v0, fp_t v1)
{
    fp_copy(v0, fp_2);
    fp_t tmp1;
    fp_add(tmp1, f->re, f->re);
    fp_copy(v1, tmp1);
    fp_t tmp2;
    fp_copy(tmp2, fp_2);
    for (uint64_t i = p_even_cofactor_minus_1_highest; i > 0; i >>= 1)
    {
        if (i & p_even_cofactor_minus_1[p_even_cofactor_minus_1_limb])
        {
            fp_mul(v0, v0, v1);
            fp_sub(v0, v0, tmp1);
            fp_sqr(v1, v1);
            fp_sub(v1, v1, tmp2);
        }
        else
        {
            fp_mul(v1, v1, v0);
            fp_sub(v1, v1, tmp1);
            fp_sqr(v0, v0);
            fp_sub(v0, v0, tmp2);
        }
    }
    for (int j = p_even_cofactor_minus_1_limb - 1; j >= 0; j--)
    {
        for (uint64_t i = 0x8000000000000000; i > 0; i >>= 1)
        {
            if (i & p_even_cofactor_minus_1[j])
            {
                fp_mul(v0, v0, v1);
                fp_sub(v0, v0, tmp1);
                fp_sqr(v1, v1);
                fp_sub(v1, v1, tmp2);
            }
            else
            {
                fp_mul(v1, v1, v0);
                fp_sub(v1, v1, tmp1);
                fp_sqr(v0, v0);
                fp_sub(v0, v0, tmp2);
            }
        }
    }
}

// compute f1^{p+1}, f2^{p+1} and f3^{p+1} with the help of Lucas sequences
void ELS(fp2_t *f1, fp2_t *f2, fp2_t *f3)
{
    fp_t tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
    Lucassequences(f1, tmp1, tmp2);
    fp_mul(tmp1, tmp1, fp_inv2);
    fp_mul(tmp2, tmp2, fp_inv2);
    fp_t re1;
    fp_copy(re1, tmp2);
    fp_mul(tmp2, tmp2, f1->re);
    fp_sub(tmp2, tmp2, tmp1);
    fp_sqr(f1->re, f1->re);
    fp_sub(f1->re, f1->re, fp_1);
    Lucassequences(f2, tmp1, tmp3);
    fp_mul(tmp1, tmp1, fp_inv2);
    fp_mul(tmp3, tmp3, fp_inv2);
    fp_t re2;
    fp_copy(re2, tmp3);
    fp_mul(tmp3, tmp3, f2->re);
    fp_sub(tmp3, tmp3, tmp1);
    fp_sqr(f2->re, f2->re);
    fp_sub(f2->re, f2->re, fp_1);
    Lucassequences(f3, tmp1, tmp4);
    fp_mul(tmp1, tmp1, fp_inv2);
    fp_mul(tmp4, tmp4, fp_inv2);
    fp_t re3;
    fp_copy(re3, tmp4);
    fp_mul(tmp4, tmp4, f3->re);
    fp_sub(tmp4, tmp4, tmp1);
    fp_sqr(f3->re, f3->re);
    fp_sub(f3->re, f3->re, fp_1);
    fp_mul(tmp1, f2->re, f3->re);
    fp_mul(tmp5, f1->re, tmp1);
    fp_copy(tmp6, f1->re);
    fp_inv(tmp5);
    fp_mul(f1->re, tmp5, tmp1);
    fp_mul(tmp5, tmp5, tmp6);
    fp_copy(tmp6, f2->re);
    fp_mul(f2->re, tmp5, f3->re);
    fp_mul(f3->re, tmp5, tmp6);
    fp_mul(tmp2, tmp2, f1->re);
    fp_mul(f1->im, f1->im, tmp2);
    fp_mul(tmp3, tmp3, f2->re);
    fp_mul(f2->im, f2->im, tmp3);
    fp_mul(tmp4, tmp4, f3->re);
    fp_mul(f3->im, f3->im, tmp4);
    fp_copy(f1->re, re1);
    fp_copy(f2->re, re2);
    fp_copy(f3->re, re3);
}

// Pairing computation
void pairing_share(fp2_t *v1, fp2_t *v2, fp2_t *v3, const jac_point_t *P, const jac_point_t *Q1, const jac_point_t *Q2, const fp2_t *WA, long f)
{
    jac_point_t P_J, PP_J;
    fp2_t lam1, lam2;
    fp2_t t1, t2, t3, t4, t5;
    P_J.x = P->x;
    P_J.y = P->y;
    P_J.z = P->z;
    fp2_t T;
    fp2_copy(&T, WA);
    for(int i = 0; i < f/2; i++)
    {
        DBL_MJ(&PP_J, &T, &lam1, &P_J);
        DBL_MJ(&P_J, &T, &lam2, &PP_J);
        line_fun_share(v1, v2, &PP_J, Q1, Q2, &lam1, &lam2);
    }
    fp2_sqr(&t1, &P_J.z);
    fp2_mul(&t2, &t1, &Q1->x);
    fp2_sub(&t2, &t2, &P_J.x);
    fp2_mul(&t3, &t1, &Q2->x);
    fp2_sub(&t3, &t3, &P_J.x);
    fp2_copy(&P_J.x, &Q1->x);
    fp2_copy(&P_J.y, &Q1->y);
    fp2_copy(&P_J.z, &Q1->z);
    fp2_copy(&T, WA);
    for(int i = 0; i < f/2; i++)
    {
        DBL_MJ(&PP_J, &T, &lam1, &P_J);
        DBL_MJ(&P_J, &T, &lam2, &PP_J);
        line_fun(v3, &PP_J, Q2, &lam1, &lam2);
    }
    fp2_sqr(&t4, &P_J.z);
    fp2_mul(&t5, &t4, &Q2->x);
    fp2_sub(&t5, &t5, &P_J.x);
    fp_neg(t1.im, t1.im);
    fp_neg(t4.im, t4.im);
    fp2_mul(&t2, &t2, &t1);
    fp2_mul(&t5, &t5, &t4);
    fp2_sqr(v1, v1);
    fp2_mul(v1, v1, &t2);
    fp2_mul(&t3, &t3, &t1);
    fp2_sqr(v2, v2);
    fp2_mul(v2, v2, &t3);
    fp2_sqr(v3, v3);
    fp2_mul(v3, v3, &t5);
    fp2_copy(&t2, v1);
    fp2_copy(&t2, v1);
    fp2_copy(&t3, v1);
    fp_neg(t2.im, t2.im);
    fp2_mul(&t1, v1, v2);
    fp2_mul(&t4, &t1, v3);
    fp2_inv(&t4);
    fp2_mul(&t5, &t1, &t4);
    fp2_mul(&t1, &t4, v3);
    fp2_copy(&t4, v3);
    fp_neg(t4.im, t4.im);
    fp2_mul(v3, &t5, &t4);
    fp2_mul(v1, &t1, v2);
    fp2_mul(v1, v1, &t2);
    fp2_copy(&t2, v2);
    fp_neg(t2.im, t2.im);
    fp2_mul(v2, &t1, &t3);
    fp2_mul(v2, v2, &t2);
    ELS(v1,v2,v3);
}


// DLP in the small group
void small_DLP(int *r, int *sign, fp2_t h, const fp_t *R, int row_length)
{
    *sign = 0;
    fp_t one;
    fp_mont_setone(one);
    if(fp_is_equal(h.re, one))
    {
        *r = 0;
        *sign = 1;
    }
    else
    {
        for(int i = 0; i < row_length; i++)
        {
            if(fp_is_equal(h.re, R[2*i]))
            {
                *r = i + 1;
                if(fp_is_equal(h.im, R[2*i+1]))
                {
                    *sign = 1;
                }
                else
                {
                    *sign = -1;
                }
            }
        }
        assert(*sign!=0);
    }
}

// DLP in the large group with Pohlig-Hellman algorithm
void PH_DLP(int_two_torsion *answer, int *ans, fp2_t h, const int Str[], int w, int column_length, int row_length, long f)
{
    fp2_t tmp1, tmp2, tmp3;
    fp_copy(tmp1.re, Lookupbasere);
    fp_copy(tmp1.im, Lookupbaseim);
    int time = column_length;
    int modw = f%w;
    int pw = 1;
    for(int ii = 0; ii <w; ii++)
    {
        pw = pw * 2;
    }
    int a = 0, sign = 0;
    struct
    {
        fp2_t v;
        int i1;
        int i2;
    }stack[column_length];
    int top = 0;
    fp2_copy(&stack[top].v, &h);
    for(int i = 0; i < modw; i++)
    {
        fp2_sqr(&stack[top].v, &stack[top].v);
    }
    stack[top].i1 = 0;
    stack[top].i2 = 0;
    int i = 0, j = 0, k = 0;
    while(k != time - 1)
    {
        while (j + k != time - 1)
        {
            fp2_t ht;
            fp2_copy(&ht, &stack[top].v);
            j = j + Str[i];
            for (int ii = 0; ii < Str[i] * w; ii++) {
                fp2_sqr(&ht, &ht);
            }
            top++;
            fp2_copy(&stack[top].v, &ht);
            stack[top].i1 = j + k;
            stack[top].i2 = Str[i];
            i++;
        }
        fp2_t hc;
        fp2_copy(&hc, &stack[top].v);
        j = j - stack[top].i2;
        top--;
        a = 0, sign = 0;
        small_DLP(&a, &sign, hc, LookuptableLR, row_length);
        if(a != 0)
        {
            for(int ii = 0; ii < top + 1; ii++)
            {
                int coord = stack[ii].i1*pw+2*a-2;
                fp_copy(tmp2.re, Lookuptable[coord]);
                fp_copy(tmp2.im, Lookuptable[coord+1]);
                if(sign == 1)
                {
                    fp_neg(tmp2.im, tmp2.im);
                }
                fp2_mul(&stack[ii].v, &stack[ii].v, &tmp2);
            }
        }
        for(int ii = 0; ii < top + 1; ii++)
        {
            stack[ii].i1++;
        }
        if(sign == 1)
        {
            ans[k] = a;
        }
        else
        {
            ans[k] = -a;
        }
        k++;
    }
    fp2_copy(&tmp3, &stack[top].v);
    small_DLP(&a, &sign, tmp3, LookuptableLR, row_length);
    if(sign == 1)
    {
        ans[k] = a;
    }
    else
    {
        ans[k] = -a;
    }
    if(modw != 0)
    {
        fp_mont_setone(tmp3.re);
        fp_set(tmp3.im, 0);
        int ind = ans[0];
        if(ind < 0)
        {
            ind = -ind;
        }
        if(ind == 0)
        {
            fp_mont_setone(tmp2.re);
            fp_set(tmp2.im, 0);
        }
        else
        {
            if(ind % 2 == 0)
            {
                fp_mont_setone(tmp2.re);
                fp_set(tmp2.im, 0);
            }
            else
            {
                fp_copy(tmp2.re, Lookupbasere);
                fp_copy(tmp2.im, Lookupbaseim);
            }
            while(ind != 1)
            {
                fp2_sqr(&tmp1, &tmp1);
                ind = ind/2;
                if(ind%2 == 1)
                {
                    fp2_mul(&tmp2, &tmp1, &tmp2);
                }
            }
            if(ans[0] > 0)
            {
                fp_neg(tmp2.im, tmp2.im);
            }
        }
        for(int ii=1; ii<time; ii++)
        {
            if(ans[ii]!=0)
            {
                if(ans[ii] > 0)
                {
                    sign = 1;
                }
                else
                {
                    sign = -1;
                }
                int coord = (ii-1)*pw+2*abs(ans[ii])-2;
                fp_copy(tmp1.re, Lookuptable[coord]);
                fp_copy(tmp1.im, Lookuptable[coord+1]);
                if(sign == 1)
                {
                    fp_neg(tmp1.im, tmp1.im);
                }
                fp2_mul(&tmp3, &tmp1, &tmp3);
            }
        }
        for(int ii=0; ii<w-modw; ii++)
        {
            fp2_sqr(&tmp3, &tmp3);
        }
        fp2_mul(&tmp3, &tmp3, &tmp2);
        fp2_mul(&tmp3, &tmp3, &h);
        small_DLP(&a, &sign, tmp3, LookuptableLR, row_length);

        int div = 1;
        for(int ii = 0; ii < w - modw; ii++)
        {
            div = div * 2;
        }
        if(sign == 1)
        {
            ans[k+1] = a/div;
        }
        else
        {
            ans[k+1] = -a/div;
        }
    }
    __int128_t x = 0;
    __int128_t pow2_f = 1;
    __int128_t pow2_64 = 1;
    for(int ii = 0; ii < f; ii++)
    {
        pow2_f = pow2_f * 2;
        if(ii == 63)
        {
            pow2_64 = pow2_f;
        }
    }
    __int128_t ip = 1;
    int modwlabel = 0;
    if(modw != 0)
    {
        modwlabel = 1;
    }
    for(int ii = 0; ii < column_length + modwlabel; ii++)
    {
        x = x + ans[ii] * ip;
        ip <<= w;
    }
    if(x < 0)
    {
        x = x + pow2_f;
    }
    *answer = x;
}

void ec_dlog_2_improved(digit_t* scalarP, digit_t* scalarQ, const ec_basis_t* PQ2, const ec_point_t* R, const ec_curve_t* curve)
{ // Optimized implementation based on Montgomery formulas using Jacobian coordinates
    int i;
    digit_t w0 = 0, z0 = 0, x0 = 0, y0 = 0, x1 = 0, y1 = 0, w1 = 0, z1 = 0, e, e1, f, f1, f2, f2div2, f22, w, z;
    digit_t fp2[NWORDS_ORDER] = {0}, xx[NWORDS_ORDER] = {0}, yy[NWORDS_ORDER] = {0}, ww[NWORDS_ORDER] = {0}, zz[NWORDS_ORDER] = {0};
    digit_t f22_p[NWORDS_ORDER] = {0}, w_p[NWORDS_ORDER] = {0}, z_p[NWORDS_ORDER] = {0}, x0_p[NWORDS_ORDER] = {0}, y0_p[NWORDS_ORDER] = {0}, w0_p[NWORDS_ORDER] = {0}, z0_p[NWORDS_ORDER] = {0};
    jac_point_t P, Q, RR, TT, R2r0, R2r, R2r1;
    jac_point_t Pe2[POWER_OF_2], Qe2[POWER_OF_2], PQe2[POWER_OF_2];
    ec_point_t Rnorm;    
    ec_curve_t curvenorm;
    ec_basis_t PQ2norm;

    f = POWER_OF_2;
    memset(scalarP, 0, NWORDS_ORDER*RADIX/8);
    memset(scalarQ, 0, NWORDS_ORDER*RADIX/8);

    // Normalize R,PQ2,curve
    fp2_t D;
    fp2_mul(&D, &PQ2->P.z, &PQ2->Q.z);
    fp2_mul(&D, &D, &PQ2->PmQ.z);
    fp2_mul(&D, &D, &R->z);
    fp2_mul(&D, &D, &curve->C);
    fp2_inv(&D);
    fp_mont_setone(Rnorm.z.re);
    fp_set(Rnorm.z.im, 0);
    fp2_copy(&PQ2norm.P.z, &Rnorm.z);
    fp2_copy(&PQ2norm.Q.z, &Rnorm.z);
    fp2_copy(&PQ2norm.PmQ.z, &Rnorm.z);
    fp2_copy(&curvenorm.C, &Rnorm.z);
    fp2_mul(&Rnorm.x, &R->x, &D);
    fp2_mul(&Rnorm.x, &Rnorm.x, &PQ2->P.z);
    fp2_mul(&Rnorm.x, &Rnorm.x, &PQ2->Q.z);
    fp2_mul(&Rnorm.x, &Rnorm.x, &PQ2->PmQ.z);
    fp2_mul(&Rnorm.x, &Rnorm.x, &curve->C);
    fp2_mul(&PQ2norm.P.x, &PQ2->P.x, &D);
    fp2_mul(&PQ2norm.P.x, &PQ2norm.P.x, &R->z);
    fp2_mul(&PQ2norm.P.x, &PQ2norm.P.x, &PQ2->Q.z);
    fp2_mul(&PQ2norm.P.x, &PQ2norm.P.x, &PQ2->PmQ.z);
    fp2_mul(&PQ2norm.P.x, &PQ2norm.P.x, &curve->C);
    fp2_mul(&PQ2norm.Q.x, &PQ2->Q.x, &D);
    fp2_mul(&PQ2norm.Q.x, &PQ2norm.Q.x, &R->z);
    fp2_mul(&PQ2norm.Q.x, &PQ2norm.Q.x, &PQ2->P.z);
    fp2_mul(&PQ2norm.Q.x, &PQ2norm.Q.x, &PQ2->PmQ.z);
    fp2_mul(&PQ2norm.Q.x, &PQ2norm.Q.x, &curve->C);
    fp2_mul(&PQ2norm.PmQ.x, &PQ2->PmQ.x, &D);
    fp2_mul(&PQ2norm.PmQ.x, &PQ2norm.PmQ.x, &R->z);
    fp2_mul(&PQ2norm.PmQ.x, &PQ2norm.PmQ.x, &PQ2->P.z);
    fp2_mul(&PQ2norm.PmQ.x, &PQ2norm.PmQ.x, &PQ2->Q.z);
    fp2_mul(&PQ2norm.PmQ.x, &PQ2norm.PmQ.x, &curve->C);
    fp2_mul(&curvenorm.A, &curve->A, &D);
    fp2_mul(&curvenorm.A, &curvenorm.A, &R->z);
    fp2_mul(&curvenorm.A, &curvenorm.A, &PQ2->P.z);
    fp2_mul(&curvenorm.A, &curvenorm.A, &PQ2->Q.z);
    fp2_mul(&curvenorm.A, &curvenorm.A, &PQ2->PmQ.z);

    recover_y(&P.y, &PQ2norm.P.x, &curvenorm);
    fp2_copy(&P.x, &PQ2norm.P.x);
    fp2_copy(&P.z, &PQ2norm.P.z);
    recover_y(&Q.y, &PQ2norm.Q.x, &curvenorm);   // TODO: THIS SECOND SQRT CAN BE ELIMINATED
    fp2_copy(&Q.x, &PQ2norm.Q.x);
    fp2_copy(&Q.z, &PQ2norm.Q.z);
    recover_y(&RR.y, &Rnorm.x, &curvenorm);
    fp2_copy(&RR.x, &Rnorm.x);
    fp2_copy(&RR.z, &Rnorm.z);
    
    jac_neg(&TT, &Q);
    ADD(&TT, &P, &TT, &curvenorm);
    if (!is_jac_xz_equal(&TT, &PQ2norm.PmQ))
        jac_neg(&Q, &Q);    
    
    fp2_t v1, v2, v3;
    fp_mont_setone(v1.re);
    fp_set(v1.im,0);
    fp_mont_setone(v2.re);
    fp_set(v2.im,0);
    fp_mont_setone(v3.re);
    fp_set(v3.im,0);
    
    fp2_t tmpA;
    jac_point_t tmpP, tmpQ, tmpRR;
    fp2_copy(&tmpP.x, &P.x);
    fp2_copy(&tmpP.y, &P.y);
    fp2_copy(&tmpP.z, &P.z);
    fp2_copy(&tmpQ.x, &Q.x);
    fp2_copy(&tmpQ.y, &Q.y);
    fp2_copy(&tmpQ.z, &Q.z);
    fp2_copy(&tmpRR.x, &RR.x);
    fp2_copy(&tmpRR.y, &RR.y);
    fp2_copy(&tmpRR.z, &RR.z);
    fp2_copy(&tmpA, &curvenorm.A);
    fp_mul(tmpA.re, tmpA.re, fp_inv3);
    fp_mul(tmpA.im, tmpA.im, fp_inv3);
    fp2_add(&tmpP.x, &tmpP.x, &tmpA);
    fp2_add(&tmpQ.x, &tmpQ.x, &tmpA);
    fp2_add(&tmpRR.x, &tmpRR.x, &tmpA);
    fp2_mul(&tmpA, &tmpA, &curvenorm.A);
    fp_sub(tmpA.re, fp_1, tmpA.re);
    fp_neg(tmpA.im, tmpA.im);

    // Pairing computation
    fp2_t tmpm;
    fp2_t tmpv;
    pairing_share(&v1,&v2,&v3,&tmpP,&tmpQ,&tmpRR,&tmpA,f);

    // Discrete logarithms over finite fields
    int win = 5;
    int column_length = f/win;
    int row_length = 1;
    for(int ii = 0; ii < win - 1; ii++)
    {
        row_length = row_length * 2;
    }
    int ans[column_length];

    int_two_torsion is1, is2, is3;
    PH_DLP(&is2, ans, v2, strategy, win, column_length, row_length,f);
    int2t_neg(&is2, is2);
    PH_DLP(&is3, ans, v3, strategy, win, column_length, row_length,f);
    PH_DLP(&is1, ans, v1, strategy, win, column_length, row_length,f);
    int2t_inv(&is1, is1);
    int2t_mul(&is2, is1, is2);
    int2t_mul(&is3, is1, is3);

    int2t_toarray(scalarQ, is2);
    int2t_toarray(scalarP, is3);

    fp_copy(fp2, TWOpFm1);   // 2^(f-1)
    // If scalarP > 2^(f-1) or (scalarQ > 2^(f-1) and (scalarP = 0 or scalarP = 2^(f-1))) then output -scalarP mod 2^f, -scalarQ mod 2^f
    if (mp_compare(scalarP, fp2, NWORDS_ORDER) == 1 ||
        (mp_compare(scalarQ, fp2, NWORDS_ORDER) == 1 && (mp_is_zero(scalarP, NWORDS_ORDER) == 1 || mp_compare(scalarP, fp2, NWORDS_ORDER) == 0))) {
        mp_shiftl(fp2, 1, NWORDS_ORDER);       // Get 2^f
        if (mp_is_zero(scalarP, NWORDS_ORDER) != 1)
            mp_sub(scalarP, fp2, scalarP, NWORDS_ORDER);
        if (mp_is_zero(scalarQ, NWORDS_ORDER) != 1)
            mp_sub(scalarQ, fp2, scalarQ, NWORDS_ORDER);
    }
}

static void ec_dlog_3_step(digit_t* x, digit_t* y, const jac_point_t* R, const int f, const int B, const jac_point_t* Pe3, const jac_point_t* Qe3, const jac_point_t* PQe3, const ec_curve_t* curve)
{ // Based on Montgomery formulas using Jacobian coordinates
    int i, j;
    digit_t value = 1;
    jac_point_t P, Q, PQp, PQ2p, P2Qp, P2Q2p, P2p, Q2p, TT, SS, Te[POWER_OF_3/2], Re[POWER_OF_3/2];    // Storage could be reduced to e points
    jac_point_t PQep0, PQ2ep0, P2Qep0, P2Q2ep0, P2ep0, Q2ep0, PQep1, PQ2ep1, P2Qep1, P2Q2ep1, P2ep1, Q2ep1;

    *x = 0;       
    *y = 0;

    copy_jac_point(&P, &Pe3[f-1]);
    copy_jac_point(&Q, &Qe3[f-1]);
    copy_jac_point(&Re[f-1], R);

    for (i = 0; i < (f-1); i++) {
        TPL(&Re[f-i-2], &Re[f-i-1], curve);
    }

    ADD(&PQp, &P, &Q, curve);                // P+Q
    ADD(&PQ2p, &PQp, &Q, curve);             // P+2Q
    ADD(&P2Qp, &PQp, &P, curve);             // P2+Q
    DBL(&P2Q2p, &PQp, curve);                // 2P+2Q
    DBL(&P2p, &P, curve);                    // 2P
    DBL(&Q2p, &Q, curve);                    // 2Q
    ADD(&PQep0, &Pe3[0], &Qe3[0], curve);    // Pe3[0] + Qe3[0]
    ADD(&PQ2ep0, &PQep0, &Qe3[0], curve);    // Pe3[0] + 2Qe3[0]
    ADD(&P2Qep0, &PQep0, &Pe3[0], curve);    // 2Pe3[0] + Qe3[0]
    DBL(&P2Q2ep0, &PQep0, curve);            // 2Pe3[0] + 2Qe3[0]
    DBL(&P2ep0, &Pe3[0], curve);             // 2Pe3[0]
    DBL(&Q2ep0, &Qe3[0], curve);             // 2Qe3[0]
    ADD(&PQep1, &Pe3[1], &Qe3[1], curve);    // Pe3[1] + Qe3[1]
    ADD(&PQ2ep1, &PQep1, &Qe3[1], curve);    // Pe3[1] + 2Qe3[1]
    ADD(&P2Qep1, &PQep1, &Pe3[1], curve);    // 2Pe3[1] + Qe3[1]
    DBL(&P2Q2ep1, &PQep1, curve);            // 2Pe3[1] + 2Qe3[1]
    DBL(&P2ep1, &Pe3[1], curve);             // 2Pe3[1]
    DBL(&Q2ep1, &Qe3[1], curve);             // 2Qe3[1]

    // Unrolling the first two iterations
    if (is_jac_equal(&Pe3[0], &Re[0])) {
        *x += 1;
        copy_jac_point(&Te[0], &P);
        for (j = 3; j <= B; j++) {
            copy_jac_point(&Te[j-1], &Pe3[j-1]);
        }
        jac_neg(&TT, &Pe3[1]);
        ADD(&Re[0], &Re[1], &TT, curve);
    } else if (is_jac_equal(&PQep0, &Re[0])) {
        *x += 1;
        *y += 1;
        copy_jac_point(&Te[0], &PQp);
        for (j = 3; j <= B; j++) {
            ADD(&Te[j-1], &Pe3[j-1], &Qe3[j-1], curve);
        }
        jac_neg(&TT, &PQep1);
        ADD(&Re[0], &Re[1], &TT, curve);
    } else if (is_jac_equal(&PQ2ep0, &Re[0])) {
        *x += 1;
        *y += 2;
        copy_jac_point(&Te[0], &PQ2p);
        for (j = 3; j <= B; j++) {
            DBL(&TT, &Qe3[j-1], curve);
            ADD(&Te[j-1], &Pe3[j-1], &TT, curve);
        }
        jac_neg(&TT, &PQ2ep1);
        ADD(&Re[0], &Re[1], &TT, curve);
    } else if (is_jac_equal(&P2ep0, &Re[0])) {
        *x += 2;
        copy_jac_point(&Te[0], &P2p);
        for (j = 3; j <= B; j++) {
            DBL(&Te[j-1], &Pe3[j-1], curve);
        }
        jac_neg(&TT, &P2ep1);
        ADD(&Re[0], &Re[1], &TT, curve);
    } else if (is_jac_equal(&P2Qep0, &Re[0])) {
        *x += 2;
        *y += 1;
        copy_jac_point(&Te[0], &P2Qp);
        for (j = 3; j <= B; j++) {
            DBL(&TT, &Pe3[j-1], curve);
            ADD(&Te[j-1], &TT, &Qe3[j-1], curve);
        }
        jac_neg(&TT, &P2Qep1);
        ADD(&Re[0], &Re[1], &TT, curve);
    } else if (is_jac_equal(&P2Q2ep0, &Re[0])) {
        *x += 2;
        *y += 2;
        copy_jac_point(&Te[0], &P2Q2p);
        for (j = 3; j <= B; j++) {
            DBL(&TT, &Pe3[j-1], curve);
            DBL(&SS, &Qe3[j-1], curve);
            ADD(&Te[j-1], &TT, &SS, curve);
        }
        jac_neg(&TT, &P2Q2ep1);
        ADD(&Re[0], &Re[1], &TT, curve);
    } else if (is_jac_equal(&Qe3[0], &Re[0])) {
        *y += 1;
        copy_jac_point(&Te[0], &Q);
        for (j = 3; j <= B; j++) {
            copy_jac_point(&Te[j-1], &Qe3[j-1]);
        }
        jac_neg(&TT, &Qe3[1]);
        ADD(&Re[0], &Re[1], &TT, curve);
    } else if (is_jac_equal(&Q2ep0, &Re[0])) {
        *y += 2;
        copy_jac_point(&Te[0], &Q2p);
        for (j = 3; j <= B; j++) {
            DBL(&Te[j-1], &Qe3[j-1], curve);
        }
        jac_neg(&TT, &Q2ep1);
        ADD(&Re[0], &Re[1], &TT, curve);
    } else {
        jac_init(&Te[0]);
        for (j = 3; j <= B; j++) {
            jac_init(&Te[j-1]);
        }
        copy_jac_point(&Re[0], &Re[1]);
    }

    // Unrolling iterations 3-B
    for (i = 3; i <= B; i++) {
        value *= 3;             
        jac_neg(&TT, &Te[i-1]);
        ADD(&TT, &Re[i-1], &TT, curve);       // TT <- Re[i-1]-T[i-1]
        if (is_jac_equal(&Pe3[0], &Re[0])) {
            *x += value;
            ADD(&Te[0], &Te[0], &Pe3[f-i+1], curve);
            for (j = i+1; j <= B; j++) {
                ADD(&Te[j-1], &Te[j-1], &Pe3[j-i+1], curve);
            }
            jac_neg(&SS, &Pe3[1]);
            ADD(&Re[0], &TT, &SS, curve);
        } else if (is_jac_equal(&PQep0, &Re[0])) {
            *x += value;
            *y += value;
            ADD(&Te[0], &Te[0], &Pe3[f-i+1], curve);
            ADD(&Te[0], &Te[0], &Qe3[f-i+1], curve);
            for (j = i+1; j <= B; j++) {
                ADD(&Te[j-1], &Te[j-1], &Pe3[j-i+1], curve);
                ADD(&Te[j-1], &Te[j-1], &Qe3[j-i+1], curve);
            }
            jac_neg(&SS, &PQep1);
            ADD(&Re[0], &TT, &SS, curve);
        } else if (is_jac_equal(&PQ2ep0, &Re[0])) {
            *x += value;
            *y += value;
            *y += value;
            ADD(&Te[0], &Te[0], &Pe3[f-i+1], curve);
            ADD(&Te[0], &Te[0], &Qe3[f-i+1], curve);
            ADD(&Te[0], &Te[0], &Qe3[f-i+1], curve);
            for (j = i+1; j <= B; j++) {
                ADD(&Te[j-1], &Te[j-1], &Pe3[j-i+1], curve);
                ADD(&Te[j-1], &Te[j-1], &Qe3[j-i+1], curve);
                ADD(&Te[j-1], &Te[j-1], &Qe3[j-i+1], curve);
            }
            jac_neg(&SS, &PQ2ep1);
            ADD(&Re[0], &TT, &SS, curve);
        } else if (is_jac_equal(&P2ep0, &Re[0])) {
            *x += value;
            *x += value;
            ADD(&Te[0], &Te[0], &Pe3[f-i+1], curve);
            ADD(&Te[0], &Te[0], &Pe3[f-i+1], curve);
            for (j = i+1; j <= B; j++) {
                ADD(&Te[j-1], &Te[j-1], &Pe3[j-i+1], curve);
                ADD(&Te[j-1], &Te[j-1], &Pe3[j-i+1], curve);
            }
            jac_neg(&SS, &P2ep1);
            ADD(&Re[0], &TT, &SS, curve);
        } else if (is_jac_equal(&P2Qep0, &Re[0])) {
            *x += value;
            *x += value;
            *y += value;
            ADD(&Te[0], &Te[0], &Pe3[f-i+1], curve);
            ADD(&Te[0], &Te[0], &Pe3[f-i+1], curve);
            ADD(&Te[0], &Te[0], &Qe3[f-i+1], curve);
            for (j = i+1; j <= B; j++) {
                ADD(&Te[j-1], &Te[j-1], &Pe3[j-i+1], curve);
                ADD(&Te[j-1], &Te[j-1], &Pe3[j-i+1], curve);
                ADD(&Te[j-1], &Te[j-1], &Qe3[j-i+1], curve);
            }
            jac_neg(&SS, &P2Qep1);
            ADD(&Re[0], &TT, &SS, curve);
        } else if (is_jac_equal(&P2Q2ep0, &Re[0])) {
            *x += value;
            *x += value;
            *y += value;
            *y += value;
            ADD(&Te[0], &Te[0], &Pe3[f-i+1], curve);
            ADD(&Te[0], &Te[0], &Pe3[f-i+1], curve);
            ADD(&Te[0], &Te[0], &Qe3[f-i+1], curve);
            ADD(&Te[0], &Te[0], &Qe3[f-i+1], curve);
            for (j = i+1; j <= B; j++) {
                ADD(&Te[j-1], &Te[j-1], &Pe3[j-i+1], curve);
                ADD(&Te[j-1], &Te[j-1], &Pe3[j-i+1], curve);
                ADD(&Te[j-1], &Te[j-1], &Qe3[j-i+1], curve);
                ADD(&Te[j-1], &Te[j-1], &Qe3[j-i+1], curve);
            }
            jac_neg(&SS, &P2Q2ep1);
            ADD(&Re[0], &TT, &SS, curve);
        } else if (is_jac_equal(&Qe3[0], &Re[0])) {
            *y += value;
            ADD(&Te[0], &Te[0], &Qe3[f-i+1], curve);
            for (j = i+1; j <= B; j++) {
                ADD(&Te[j-1], &Te[j-1], &Qe3[j-i+1], curve);
            }
            jac_neg(&SS, &Qe3[1]);
            ADD(&Re[0], &TT, &SS, curve);
        } else if (is_jac_equal(&Q2ep0, &Re[0])) {
            *y += value;
            *y += value;
            ADD(&Te[0], &Te[0], &Qe3[f-i+1], curve);
            ADD(&Te[0], &Te[0], &Qe3[f-i+1], curve);
            for (j = i+1; j <= B; j++) {
                ADD(&Te[j-1], &Te[j-1], &Qe3[j-i+1], curve);
                ADD(&Te[j-1], &Te[j-1], &Qe3[j-i+1], curve);
            }
            jac_neg(&SS, &Q2ep1);
            ADD(&Re[0], &TT, &SS, curve);
        } else {
            copy_jac_point(&Re[0], &TT);
        }
    }

    // Main Loop
    for (i = B; i < f; i++) {
        value *= 3;             
        if (is_jac_equal(&Pe3[0], &Re[0])) {
            *x += value;
            copy_jac_point(&TT, &P);
            for (j = 0; j < (i-1); j++) {
                TPL(&TT, &TT, curve);
            }
            ADD(&Te[0], &Te[0], &TT, curve);
        } else if (is_jac_equal(&PQep0, &Re[0])) {
            *x += value;
            *y += value;
            copy_jac_point(&TT, &PQp);
            for (j = 0; j < (i-1); j++) {
                TPL(&TT, &TT, curve);
            }
            ADD(&Te[0], &Te[0], &TT, curve);
        } else if (is_jac_equal(&PQ2ep0, &Re[0])) {
            *x += value;
            *y += value;
            *y += value;
            copy_jac_point(&TT, &PQ2p);
            for (j = 0; j < (i-1); j++) {
                TPL(&TT, &TT, curve);
            }
            ADD(&Te[0], &Te[0], &TT, curve);
        } else if (is_jac_equal(&P2ep0, &Re[0])) {
            *x += value;
            *x += value;
            copy_jac_point(&TT, &P2p);
            for (j = 0; j < (i-1); j++) {
                TPL(&TT, &TT, curve);
            }
            ADD(&Te[0], &Te[0], &TT, curve);
        } else if (is_jac_equal(&P2Qep0, &Re[0])) {
            *x += value;
            *x += value;
            *y += value;
            copy_jac_point(&TT, &P2Qp);
            for (j = 0; j < (i-1); j++) {
                TPL(&TT, &TT, curve);
            }
            ADD(&Te[0], &Te[0], &TT, curve);
        } else if (is_jac_equal(&P2Q2ep0, &Re[0])) {
            *x += value;
            *x += value;
            *y += value;
            *y += value;
            copy_jac_point(&TT, &P2Q2p);
            for (j = 0; j < (i-1); j++) {
                TPL(&TT, &TT, curve);
            }
            ADD(&Te[0], &Te[0], &TT, curve);
        } else if (is_jac_equal(&Qe3[0], &Re[0])) {
            *y += value;
            copy_jac_point(&TT, &Q);
            for (j = 0; j < (i-1); j++) {
                TPL(&TT, &TT, curve);
            }
            ADD(&Te[0], &Te[0], &TT, curve);
        } else if (is_jac_equal(&Q2ep0, &Re[0])) {
            *y += value;
            *y += value;
            copy_jac_point(&TT, &Q2p);
            for (j = 0; j < (i-1); j++) {
                TPL(&TT, &TT, curve);
            }
            ADD(&Te[0], &Te[0], &TT, curve);
        }
        jac_neg(&TT, &Te[0]);
        ADD(&Re[0], R, &TT, curve);
        for (j = 0; j < (f-i-1); j++) {
            TPL(&Re[0], &Re[0], curve);
        }
    }
    value *= 3;           
    if (is_jac_equal(&Pe3[0], &Re[0])) {
        *x += value;
    } else if (is_jac_equal(&PQep0, &Re[0])) {
        *x += value;
        *y += value;
    } else if (is_jac_equal(&PQ2ep0, &Re[0])) {
        *x += value;
        *y += value;
        *y += value;
    } else if (is_jac_equal(&P2ep0, &Re[0])) {
        *x += value;
        *x += value;
    } else if (is_jac_equal(&P2Qep0, &Re[0])) {
        *x += value;
        *x += value;
        *y += value;
    } else if (is_jac_equal(&P2Q2ep0, &Re[0])) {
        *x += value;
        *x += value;
        *y += value;
        *y += value;
    } else if (is_jac_equal(&Qe3[0], &Re[0])) {
        *y += value;
    } else if (is_jac_equal(&Q2ep0, &Re[0])) {
        *y += value;
        *y += value;
    }
}

void ec_dlog_3(digit_t* scalarP, digit_t* scalarQ, const ec_basis_t* PQ3, const ec_point_t* R, const ec_curve_t* curve)
{ // Optimized implementation based on Montgomery formulas using Jacobian coordinates
    int i;
    digit_t w0 = 0, z0 = 0, x0 = 0, y0 = 0, x1 = 0, y1 = 0, w1 = 0, z1 = 0, e, f, f1, f1div2;
    digit_t fp2[NWORDS_ORDER] = {0}, xx[NWORDS_ORDER] = {0}, yy[NWORDS_ORDER] = {0}, ww[NWORDS_ORDER] = {0}, zz[NWORDS_ORDER] = {0};
    jac_point_t P, Q, RR, TT, R2r0;
    jac_point_t Pe3[POWER_OF_3], Qe3[POWER_OF_3], PQe3[POWER_OF_3];
    ec_point_t Rnorm;
    ec_curve_t curvenorm;
    ec_basis_t PQ3norm;

    f = POWER_OF_3;
    memset(scalarP, 0, NWORDS_ORDER*RADIX/8);
    memset(scalarQ, 0, NWORDS_ORDER*RADIX/8);

    // Normalize R,PQ3,curve
    fp2_t D;
    fp2_mul(&D, &PQ3->P.z, &PQ3->Q.z);
    fp2_mul(&D, &D, &PQ3->PmQ.z);
    fp2_mul(&D, &D, &R->z);
    fp2_mul(&D, &D, &curve->C);
    fp2_inv(&D);
    fp_mont_setone(Rnorm.z.re);
    fp_set(Rnorm.z.im, 0);
    fp2_copy(&PQ3norm.P.z, &Rnorm.z);
    fp2_copy(&PQ3norm.Q.z, &Rnorm.z);
    fp2_copy(&PQ3norm.PmQ.z, &Rnorm.z);
    fp2_copy(&curvenorm.C, &Rnorm.z);
    fp2_mul(&Rnorm.x, &R->x, &D);
    fp2_mul(&Rnorm.x, &Rnorm.x, &PQ3->P.z);
    fp2_mul(&Rnorm.x, &Rnorm.x, &PQ3->Q.z);
    fp2_mul(&Rnorm.x, &Rnorm.x, &PQ3->PmQ.z);
    fp2_mul(&Rnorm.x, &Rnorm.x, &curve->C);
    fp2_mul(&PQ3norm.P.x, &PQ3->P.x, &D);
    fp2_mul(&PQ3norm.P.x, &PQ3norm.P.x, &R->z);
    fp2_mul(&PQ3norm.P.x, &PQ3norm.P.x, &PQ3->Q.z);
    fp2_mul(&PQ3norm.P.x, &PQ3norm.P.x, &PQ3->PmQ.z);
    fp2_mul(&PQ3norm.P.x, &PQ3norm.P.x, &curve->C);
    fp2_mul(&PQ3norm.Q.x, &PQ3->Q.x, &D);
    fp2_mul(&PQ3norm.Q.x, &PQ3norm.Q.x, &R->z);
    fp2_mul(&PQ3norm.Q.x, &PQ3norm.Q.x, &PQ3->P.z);
    fp2_mul(&PQ3norm.Q.x, &PQ3norm.Q.x, &PQ3->PmQ.z);
    fp2_mul(&PQ3norm.Q.x, &PQ3norm.Q.x, &curve->C);
    fp2_mul(&PQ3norm.PmQ.x, &PQ3->PmQ.x, &D);
    fp2_mul(&PQ3norm.PmQ.x, &PQ3norm.PmQ.x, &R->z);
    fp2_mul(&PQ3norm.PmQ.x, &PQ3norm.PmQ.x, &PQ3->P.z);
    fp2_mul(&PQ3norm.PmQ.x, &PQ3norm.PmQ.x, &PQ3->Q.z);
    fp2_mul(&PQ3norm.PmQ.x, &PQ3norm.PmQ.x, &curve->C);
    fp2_mul(&curvenorm.A, &curve->A, &D);
    fp2_mul(&curvenorm.A, &curvenorm.A, &R->z);
    fp2_mul(&curvenorm.A, &curvenorm.A, &PQ3->P.z);
    fp2_mul(&curvenorm.A, &curvenorm.A, &PQ3->Q.z);
    fp2_mul(&curvenorm.A, &curvenorm.A, &PQ3->PmQ.z);

    recover_y(&P.y, &PQ3norm.P.x, &curvenorm);
    fp2_copy(&P.x, &PQ3norm.P.x);
    fp2_copy(&P.z, &PQ3norm.P.z);
    recover_y(&Q.y, &PQ3norm.Q.x, &curvenorm);      // TODO: THIS SECOND SQRT CAN BE ELIMINATED
    fp2_copy(&Q.x, &PQ3norm.Q.x);
    fp2_copy(&Q.z, &PQ3norm.Q.z);
    recover_y(&RR.y, &Rnorm.x, &curvenorm);
    fp2_copy(&RR.x, &Rnorm.x);
    fp2_copy(&RR.z, &Rnorm.z);

    jac_neg(&TT, &Q);
    ADD(&TT, &P, &TT, &curvenorm);
    if (!is_jac_xz_equal(&TT, &PQ3norm.PmQ))
        jac_neg(&Q, &Q);

    // Computing torsion-2^f points, multiples of P, Q and P+Q 
    copy_jac_point(&Pe3[POWER_OF_3-1], &P);
    copy_jac_point(&Qe3[POWER_OF_3-1], &Q);
    ADD(&PQe3[POWER_OF_3-1], &P, &Q, &curvenorm);       // P+Q

    for (i = 0; i < (POWER_OF_3-1); i++) {
        TPL(&Pe3[POWER_OF_3-i-2], &Pe3[POWER_OF_3-i-1], &curvenorm);
        TPL(&Qe3[POWER_OF_3-i-2], &Qe3[POWER_OF_3-i-1], &curvenorm);
        TPL(&PQe3[POWER_OF_3-i-2], &PQe3[POWER_OF_3-i-1], &curvenorm);
    }

    e = f >> 1;
    f1 = f-e;
    copy_jac_point(&TT, &RR);
    for (i = 0; i < (f-f1); i++) {
        TPL(&TT, &TT, &curvenorm);
    }
    // w0, z0 <- dlog3(3^(f-f1)*R, f1, f1 div 2)
    f1div2 = f1 >> 1;
    ec_dlog_3_step(&w0, &z0, &TT, (int)f1, (int)f1div2, &Pe3[0], &Qe3[0], &PQe3[0], &curvenorm);

    // R2r0 <- R2 - (w0*Pe3[f-1] + z0*Qe3[f-1])
    DBLMUL(&R2r0, &Pe3[f-1], w0, &Qe3[f-1], z0, &curvenorm);
    jac_neg(&R2r0, &R2r0);
    ADD(&R2r0, &RR, &R2r0, &curvenorm);
    ec_dlog_3_step(&x0, &y0, &R2r0, (int)e, (int)f1div2, &Pe3[0], &Qe3[0], &PQe3[0], &curvenorm);

    xx[0] = x0;
    yy[0] = y0;
    ww[0] = w0;
    zz[0] = z0;
    mp_mul2(&xx[0], THREEpE, &xx[0]);       
    mp_add(scalarP, &xx[0], &ww[0], NWORDS_ORDER);
    mp_mul2(&yy[0], THREEpE, &yy[0]);       
    mp_add(scalarQ, &yy[0], &zz[0], NWORDS_ORDER);

    // If scalarP > Floor(3^f/2) or (scalarQ > Floor(3^f/2) and scalarP = 0) then output -scalarP mod 3^f, -scalarQ mod 3^f
    if (mp_compare(scalarP, THREEpFdiv2, NWORDS_ORDER) == 1 ||
       (mp_compare(scalarQ, THREEpFdiv2, NWORDS_ORDER) == 1 && (mp_is_zero(scalarP, NWORDS_ORDER) == 1))) {
        if (mp_is_zero(scalarP, NWORDS_ORDER) != 1)
            mp_sub(scalarP, THREEpF, scalarP, NWORDS_ORDER);
        if (mp_is_zero(scalarQ, NWORDS_ORDER) != 1)
            mp_sub(scalarQ, THREEpF, scalarQ, NWORDS_ORDER);
    }
}


// WRAPPERS

void ec_dbl(ec_point_t* res, const ec_curve_t* curve, const ec_point_t* P){
    xDBL(res, P, (ec_point_t const*)curve);
}
void ec_mul(ec_point_t* res, const ec_curve_t* curve, const digit_t* scalar, const ec_point_t* P){
    xMUL(res, P, scalar, curve);
}
void ec_biscalar_mul(ec_point_t* res, const ec_curve_t* curve, const digit_t* scalarP, const digit_t* scalarQ, const ec_basis_t* PQ){
    xDBLMUL(res, &PQ->P, scalarP, &PQ->Q, scalarQ, &PQ->PmQ, curve);
}
