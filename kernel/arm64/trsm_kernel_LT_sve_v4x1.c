/***************************************************************************
(c) RIKEN 2024, 2024. All rights reserved. trsm_kernel_LT_sve_v4x1.c 0.3.26
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

   1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the
      distribution.
   3. Neither the name of the OpenBLAS project nor the names of
      its contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/
/*********************************************************************/
/* Copyright 2009, 2010 The University of Texas at Austin.           */
/* All rights reserved.                                              */
/*                                                                   */
/* Redistribution and use in source and binary forms, with or        */
/* without modification, are permitted provided that the following   */
/* conditions are met:                                               */
/*                                                                   */
/*   1. Redistributions of source code must retain the above         */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer.                                                  */
/*                                                                   */
/*   2. Redistributions in binary form must reproduce the above      */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer in the documentation and/or other materials       */
/*      provided with the distribution.                              */
/*                                                                   */
/*    THIS  SOFTWARE IS PROVIDED  BY THE  UNIVERSITY OF  TEXAS AT    */
/*    AUSTIN  ``AS IS''  AND ANY  EXPRESS OR  IMPLIED WARRANTIES,    */
/*    INCLUDING, BUT  NOT LIMITED  TO, THE IMPLIED  WARRANTIES OF    */
/*    MERCHANTABILITY  AND FITNESS FOR  A PARTICULAR  PURPOSE ARE    */
/*    DISCLAIMED.  IN  NO EVENT SHALL THE UNIVERSITY  OF TEXAS AT    */
/*    AUSTIN OR CONTRIBUTORS BE  LIABLE FOR ANY DIRECT, INDIRECT,    */
/*    INCIDENTAL,  SPECIAL, EXEMPLARY,  OR  CONSEQUENTIAL DAMAGES    */
/*    (INCLUDING, BUT  NOT LIMITED TO,  PROCUREMENT OF SUBSTITUTE    */
/*    GOODS  OR  SERVICES; LOSS  OF  USE,  DATA,  OR PROFITS;  OR    */
/*    BUSINESS INTERRUPTION) HOWEVER CAUSED  AND ON ANY THEORY OF    */
/*    LIABILITY, WHETHER  IN CONTRACT, STRICT  LIABILITY, OR TORT    */
/*    (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY WAY OUT    */
/*    OF  THE  USE OF  THIS  SOFTWARE,  EVEN  IF ADVISED  OF  THE    */
/*    POSSIBILITY OF SUCH DAMAGE.                                    */
/*                                                                   */
/* The views and conclusions contained in the software and           */
/* documentation are those of the authors and should not be          */
/* interpreted as representing official policies, either expressed   */
/* or implied, of The University of Texas at Austin.                 */
/*********************************************************************/

#include "common.h"
#include "arm_sve.h"

static FLOAT dm1 = -1.;

#ifdef CONJ
#define GEMM_KERNEL   GEMM_KERNEL_L
#else
#define GEMM_KERNEL   GEMM_KERNEL_N
#endif

#ifndef COMPLEX

static inline void solve(BLASLONG m, BLASLONG n, FLOAT *a, FLOAT *b, BLASLONG ldb, FLOAT *c, BLASLONG ldc) {

  FLOAT aa, bb;
  FLOAT *bw;

  int i, j, k;

  for (i = 0; i < m; i++) {

    bw=b+i;
    aa = *(a + i);

    for (j = 0; j < n; j ++) {
      bb = *(c + i + j * ldc);
      bb *= aa;
      *bw            = bb;
      *(c + i + j * ldc) = bb;
      bw+=ldb;

      for (k = i + 1; k < m; k ++){
	*(c + k + j * ldc) -= bb * *(a + k);
      }

    }
    a += m;
  }
}

#else
#error not supported
#endif


int CNAME(BLASLONG m, BLASLONG n, BLASLONG k, FLOAT dummy1,
#ifdef COMPLEX
#error not supported
#endif
	   FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc, BLASLONG offset){

  FLOAT *aa, *cc;
  BLASLONG  kk;
  BLASLONG i, j;

  FLOAT *bp,*bw;
  BLASLONG ik,in,iw;

  bw=malloc(k*GEMM_UNROLL_N*SIZE*COMPSIZE);
  if(bw==NULL){
    fprintf(stderr, "OpenBLAS: malloc failed in %s\n", __func__);
    exit(1);
  }

#ifdef DOUBLE
#error not supported
#else
  int sve_size = svcntw();
#endif
  int unroll_size = sve_size*4;

#if 0
  fprintf(stderr, "TRSM KERNEL LT : m = %3ld  n = %3ld  k = %3ld offset = %3ld\n",
	  m, n, k, offset);
#endif

  j = n;

  while (j>=GEMM_UNROLL_N) {
    kk = offset;
    aa = a;
    cc = c;

    i = unroll_size;

    while (i <= m) {
      if (kk > 0) {

        iw=0;
        for(in=0;in<GEMM_UNROLL_N;in++) {
          bp=b + in * k * COMPSIZE;
          for(ik=0;ik<kk*COMPSIZE;ik++) bw[iw++]=bp[ik];
        }

        GEMM_KERNEL(unroll_size, GEMM_UNROLL_N, kk, dm1,
#ifdef COMPLEX
#error not supported
#endif
            aa,
            bw,
            cc,
            ldc);
      }

      solve(unroll_size, GEMM_UNROLL_N,
          aa + kk * unroll_size * COMPSIZE,
          b  + kk * COMPSIZE, k,  cc, ldc);

      aa += unroll_size * k * COMPSIZE;
      cc += unroll_size     * COMPSIZE;
      kk += unroll_size;
      i += unroll_size;
    }

    i = m % unroll_size;
    if (i) {
      if (kk > 0) {

        iw=0;
        for(in=0;in<GEMM_UNROLL_N;in++) {
          bp=b + in * k * COMPSIZE;
          for(ik=0;ik<kk*COMPSIZE;ik++) bw[iw++]=bp[ik];
        }

        GEMM_KERNEL(i, GEMM_UNROLL_N, kk, dm1,
#ifdef COMPLEX
#error not supported
#endif
            aa, bw, cc, ldc);
      }

      solve(i, GEMM_UNROLL_N,
          aa + kk * i * COMPSIZE,
          b  + kk * COMPSIZE, k, cc, ldc);

      aa += i * k * COMPSIZE;
      cc += i     * COMPSIZE;
      kk += i;

    }
    b += GEMM_UNROLL_N * k   * COMPSIZE;
    c += GEMM_UNROLL_N * ldc * COMPSIZE;
    j-=GEMM_UNROLL_N;
  }

  while (j>0) {
    kk = offset;
    aa = a;
    cc = c;

    i = unroll_size;

    while (i <= m) {
      if (kk > 0) {
        GEMM_KERNEL(unroll_size, 1, kk, dm1,
#ifdef COMPLEX
#error not supported
#endif
            aa,
            b,
            cc,
            ldc);
      }

      solve(unroll_size, 1,
          aa + kk * unroll_size * COMPSIZE,
          b  + kk * COMPSIZE, k,  cc, ldc);

      aa += unroll_size * k * COMPSIZE;
      cc += unroll_size     * COMPSIZE;
      kk += unroll_size;
      i += unroll_size;
    }

    i = m % unroll_size;
    if (i) {
      if (kk > 0) {
        GEMM_KERNEL(i, 1, kk, dm1,
#ifdef COMPLEX
#error not supported
#endif
            aa, b, cc, ldc);
      }

      solve(i, 1,
          aa + kk * i * COMPSIZE,
          b  + kk * COMPSIZE, k, cc, ldc);

      aa += i * k * COMPSIZE;
      cc += i     * COMPSIZE;
      kk += i;

    }
    b += k   * COMPSIZE;
    c += ldc * COMPSIZE;
    j --;
  }

  free(bw);

  return 0;
}
