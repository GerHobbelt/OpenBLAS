/***************************************************************************
(c) RIKEN 2024, 2024. All rights reserved. trsm_kernel_RT_sve_v4x1.c 0.3.26
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
#error not supported
#else
#define GEMM_KERNEL   GEMM_KERNEL_N
#endif

#ifndef COMPLEX

static inline void solve(BLASLONG m, BLASLONG n, FLOAT *a, FLOAT *b, BLASLONG ldb,  FLOAT *c, BLASLONG ldc) {

  FLOAT aa,  bb;

  int i, j, k;

  a += (n - 1) * m;
  b += n - 1;

  for (i = n - 1; i >= 0; i--) {

    bb = *(b + i * ldb);

    for (j = 0; j < m; j ++) {
      aa = *(c + j + i * ldc);
      aa *= bb;
      *a   = aa;
      *(c + j + i * ldc) = aa;
      a ++;

      for (k = 0; k < i; k ++){
	*(c + j + k * ldc) -= aa * *(b + k * ldb );
      }

    }
    b --;
    a -= 2 * m;
  }

}

#else
#error not supported
#endif

int CNAME(BLASLONG m, BLASLONG n, BLASLONG k,  FLOAT dummy1,
#ifdef COMPLEX
#error not supported
#endif
    FLOAT *a, FLOAT *b, FLOAT *c, BLASLONG ldc, BLASLONG offset){

  BLASLONG i, j;
  FLOAT *aa, *cc;
  BLASLONG  kk;
#ifdef DOUBLE
#error not supported
#else
  int sve_size = svcntw();
#endif
  int unroll_size = sve_size * 4;

  kk = n - offset;
  c += n * ldc * COMPSIZE;
  b += n * k   * COMPSIZE;

  j = n;
  if (j > 0) {
    do {
      aa  = a;
      b -= k  * COMPSIZE;
      c -= ldc* COMPSIZE;
      cc  = c;

      i = unroll_size;
      if (i <= m) {

        do {
          if (k - kk > 0) {
            GEMM_KERNEL(unroll_size, 1, k - kk, dm1,
#ifdef COMPLEX
#error not supported
#endif
                aa + unroll_size * kk * COMPSIZE,
                b  + kk * COMPSIZE,
                cc,
                ldc);
          }

          solve(unroll_size, 1,
              aa + (kk - 1) * unroll_size * COMPSIZE,
              b  + (kk - 1) * COMPSIZE, k,
              cc, ldc);

          aa += unroll_size * k * COMPSIZE;
          cc += unroll_size     * COMPSIZE;
          i += unroll_size;
        } while (i <= m);
      }

      i = m % unroll_size;
      if (i) {
        if (k - kk > 0) {
          GEMM_KERNEL(i, 1, k - kk, dm1,
#ifdef COMPLEX
#error not supported
#endif
              aa + i * kk * COMPSIZE,
              b  + kk * COMPSIZE,
              cc, ldc);
        }

        solve(i, 1,
            aa + (kk - 1) * i * COMPSIZE,
            b  + (kk - 1) * COMPSIZE, k,
            cc, ldc);

        aa += i * k * COMPSIZE;
        cc += i     * COMPSIZE;

      }
      kk--;
      j--;
    } while (j > 0);
  }

  return 0;
}


