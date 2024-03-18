/***************************************************************************
(c) RIKEN 2024, 2024. All rights reserved. trsm_utcopy_sve_v4.c 0.3.26
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
/* Copyright 2023 The OpenBLAS Project                               */
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

#include <stdio.h>
#include "common.h"
#include "arm_sve.h"

#ifndef UNIT
#define INV(a) (ONE / (a))
#else
#define INV(a) (ONE)
#endif

int CNAME(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, BLASLONG offset, FLOAT *b){

  BLASLONG i, ii, jj;

  FLOAT *ao;

  jj = offset;
#ifdef DOUBLE
#error DOUBLE is not supported
#else
  uint32_t sve_size = svcntw(),sve_size2=sve_size*2,sve_size3=sve_size*3;
  int32_t N = n;
  int32_t js = 0;
  svbool_t pn1 = svwhilelt_b32((uint64_t)js, (uint64_t)n);
  svbool_t pn2 = svwhilelt_b32((uint64_t)js+sve_size, (uint64_t)n);
  svbool_t pn3 = svwhilelt_b32((uint64_t)js+sve_size2, (uint64_t)n);
  svbool_t pn4 = svwhilelt_b32((uint64_t)js+sve_size3, (uint64_t)n);
  int n_active = svcntp_b32(svptrue_b32(), pn1)+
                 svcntp_b32(svptrue_b32(), pn2)+
                 svcntp_b32(svptrue_b32(), pn3)+
                 svcntp_b32(svptrue_b32(), pn4);

#endif
  do {

    ao = a;

    i = 0;
    ii = 0;
    do {

      if (ii == jj) {
        for (int j = 0; j < n_active; j++) {
          for (int k = 0; k < j; k++) {
            *(b + j * n_active + k) = *(ao + j * lda + k);
          }
          *(b + j * n_active + j) = INV(*(ao + j * lda + j));
        }
        ao += lda * n_active;
        b += n_active * n_active;
        i += n_active;
        ii += n_active;
      } else {
        if (ii > jj) {
#ifdef DOUBLE
#error DOUBLE is not supported
#else
          svfloat32_t aj_vec1 = svld1(pn1, ao);
          svfloat32_t aj_vec2 = svld1_vnum(pn2, ao,1);
          svfloat32_t aj_vec3 = svld1_vnum(pn3, ao,2);
          svfloat32_t aj_vec4 = svld1_vnum(pn4, ao,3);
#endif
          svst1(pn1, b, aj_vec1);
          svst1_vnum(pn2, b,1, aj_vec2);
          svst1_vnum(pn3, b,2, aj_vec3);
          svst1_vnum(pn4, b,3, aj_vec4);
        }
        ao += lda;
        b += n_active;
        i ++;
        ii ++;
      } 
    } while (i < m);


    a += n_active;
    jj += n_active;

    js += n_active;
#ifdef DOUBLE
#error DOUBLE is not supported
#else
        pn1 = svwhilelt_b32((uint64_t)js, (uint64_t)N);
        pn2 = svwhilelt_b32((uint64_t)js+sve_size, (uint64_t)N);
        pn3 = svwhilelt_b32((uint64_t)js+sve_size2, (uint64_t)N);
        pn4 = svwhilelt_b32((uint64_t)js+sve_size3, (uint64_t)N);
        n_active = svcntp_b32(svptrue_b32(), pn1)+
                   svcntp_b32(svptrue_b32(), pn2)+
                   svcntp_b32(svptrue_b32(), pn3)+
                   svcntp_b32(svptrue_b32(), pn4);
  } while (n_active>0);
#endif

return 0;
}
