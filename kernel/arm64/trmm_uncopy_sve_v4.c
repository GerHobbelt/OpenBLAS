/***************************************************************************
(c) RIKEN 2024, 2024. All rights reserved. trmm_uncopy_sve_v4.c 0.3.26
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

#include <stdio.h>
#include "common.h"

#ifdef __ARM_FEATURE_SVE
#include <arm_sve.h>
#endif

int CNAME(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, BLASLONG posX, BLASLONG posY, FLOAT *b){

    BLASLONG i, js;
    BLASLONG X;

    js = 0;
    FLOAT *ao;
#ifdef DOUBLE
#error DOUBLE is not supported
#else
    uint32_t sve_size = svcntw(),sve_size2=sve_size*2,sve_size3=sve_size*3;
    svint32_t index1 = svindex_s32(0, lda);
    svint32_t index2 = svindex_s32(lda*sve_size, lda);
    svint32_t index3 = svindex_s32(lda*sve_size2, lda);
    svint32_t index4 = svindex_s32(lda*sve_size3, lda);
    svbool_t pn1 = svwhilelt_b32((uint64_t)js, (uint64_t)n);
    svbool_t pn2 = svwhilelt_b32((uint64_t)js+sve_size, (uint64_t)n);
    svbool_t pn3 = svwhilelt_b32((uint64_t)js+sve_size2, (uint64_t)n);
    svbool_t pn4 = svwhilelt_b32((uint64_t)js+sve_size3, (uint64_t)n);
    int n_active = svcntp_b32(svptrue_b32(), pn1)+
                   svcntp_b32(svptrue_b32(), pn2)+
                   svcntp_b32(svptrue_b32(), pn3)+
                   svcntp_b32(svptrue_b32(), pn4);
#endif
    do
    {
        X = posX;

        if (posX <= posY) {
            ao = a + posX + posY * lda;
        } else {
            ao = a + posY + posX * lda;
        }

        i = 0;
        do 
        {
            if (X < posY) {
#ifdef DOUBLE
#error DOUBLE is not supported
#else
                svfloat32_t aj_vec1 = svld1_gather_index(pn1, ao, index1);
                svfloat32_t aj_vec2 = svld1_gather_index(pn2, ao, index2);
                svfloat32_t aj_vec3 = svld1_gather_index(pn3, ao, index3);
                svfloat32_t aj_vec4 = svld1_gather_index(pn4, ao, index4);
#endif
                svst1(pn1, b, aj_vec1);
                svst1_vnum(pn2, b,1, aj_vec2);
                svst1_vnum(pn3, b,2, aj_vec3);
                svst1_vnum(pn4, b,3, aj_vec4);
                ao ++;
                b += n_active;
                X ++;
                i ++;
            } else 
                if (X > posY) {
                    ao += lda;
                    b += n_active;
                    X ++;
                    i ++;
                } else {
                    /* I did not find a way to unroll this while preserving vector-length-agnostic code. */
#ifdef UNIT
                    int temp = 0;
                    for (int j = 0; j < n_active; j++) {
                        for (int k = 0 ; k < j; k++) {
                            b[temp++] = ZERO;
                        }
                        b[temp++] = ONE;
                        for (int k = j+1; k < n_active; k++) {
                            b[temp++] = *(ao+k*lda+j);
                        }
                    }
#else 
                    int temp = 0;
                    for (int j = 0; j < n_active; j++) {
                        for (int k = 0 ; k < j; k++) {
                            b[temp++] = ZERO;
                        }
                        for (int k = j; k < n_active; k++) {
                            b[temp++] = *(ao+k*lda+j);
                        }
                    }
#endif
                    ao += n_active;
                    b += n_active*n_active;
                    X += n_active;
                    i += n_active;
                }
        } while (i < m);

        posY += n_active;
        js += n_active;
#ifdef DOUBLE
#error DOUBLE is not supported
#else
        pn1 = svwhilelt_b32((uint64_t)js, (uint64_t)n);
        pn2 = svwhilelt_b32((uint64_t)js+sve_size, (uint64_t)n);
        pn3 = svwhilelt_b32((uint64_t)js+sve_size2, (uint64_t)n);
        pn4 = svwhilelt_b32((uint64_t)js+sve_size3, (uint64_t)n);
        n_active = svcntp_b32(svptrue_b32(), pn1)+
                       svcntp_b32(svptrue_b32(), pn2)+
                       svcntp_b32(svptrue_b32(), pn3)+
                       svcntp_b32(svptrue_b32(), pn4);
    } while (n_active>0);
#endif

    return 0;
}
