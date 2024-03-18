/***************************************************************************
(c) RIKEN 2024, 2024. All rights reserved. symm_ucopy_sve_v4.c 0.3.26
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
#include <arm_sve.h>

int CNAME(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, BLASLONG posX, BLASLONG posY, FLOAT *b){

  BLASLONG i, offset;
#if defined(DOUBLE)

#error DOUBLE is not supported

#else
  uint32_t sve_size = svcntw(),sve_size2=sve_size*2,sve_size3=sve_size*3;
  uint32_t unroll_size = sve_size*4;
  svint32_t posY_vec = svdup_s32(posY);
  svint32_t posX_vec = svdup_s32(posX);
  svint32_t lda_vec = svdup_s32(lda);
  svint32_t one_vec = svdup_s32(1);

  int32_t N = n;
  int32_t j = 0;

  svbool_t pg1 = svwhilelt_b32((uint32_t)j, (uint32_t)N);
  svbool_t pg2 = svwhilelt_b32((uint32_t)j+sve_size, (uint32_t)N);
  svbool_t pg3 = svwhilelt_b32((uint32_t)j+sve_size2, (uint32_t)N);
  svbool_t pg4 = svwhilelt_b32((uint32_t)j+sve_size3, (uint32_t)N);

  int32_t active = svcntp_b32(svptrue_b32(), pg1)+
                   svcntp_b32(svptrue_b32(), pg2)+
                   svcntp_b32(svptrue_b32(), pg3)+
                   svcntp_b32(svptrue_b32(), pg4);
  svint32_t index_neg1 = svindex_s32(0, -1);
  svint32_t index_neg2 = svindex_s32(-sve_size, -1);
  svint32_t index_neg3 = svindex_s32(-sve_size2, -1);
  svint32_t index_neg4 = svindex_s32(-sve_size3, -1);
  svint32_t index1 = svindex_s32(0, 1);
  svint32_t index2 = svindex_s32(sve_size, 1);
  svint32_t index3 = svindex_s32(sve_size2, 1);
  svint32_t index4 = svindex_s32(sve_size3, 1);
  do {
    offset = posX - posY;
    svint32_t vec_off = svdup_s32(offset);
    svbool_t cmp1 = svcmpgt(pg1, vec_off, index_neg1);
    svbool_t cmp2 = svcmpgt(pg2, vec_off, index_neg2);
    svbool_t cmp3 = svcmpgt(pg3, vec_off, index_neg3);
    svbool_t cmp4 = svcmpgt(pg4, vec_off, index_neg4);

    svint32_t temp_1 = svadd_z(pg1, posX_vec, index1);
    svint32_t temp_2 = svadd_z(pg2, posX_vec, index2);
    svint32_t temp_3 = svadd_z(pg3, posX_vec, index3);
    svint32_t temp_4 = svadd_z(pg4, posX_vec, index4);
    svint32_t temp1_1 = svmla_z(pg1, temp_1, posY_vec, lda_vec);
    svint32_t temp1_2 = svmla_z(pg2, temp_2, posY_vec, lda_vec);
    svint32_t temp1_3 = svmla_z(pg3, temp_3, posY_vec, lda_vec);
    svint32_t temp1_4 = svmla_z(pg4, temp_4, posY_vec, lda_vec);
    svint32_t temp2_1 = svmla_z(pg1, posY_vec, temp_1, lda);
    svint32_t temp2_2 = svmla_z(pg2, posY_vec, temp_2, lda);
    svint32_t temp2_3 = svmla_z(pg3, posY_vec, temp_3, lda);
    svint32_t temp2_4 = svmla_z(pg4, posY_vec, temp_4, lda);
    svint32_t gat_ind1 = svsel(cmp1, temp2_1, temp1_1);
    svint32_t gat_ind2 = svsel(cmp2, temp2_2, temp1_2);
    svint32_t gat_ind3 = svsel(cmp3, temp2_3, temp1_3);
    svint32_t gat_ind4 = svsel(cmp4, temp2_4, temp1_4);

    i = m;
    while (i>0) {
        svfloat32_t data_vec1 = svld1_gather_index(pg1, a, gat_ind1);
        svfloat32_t data_vec2 = svld1_gather_index(pg2, a, gat_ind2);
        svfloat32_t data_vec3 = svld1_gather_index(pg3, a, gat_ind3);
        svfloat32_t data_vec4 = svld1_gather_index(pg4, a, gat_ind4);

        gat_ind1 = svadd_m(cmp1, gat_ind1, one_vec);
        gat_ind2 = svadd_m(cmp2, gat_ind2, one_vec);
        gat_ind3 = svadd_m(cmp3, gat_ind3, one_vec);
        gat_ind4 = svadd_m(cmp4, gat_ind4, one_vec);

        gat_ind1 = svadd_m(svnot_z(pg1, cmp1) , gat_ind1, lda_vec);
        gat_ind2 = svadd_m(svnot_z(pg2, cmp2) , gat_ind2, lda_vec);
        gat_ind3 = svadd_m(svnot_z(pg3, cmp3) , gat_ind3, lda_vec);
        gat_ind4 = svadd_m(svnot_z(pg4, cmp4) , gat_ind4, lda_vec);

        svst1(pg1, b, data_vec1);
        svst1_vnum(pg2, b, 1, data_vec2);
        svst1_vnum(pg3, b, 2, data_vec3);
        svst1_vnum(pg4, b, 3, data_vec4);

        b += active;
        offset --;
        vec_off = svsub_z(pg1, vec_off, one_vec);
        cmp1 = svcmpgt(pg1, vec_off, index_neg1);
        cmp2 = svcmpgt(pg2, vec_off, index_neg2);
        cmp3 = svcmpgt(pg3, vec_off, index_neg3);
        cmp4 = svcmpgt(pg4, vec_off, index_neg4);
        
        i--;
    }

    posX += unroll_size;
    posX_vec = svdup_s32(posX);
    j += unroll_size;
    pg1 = svwhilelt_b32((uint32_t)j, (uint32_t)N);
    pg2 = svwhilelt_b32((uint32_t)j+sve_size, (uint32_t)N);
    pg3 = svwhilelt_b32((uint32_t)j+sve_size2, (uint32_t)N);
    pg4 = svwhilelt_b32((uint32_t)j+sve_size3, (uint32_t)N);
    active = svcntp_b32(svptrue_b32(), pg1)+
             svcntp_b32(svptrue_b32(), pg2)+
             svcntp_b32(svptrue_b32(), pg3)+
             svcntp_b32(svptrue_b32(), pg4);
  } while (active>0);

#endif

  return 0;
}
