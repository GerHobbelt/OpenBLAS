/***************************************************************************
(c) RIKEN 2024, 2024. All rights reserved. sgemm_ncopy_sve_1.c 0.3.26
Copyright 2024 FUJITSU limited
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

#include "common.h"
#include "def_sve_asm.h"
#include <arm_sve.h>

int CNAME(BLASLONG m, BLASLONG n, FLOAT *a, BLASLONG lda, FLOAT *b){
  BLASLONG im, in;
  FLOAT *a_offset, *a_offset_p;
  FLOAT *b_offset, *b_offset_p;
  BLASLONG nvl=svcntw();
  BLASLONG nul=nvl*4;

  a_offset = a;
  a_offset_p = a+lda*2;
  b_offset = b;
  b_offset_p = b+m*2;

  PTRUE_PS(p0);
  for(in=0;in<n;in++) {
    for(im=0;im<m-nul+1;im+=nul) {
      LD1W_ZXI(z0,p0,a_offset,0);
      LD1W_ZXI(z1,p0,a_offset,1);
      LD1W_ZXI(z2,p0,a_offset,2);
      LD1W_ZXI(z3,p0,a_offset,3);

      ST1W_ZXI(z0,p0,b_offset,0);
      ST1W_ZXI(z1,p0,b_offset,1);
      ST1W_ZXI(z2,p0,b_offset,2);
      ST1W_ZXI(z3,p0,b_offset,3);
      PRFM_XI(PLDL2KEEP,a_offset_p,0);
      PRFM_XI(PSTL2KEEP,b_offset_p,0);
      a_offset+=nul;
      b_offset+=nul;
      a_offset_p+=nul;
      b_offset_p+=nul;
    }
    if(m-im>=nvl) {
      for(;im<m-nvl+1;im+=nvl) {
        LD1W_ZXI(z0,p0,a_offset,0);
        ST1W_ZXI(z0,p0,b_offset,0);
        a_offset+=nvl;
        b_offset+=nvl;
        a_offset_p+=nvl;
        b_offset_p+=nvl;
      }
    }
    if(m-im>0) {
      WHILELT_PSX(p1,im,m);
      LD1W_ZXI(z0,p1,a_offset,0);
      ST1W_ZXI(z0,p1,b_offset,0);
      a_offset+=m-im;
      b_offset+=m-im;
      a_offset_p+=m-im;
      b_offset_p+=m-im;
      PRFM_XI(PLDL2KEEP,a_offset_p,-1);
      PRFM_XI(PSTL2KEEP,b_offset_p,-1);
    }
    a_offset+=lda-m; 
    a_offset_p+=lda-m; 
  }

  return 0;
}
