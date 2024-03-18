/***************************************************************************
Copyright 2024 RIKEN
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

#ifndef __DEF_SVE_ASM_FLOAT
#define __DEF_SVE_ASM_FLOAT

/* fmla(predicated) zds3=zds3+zs1*zs2 p/M */
#define FMLA_ZSP(zds3,pg,zs1,zs2) FMLA_ZSP_base(zds3,pg,zs1,zs2)
#define FMLA_ZSP_base(zds3,pg,zs1,zs2) \
  __asm__ __volatile__("\t\t\tfmla "#zds3".s,"#pg"/M,"#zs1".s,"#zs2".s":::#zs1,#zs2,#zds3,#pg);

/* fmad(predicated) zds1=zs3+zds1*zs2  p/M */
#define FMAD_ZSP(zds1,pg,zs2,zs3) FMAD_ZSP_base(zds1,pg,zs2,zs3)
#define FMAD_ZSP_base(zds1,pg,zs2,zs3) \
	  __asm__ __volatile__("\t\t\tfmad "#zds1".s,"#pg"/M,"#zs2".s,"#zs3".s":::#zds1,#zs2,#zs3,#pg);

/* fmul(unpredicated) zd=zs1*zs2 */
#define FMUL_ZS(zd,zs1,zs2) FMUL_ZS_base(zd,zs1,zs2)
#define FMUL_ZS_base(zd,zs1,zs2) \
  __asm__ __volatile__("\t\t\tfmul "#zd".s,"#zs1".s,"#zs2".s":::#zd,#zs1,#zs2);

/* fmul(predicated) zds1*=zs2  p/M */
#define FMUL_ZSP(zds1,pg,zs2) FMUL_ZSP_base(zds1,pg,zs2)
#define FMUL_ZSP_base(zds1,pg,zs2) \
  __asm__ __volatile__("\t\t\tfmul "#zds1".s,"#pg"/M,"#zds1".s,"#zs2".s":::#zds1,#zs2,#pg);

#define LD1W_ZXI(zt,pg,x1,imm) LD1W_ZXI_base(zt,pg,x1,imm)
#define LD1W_ZXI_base(zt,pg,x1,imm) \
__asm__ __volatile__("\t\t\tld1w {"#zt".s},"#pg"/Z,[%0,#"#imm",MUL VL]"::"r"(x1):#zt,#pg);

#define LD1RW_ZXI(zt,pg,x1,imm) LD1RW_ZXI_base(zt,pg,x1,imm)
#define LD1RW_ZXI_base(zt,pg,x1,imm) \
__asm__ __volatile__("\t\t\tld1rw "#zt".s,"#pg"/Z,[%0,#"#imm"]"::"r"(x1):#zt,#pg);

#define ST1W_ZXI(zt,pg,x1,imm) ST1W_ZXI_base(zt,pg,x1,imm)
#define ST1W_ZXI_base(zt,pg,x1,imm) \
__asm__ __volatile__("\t\t\tst1w {"#zt".s},"#pg",[%0,#"#imm",MUL VL]"::"r"(x1):"memory",#zt,#pg);

#define FNEG_ZSP(zd1,pg,zs2) FNEG_ZSP_base(zd1,pg,zs2)
#define FNEG_ZSP_base(zd1,pg,zs2) \
  __asm__ __volatile__("\t\t\tfneg "#zd1".s,"#pg"/M,"#zs2".s":::#zd1,#zs2,#pg);

#endif  /* __DEF_SVE_ASM_FLOAT */



