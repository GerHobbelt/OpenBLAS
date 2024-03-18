/***************************************************************************
(c) RIKEN 2024, 2024. All rights reserved. sgemm_kernel_sve_v4x5_1.c 0.3.26
Copyright 2017,2024 FUJITSU limited
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

/*
   A
   z0
   z1

   B
            |  0    +1   +2   +3   +4
  ----------+--------------------------
   0        | z22  z23  z24  z25  z26
   +nvl*1   | z27  z28  z29  z30  z31

   C
           |  0   +1   +2   +3   +4
   --------+-------------------------
   0       | z2   z3   z4   z5   z6
   +nvl*1  | z7   z8   z9   z10  z11
   +nvl*2  | z12  z13  z14  z15  z16
   +nvl*3  | z17  z18  z19  z20  z21

  ALPHA : z0

*/

/* clear C */
#define CZERO_11   DUP_ZSI(z2 ,0)
#define CZERO_12   CZERO_11;DUP_ZSI(z3 ,0)
#define CZERO_13   CZERO_12;DUP_ZSI(z4 ,0)
#define CZERO_14   CZERO_13;DUP_ZSI(z5 ,0)
#define CZERO_15   CZERO_14;DUP_ZSI(z6 ,0)

#define CZERO_21   CZERO_11;DUP_ZSI(z7 ,0)
#define CZERO_22   CZERO_12;DUP_ZSI(z7 ,0);DUP_ZSI(z8 ,0)
#define CZERO_23   CZERO_13;DUP_ZSI(z7 ,0);DUP_ZSI(z8 ,0);DUP_ZSI(z9 ,0)
#define CZERO_24   CZERO_14;DUP_ZSI(z7 ,0);DUP_ZSI(z8 ,0);DUP_ZSI(z9 ,0);DUP_ZSI(z10 ,0)
#define CZERO_25   CZERO_15;DUP_ZSI(z7 ,0);DUP_ZSI(z8 ,0);DUP_ZSI(z9 ,0);DUP_ZSI(z10 ,0);DUP_ZSI(z11 ,0)

#define CZERO_31   CZERO_21;DUP_ZSI(z12 ,0)
#define CZERO_32   CZERO_22;DUP_ZSI(z12 ,0);DUP_ZSI(z13 ,0)
#define CZERO_33   CZERO_23;DUP_ZSI(z12 ,0);DUP_ZSI(z13 ,0);DUP_ZSI(z14 ,0)
#define CZERO_34   CZERO_24;DUP_ZSI(z12 ,0);DUP_ZSI(z13 ,0);DUP_ZSI(z14 ,0);DUP_ZSI(z15 ,0)
#define CZERO_35   CZERO_25;DUP_ZSI(z12 ,0);DUP_ZSI(z13 ,0);DUP_ZSI(z14 ,0);DUP_ZSI(z15 ,0);DUP_ZSI(z16 ,0)

#define CZERO_41   CZERO_31;DUP_ZSI(z17 ,0)
#define CZERO_42   CZERO_32;DUP_ZSI(z17 ,0);DUP_ZSI(z18 ,0)
#define CZERO_43   CZERO_33;DUP_ZSI(z17 ,0);DUP_ZSI(z18 ,0);DUP_ZSI(z19 ,0)
#define CZERO_44   CZERO_34;DUP_ZSI(z17 ,0);DUP_ZSI(z18 ,0);DUP_ZSI(z19 ,0);DUP_ZSI(z20 ,0)
#define CZERO_45   CZERO_35;DUP_ZSI(z17 ,0);DUP_ZSI(z18 ,0);DUP_ZSI(z19 ,0);DUP_ZSI(z20 ,0);DUP_ZSI(z21 ,0)

#define CALPHA_41 LD1W_ZXI(z22,p0,cq0,0); LD1W_ZXI(z23,p0,cq0,1); LD1W_ZXI(z24,p0,cq0,2); LD1W_ZXI(z25,p0,cq0,3); FMAD_ZSP(z2,p0,z0,z22); FMAD_ZSP(z7,p0,z0,z23); FMAD_ZSP(z12,p0,z0,z24); FMAD_ZSP(z17,p0,z0,z25); ST1W_ZXI(z2,p0,cq0,0); ST1W_ZXI(z7,p0,cq0,1); ST1W_ZXI(z12,p0,cq0,2); ST1W_ZXI(z17,p0,cq0,3);

#define CALPHA_42 CALPHA_41;LD1W_ZXI(z26,p0,cq1,0); LD1W_ZXI(z27,p0,cq1,1); LD1W_ZXI(z28,p0,cq1,2); LD1W_ZXI(z29,p0,cq1,3); FMAD_ZSP(z3,p0,z0,z26); FMAD_ZSP(z8,p0,z0,z27); FMAD_ZSP(z13,p0,z0,z28); FMAD_ZSP(z18,p0,z0,z29); ST1W_ZXI(z3,p0,cq1,0); ST1W_ZXI(z8,p0,cq1,1); ST1W_ZXI(z13,p0,cq1,2); ST1W_ZXI(z18,p0,cq1,3);

#define CALPHA_43 CALPHA_42;LD1W_ZXI(z22,p0,cq2,0); LD1W_ZXI(z23,p0,cq2,1); LD1W_ZXI(z24,p0,cq2,2); LD1W_ZXI(z25,p0,cq2,3); FMAD_ZSP(z4,p0,z0,z22); FMAD_ZSP(z9,p0,z0,z23); FMAD_ZSP(z14,p0,z0,z24); FMAD_ZSP(z19,p0,z0,z25); ST1W_ZXI(z4,p0,cq2,0); ST1W_ZXI(z9,p0,cq2,1); ST1W_ZXI(z14,p0,cq2,2); ST1W_ZXI(z19,p0,cq2,3);

#define CALPHA_44 CALPHA_43;LD1W_ZXI(z26,p0,cq3,0); LD1W_ZXI(z27,p0,cq3,1); LD1W_ZXI(z28,p0,cq3,2); LD1W_ZXI(z29,p0,cq3,3); FMAD_ZSP(z5,p0,z0,z26); FMAD_ZSP(z10,p0,z0,z27); FMAD_ZSP(z15,p0,z0,z28); FMAD_ZSP(z20,p0,z0,z29); ST1W_ZXI(z5,p0,cq3,0); ST1W_ZXI(z10,p0,cq3,1); ST1W_ZXI(z15,p0,cq3,2); ST1W_ZXI(z20,p0,cq3,3);

#define CALPHA_45 CALPHA_44; LD1W_ZXI(z22,p0,cq4,0); LD1W_ZXI(z23,p0,cq4,1); LD1W_ZXI(z24,p0,cq4,2); LD1W_ZXI(z25,p0,cq4,3); FMAD_ZSP(z6,p0,z0,z22); FMAD_ZSP(z11,p0,z0,z23); FMAD_ZSP(z16,p0,z0,z24); FMAD_ZSP(z21,p0,z0,z25); ST1W_ZXI(z6,p0,cq4,0); ST1W_ZXI(z11,p0,cq4,1); ST1W_ZXI(z16,p0,cq4,2); ST1W_ZXI(z21,p0,cq4,3);


#define CALPHA_P_11 LD1W_ZXI(z22,p1,cq0,0); FMAD_ZSP(z2,p1,z0,z22);ST1W_ZXI(z2,p1,cq0,0);
#define CALPHA_P_21 LD1W_ZXI(z22,p0,cq0,0); LD1W_ZXI(z23,p1,cq0,1); FMAD_ZSP(z2,p0,z0,z22); FMAD_ZSP(z7,p1,z0,z23); ST1W_ZXI(z2,p0,cq0,0); ST1W_ZXI(z7,p1,cq0,1);
#define CALPHA_P_31 LD1W_ZXI(z22,p0,cq0,0); LD1W_ZXI(z23,p0,cq0,1); LD1W_ZXI(z24,p1,cq0,2); FMAD_ZSP(z2,p0,z0,z22); FMAD_ZSP(z7,p0,z0,z23); FMAD_ZSP(z12,p1,z0,z24); ST1W_ZXI(z2,p0,cq0,0); ST1W_ZXI(z7,p0,cq0,1); ST1W_ZXI(z12,p1,cq0,2);
#define CALPHA_P_41 LD1W_ZXI(z22,p0,cq0,0); LD1W_ZXI(z23,p0,cq0,1); LD1W_ZXI(z24,p0,cq0,2); LD1W_ZXI(z25,p1,cq0,3); FMAD_ZSP(z2,p0,z0,z22); FMAD_ZSP(z7,p0,z0,z23); FMAD_ZSP(z12,p0,z0,z24); FMAD_ZSP(z17,p1,z0,z25); ST1W_ZXI(z2,p0,cq0,0); ST1W_ZXI(z7,p0,cq0,1); ST1W_ZXI(z12,p0,cq0,2); ST1W_ZXI(z17,p1,cq0,3);

#define CALPHA_P_12 CALPHA_P_11;LD1W_ZXI(z26,p1,cq1,0); FMAD_ZSP(z3,p1,z0,z26); ST1W_ZXI(z3,p1,cq1,0);
#define CALPHA_P_22 CALPHA_P_21;LD1W_ZXI(z26,p0,cq1,0); LD1W_ZXI(z27,p1,cq1,1); FMAD_ZSP(z3,p0,z0,z26); FMAD_ZSP(z8,p1,z0,z27); ST1W_ZXI(z3,p0,cq1,0); ST1W_ZXI(z8,p1,cq1,1);
#define CALPHA_P_32 CALPHA_P_31;LD1W_ZXI(z26,p0,cq1,0); LD1W_ZXI(z27,p0,cq1,1); LD1W_ZXI(z28,p1,cq1,2); FMAD_ZSP(z3,p0,z0,z26); FMAD_ZSP(z8,p0,z0,z27); FMAD_ZSP(z13,p1,z0,z28); ST1W_ZXI(z3,p0,cq1,0); ST1W_ZXI(z8,p0,cq1,1); ST1W_ZXI(z13,p1,cq1,2);
#define CALPHA_P_42 CALPHA_P_41;LD1W_ZXI(z26,p0,cq1,0); LD1W_ZXI(z27,p0,cq1,1); LD1W_ZXI(z28,p0,cq1,2); LD1W_ZXI(z29,p1,cq1,3); FMAD_ZSP(z3,p0,z0,z26); FMAD_ZSP(z8,p0,z0,z27); FMAD_ZSP(z13,p0,z0,z28); FMAD_ZSP(z18,p1,z0,z29); ST1W_ZXI(z3,p0,cq1,0); ST1W_ZXI(z8,p0,cq1,1); ST1W_ZXI(z13,p0,cq1,2); ST1W_ZXI(z18,p1,cq1,3);

#define CALPHA_P_13 CALPHA_P_12;LD1W_ZXI(z22,p1,cq2,0); FMAD_ZSP(z4,p1,z0,z22); ST1W_ZXI(z4,p1,cq2,0);
#define CALPHA_P_23 CALPHA_P_22;LD1W_ZXI(z22,p0,cq2,0); LD1W_ZXI(z23,p1,cq2,1); FMAD_ZSP(z4,p0,z0,z22); FMAD_ZSP(z9,p1,z0,z23); ST1W_ZXI(z4,p0,cq2,0); ST1W_ZXI(z9,p1,cq2,1);
#define CALPHA_P_33 CALPHA_P_32;LD1W_ZXI(z22,p0,cq2,0); LD1W_ZXI(z23,p0,cq2,1); LD1W_ZXI(z24,p1,cq2,2); FMAD_ZSP(z4,p0,z0,z22); FMAD_ZSP(z9,p0,z0,z23); FMAD_ZSP(z14,p1,z0,z24); ST1W_ZXI(z4,p0,cq2,0); ST1W_ZXI(z9,p0,cq2,1); ST1W_ZXI(z14,p1,cq2,2);
#define CALPHA_P_43 CALPHA_P_42;LD1W_ZXI(z22,p0,cq2,0); LD1W_ZXI(z23,p0,cq2,1); LD1W_ZXI(z24,p0,cq2,2); LD1W_ZXI(z25,p1,cq2,3); FMAD_ZSP(z4,p0,z0,z22); FMAD_ZSP(z9,p0,z0,z23); FMAD_ZSP(z14,p0,z0,z24); FMAD_ZSP(z19,p1,z0,z25); ST1W_ZXI(z4,p0,cq2,0); ST1W_ZXI(z9,p0,cq2,1); ST1W_ZXI(z14,p0,cq2,2); ST1W_ZXI(z19,p1,cq2,3);

#define CALPHA_P_14 CALPHA_P_13;LD1W_ZXI(z26,p1,cq3,0); FMAD_ZSP(z5,p1,z0,z26); ST1W_ZXI(z5,p1,cq3,0);
#define CALPHA_P_24 CALPHA_P_23;LD1W_ZXI(z26,p0,cq3,0); LD1W_ZXI(z27,p1,cq3,1); FMAD_ZSP(z5,p0,z0,z26); FMAD_ZSP(z10,p1,z0,z27); ST1W_ZXI(z5,p0,cq3,0); ST1W_ZXI(z10,p1,cq3,1);
#define CALPHA_P_34 CALPHA_P_33;LD1W_ZXI(z26,p0,cq3,0); LD1W_ZXI(z27,p0,cq3,1); LD1W_ZXI(z28,p1,cq3,2); FMAD_ZSP(z5,p0,z0,z26); FMAD_ZSP(z10,p0,z0,z27); FMAD_ZSP(z15,p1,z0,z28); ST1W_ZXI(z5,p0,cq3,0); ST1W_ZXI(z10,p0,cq3,1); ST1W_ZXI(z15,p1,cq3,2);
#define CALPHA_P_44 CALPHA_P_43;LD1W_ZXI(z26,p0,cq3,0); LD1W_ZXI(z27,p0,cq3,1); LD1W_ZXI(z28,p0,cq3,2); LD1W_ZXI(z29,p1,cq3,3); FMAD_ZSP(z5,p0,z0,z26); FMAD_ZSP(z10,p0,z0,z27); FMAD_ZSP(z15,p0,z0,z28); FMAD_ZSP(z20,p1,z0,z29); ST1W_ZXI(z5,p0,cq3,0); ST1W_ZXI(z10,p0,cq3,1); ST1W_ZXI(z15,p0,cq3,2); ST1W_ZXI(z20,p1,cq3,3);

#define CALPHA_P_15 CALPHA_P_14; LD1W_ZXI(z22,p1,cq4,0); FMAD_ZSP(z6,p1,z0,z22); ST1W_ZXI(z6,p1,cq4,0);
#define CALPHA_P_25 CALPHA_P_24; LD1W_ZXI(z22,p0,cq4,0); LD1W_ZXI(z23,p1,cq4,1); FMAD_ZSP(z6,p0,z0,z22); FMAD_ZSP(z11,p1,z0,z23); ST1W_ZXI(z6,p0,cq4,0); ST1W_ZXI(z11,p1,cq4,1);
#define CALPHA_P_35 CALPHA_P_34; LD1W_ZXI(z22,p0,cq4,0); LD1W_ZXI(z23,p0,cq4,1); LD1W_ZXI(z24,p1,cq4,2); FMAD_ZSP(z6,p0,z0,z22); FMAD_ZSP(z11,p0,z0,z23); FMAD_ZSP(z16,p1,z0,z24); ST1W_ZXI(z6,p0,cq4,0); ST1W_ZXI(z11,p0,cq4,1); ST1W_ZXI(z16,p1,cq4,2);
#define CALPHA_P_45 CALPHA_P_44; LD1W_ZXI(z22,p0,cq4,0); LD1W_ZXI(z23,p0,cq4,1); LD1W_ZXI(z24,p0,cq4,2); LD1W_ZXI(z25,p1,cq4,3); FMAD_ZSP(z6,p0,z0,z22); FMAD_ZSP(z11,p0,z0,z23); FMAD_ZSP(z16,p0,z0,z24); FMAD_ZSP(z21,p1,z0,z25); ST1W_ZXI(z6,p0,cq4,0); ST1W_ZXI(z11,p0,cq4,1); ST1W_ZXI(z16,p0,cq4,2); ST1W_ZXI(z21,p1,cq4,3);

int CNAME(BLASLONG m,BLASLONG n,BLASLONG k,FLOAT alpha0, FLOAT* ap,FLOAT* b,FLOAT* cp, BLASLONG ldc )
{
	BLASLONG ldb=k;
	BLASLONG im,in,ik,mm;
	FLOAT alpha=alpha0;
	FLOAT *alphap=&alpha;
	FLOAT *a,*c;
	FLOAT *aa;
	FLOAT *bq0,*bq1,*bq2,*bq3,*bq4;
	FLOAT *bs0,*bs1,*bs2,*bs3,*bs4;
	FLOAT *cq0,*cq1,*cq2,*cq3,*cq4;
	BLASLONG nul;

	a=ap;
	c=cp;
	BLASLONG nvl=svcntw();
        PTRUE_PS(p0);

	nul=nvl*4;
	cq0=c;
	cq1=c+ldc;
	cq2=c+ldc*2;
	cq3=c+ldc*3;
	cq4=c+ldc*4;

	bs0=b;
	bs1=b+ldb;
	bs2=b+ldb*2;
	bs3=b+ldb*3;
	bs4=b+ldb*4;

	in=0;
	for(;in<n-5+1;in+=5) {
		aa=a;
		im=0;
		for(;im<m-nul+1;im+=nul) {
			PRFM_XI(PSTL1KEEP,cq0,256*2);
			PRFM_XI(PSTL1KEEP,cq1,256*2);
			PRFM_XI(PSTL1KEEP,cq2,256*2);
			PRFM_XI(PSTL1KEEP,cq3,256*2);
			PRFM_XI(PSTL1KEEP,cq4,256*2);

			CZERO_45;

			PRFM_XI(PLDL1KEEP,aa,0);
			PRFM_XI(PLDL1KEEP,aa,256*1);
			PRFM_XI(PLDL1KEEP,aa,256*2);
			PRFM_XI(PLDL1KEEP,aa,256*3);
			PRFM_XI(PLDL1KEEP,aa,256*4);
			PRFM_XI(PLDL1KEEP,aa,256*5);
			PRFM_XI(PLDL1KEEP,aa,256*6);
			PRFM_XI(PLDL1KEEP,aa,256*7);
			PRFM_XI(PLDL1KEEP,aa,256*8);


			bq0=bs0;
			bq1=bs1;
			bq2=bs2;
			bq3=bs3;
			bq4=bs4;

                        if(im<k) {
                                PRFM_XXLSL3(PLDL2KEEP,bq0+ldb*5,im/2);
                                PRFM_XXLSL3(PLDL2KEEP,bq1+ldb*5,im/2);
                                PRFM_XXLSL3(PLDL2KEEP,bq2+ldb*5,im/2);
                                PRFM_XXLSL3(PLDL2KEEP,bq3+ldb*5,im/2);
                                PRFM_XXLSL3(PLDL2KEEP,bq4+ldb*5,im/2);
                        }
			ik=0;
			if(k>=2) { 
                                if(im==0) {
					LD1RW_ZXI(z22,p0,bq0,0);
					LD1RW_ZXI(z23,p0,bq1,0);
					LD1RW_ZXI(z24,p0,bq2,0);
					LD1RW_ZXI(z25,p0,bq3,0);
					LD1RW_ZXI(z26,p0,bq4,0);

					LD1RW_ZXI(z27,p0,bq0,SIZE);
					LD1RW_ZXI(z28,p0,bq1,SIZE);

                                        PRFM_XI(PLDL1KEEP,bq0,256);
                                        PRFM_XI(PLDL1KEEP,bq1,256);
                                        PRFM_XI(PLDL1KEEP,bq2,256);
                                        PRFM_XI(PLDL1KEEP,bq3,256);
                                        PRFM_XI(PLDL1KEEP,bq4,256);

					LD1W_ZXI(z0,p0,aa,0);

					ik+=2;

					for(;ik<k-1;) {
						LD1W_ZXI(z1,p0,aa,1);
						PRFM_XI(PLDL1KEEP,aa,256*9);

						FMLA_ZSP( z2,p0,z0,z22);
						FMLA_ZSP( z3,p0,z0,z23);

						LD1RW_ZXI(z29,p0,bq2,SIZE);
						LD1RW_ZXI(z30,p0,bq3,SIZE);

						FMLA_ZSP( z4,p0,z0,z24);
						FMLA_ZSP( z5,p0,z0,z25);
						FMLA_ZSP( z6,p0,z0,z26);
						FMLA_ZSP( z7,p0,z1,z22);

						LD1W_ZXI(z0,p0,aa,2);
						LD1RW_ZXI(z31,p0,bq4,SIZE);

						FMLA_ZSP( z8,p0,z1,z23);
						FMLA_ZSP( z9,p0,z1,z24);

                                                PRFM_XI(PLDL1KEEP,bq0,256);
                                                PRFM_XI(PLDL1KEEP,bq1,256);

						FMLA_ZSP( z10,p0,z1,z25);
						FMLA_ZSP( z11,p0,z1,z26);

						LD1W_ZXI(z1,p0,aa,3);
						bq0+=2;

						FMLA_ZSP( z12,p0,z0,z22);
						FMLA_ZSP( z13,p0,z0,z23);

						bq1+=2;
						bq2+=2;

						FMLA_ZSP( z14,p0,z0,z24);
						FMLA_ZSP( z15,p0,z0,z25);
						FMLA_ZSP( z16,p0,z0,z26);
						FMLA_ZSP( z17,p0,z1,z22);

						LD1W_ZXI(z0,p0,aa,4);
						bq3+=2;

						FMLA_ZSP( z18,p0,z1,z23);
						FMLA_ZSP( z19,p0,z1,z24);

						LD1RW_ZXI(z22,p0,bq0,0);
						LD1RW_ZXI(z23,p0,bq1,0);

						FMLA_ZSP( z20,p0,z1,z25);
						FMLA_ZSP( z21,p0,z1,z26);

						LD1W_ZXI(z1,p0,aa,5);
						PRFM_XI(PLDL1KEEP,aa,256*10);

						FMLA_ZSP( z2,p0,z0,z27);
						FMLA_ZSP( z3,p0,z0,z28);

						LD1RW_ZXI(z24,p0,bq2,0);
						LD1RW_ZXI(z25,p0,bq3,0);

						FMLA_ZSP( z4,p0,z0,z29);
						FMLA_ZSP( z5,p0,z0,z30);
						FMLA_ZSP( z6,p0,z0,z31);
						FMLA_ZSP( z7,p0,z1,z27);

						LD1W_ZXI(z0,p0,aa,6);
						LD1RW_ZXI(z26,p0,bq4,SIZE*2);

						FMLA_ZSP( z8,p0,z1,z28);
						FMLA_ZSP( z9,p0,z1,z29);

                                                PRFM_XI(PLDL1KEEP,bq2,256);
                                                PRFM_XI(PLDL1KEEP,bq3,256);

						FMLA_ZSP( z10,p0,z1,z30);
						FMLA_ZSP( z11,p0,z1,z31);

						LD1W_ZXI(z1,p0,aa,7);
						aa+=nul*2;

						FMLA_ZSP( z12,p0,z0,z27);
						FMLA_ZSP( z13,p0,z0,z28);


                                                PRFM_XI(PLDL1KEEP,bq4,256);

						FMLA_ZSP( z14,p0,z0,z29);
						FMLA_ZSP( z15,p0,z0,z30);
						FMLA_ZSP( z16,p0,z0,z31);
						FMLA_ZSP( z17,p0,z1,z27);

						LD1W_ZXI(z0,p0,aa,0);
						bq4+=2;
						ik+=2;

						FMLA_ZSP( z18,p0,z1,z28);
						FMLA_ZSP( z19,p0,z1,z29);

						LD1RW_ZXI(z27,p0,bq0,SIZE);
						LD1RW_ZXI(z28,p0,bq1,SIZE);

						FMLA_ZSP( z20,p0,z1,z30);
						FMLA_ZSP( z21,p0,z1,z31);
					}
				}
				else {

					LD1RW_ZXI(z22,p0,bq0,0);
					LD1RW_ZXI(z23,p0,bq1,0);
					LD1RW_ZXI(z24,p0,bq2,0);
					LD1RW_ZXI(z25,p0,bq3,0);
					LD1RW_ZXI(z26,p0,bq4,0);

					LD1RW_ZXI(z27,p0,bq0,SIZE);
					LD1RW_ZXI(z28,p0,bq1,SIZE);

					LD1W_ZXI(z0,p0,aa,0);

					ik+=2;

					for(;ik<k-1;) {
						LD1W_ZXI(z1,p0,aa,1);
						PRFM_XI(PLDL1KEEP,aa,256*9);

						FMLA_ZSP( z2,p0,z0,z22);
						FMLA_ZSP( z3,p0,z0,z23);

						LD1RW_ZXI(z29,p0,bq2,SIZE);
						LD1RW_ZXI(z30,p0,bq3,SIZE);

						FMLA_ZSP( z4,p0,z0,z24);
						FMLA_ZSP( z5,p0,z0,z25);
						FMLA_ZSP( z6,p0,z0,z26);
						FMLA_ZSP( z7,p0,z1,z22);

						LD1W_ZXI(z0,p0,aa,2);
						LD1RW_ZXI(z31,p0,bq4,SIZE);

						FMLA_ZSP( z8,p0,z1,z23);
						FMLA_ZSP( z9,p0,z1,z24);
						FMLA_ZSP( z10,p0,z1,z25);
						FMLA_ZSP( z11,p0,z1,z26);

						LD1W_ZXI(z1,p0,aa,3);
						bq0+=2;

						FMLA_ZSP( z12,p0,z0,z22);
						FMLA_ZSP( z13,p0,z0,z23);

						bq1+=2;
						bq2+=2;

						FMLA_ZSP( z14,p0,z0,z24);
						FMLA_ZSP( z15,p0,z0,z25);
						FMLA_ZSP( z16,p0,z0,z26);
						FMLA_ZSP( z17,p0,z1,z22);

						LD1W_ZXI(z0,p0,aa,4);
						bq3+=2;

						FMLA_ZSP( z18,p0,z1,z23);
						FMLA_ZSP( z19,p0,z1,z24);

						LD1RW_ZXI(z22,p0,bq0,0);
						LD1RW_ZXI(z23,p0,bq1,0);

						FMLA_ZSP( z20,p0,z1,z25);
						FMLA_ZSP( z21,p0,z1,z26);

						LD1W_ZXI(z1,p0,aa,5);
						PRFM_XI(PLDL1KEEP,aa,256*10);

						FMLA_ZSP( z2,p0,z0,z27);
						FMLA_ZSP( z3,p0,z0,z28);

						LD1RW_ZXI(z24,p0,bq2,0);
						LD1RW_ZXI(z25,p0,bq3,0);

						FMLA_ZSP( z4,p0,z0,z29);
						FMLA_ZSP( z5,p0,z0,z30);
						FMLA_ZSP( z6,p0,z0,z31);
						FMLA_ZSP( z7,p0,z1,z27);

						LD1W_ZXI(z0,p0,aa,6);
						LD1RW_ZXI(z26,p0,bq4,SIZE*2);

						FMLA_ZSP( z8,p0,z1,z28);
						FMLA_ZSP( z9,p0,z1,z29);
						FMLA_ZSP( z10,p0,z1,z30);
						FMLA_ZSP( z11,p0,z1,z31);

						LD1W_ZXI(z1,p0,aa,7);
						aa+=nul*2;

						FMLA_ZSP( z12,p0,z0,z27);
						FMLA_ZSP( z13,p0,z0,z28);

						bq4+=2;
						ik+=2;

						FMLA_ZSP( z14,p0,z0,z29);
						FMLA_ZSP( z15,p0,z0,z30);
						FMLA_ZSP( z16,p0,z0,z31);
						FMLA_ZSP( z17,p0,z1,z27);

						LD1W_ZXI(z0,p0,aa,0);
						NOP;

						FMLA_ZSP( z18,p0,z1,z28);
						FMLA_ZSP( z19,p0,z1,z29);

						LD1RW_ZXI(z27,p0,bq0,SIZE);
						LD1RW_ZXI(z28,p0,bq1,SIZE);

						FMLA_ZSP( z20,p0,z1,z30);
						FMLA_ZSP( z21,p0,z1,z31);
					}
				}
				LD1W_ZXI(z1,p0,aa,1);

				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z3,p0,z0,z23);

				LD1RW_ZXI(z29,p0,bq2,SIZE);
				LD1RW_ZXI(z30,p0,bq3,SIZE);

				FMLA_ZSP( z4,p0,z0,z24);
				FMLA_ZSP( z5,p0,z0,z25);
				FMLA_ZSP( z6,p0,z0,z26);

				LD1W_ZXI(z0,p0,aa,2);
				LD1RW_ZXI(z31,p0,bq4,SIZE);

				FMLA_ZSP( z7 ,p0,z1,z22);
				FMLA_ZSP( z8 ,p0,z1,z23);
				FMLA_ZSP( z9 ,p0,z1,z24);
				FMLA_ZSP( z10,p0,z1,z25);
				FMLA_ZSP( z11,p0,z1,z26);

				LD1W_ZXI(z1,p0,aa,3);

				FMLA_ZSP( z12,p0,z0,z22);
				FMLA_ZSP( z13,p0,z0,z23);
				FMLA_ZSP( z14,p0,z0,z24);
				FMLA_ZSP( z15,p0,z0,z25);
				FMLA_ZSP( z16,p0,z0,z26);

				LD1W_ZXI(z0,p0,aa,4);

				FMLA_ZSP( z17,p0,z1,z22);
				FMLA_ZSP( z18,p0,z1,z23);
				FMLA_ZSP( z19,p0,z1,z24);
				FMLA_ZSP( z20,p0,z1,z25);
				FMLA_ZSP( z21,p0,z1,z26);

				LD1W_ZXI(z1,p0,aa,5);

				FMLA_ZSP( z2,p0,z0,z27);
				FMLA_ZSP( z3,p0,z0,z28);
				FMLA_ZSP( z4,p0,z0,z29);
				FMLA_ZSP( z5,p0,z0,z30);
				FMLA_ZSP( z6,p0,z0,z31);


				LD1W_ZXI(z0,p0,aa,6);

				FMLA_ZSP( z7,p0,z1,z27);
				FMLA_ZSP( z8,p0,z1,z28);
				FMLA_ZSP( z9,p0,z1,z29);
				FMLA_ZSP( z10,p0,z1,z30);
				FMLA_ZSP( z11,p0,z1,z31);

				LD1W_ZXI(z1,p0,aa,7);
				bq0+=2;
				bq1+=2;
				bq2+=2;
				bq3+=2;
				bq4+=2;
				aa+=nul*2;

				if(k-ik>=1) {
					LD1RW_ZXI(z22,p0,bq0,0);
					LD1RW_ZXI(z23,p0,bq1,0);

					FMLA_ZSP( z12,p0,z0,z27);
					FMLA_ZSP( z13,p0,z0,z28);

					LD1RW_ZXI(z24,p0,bq2,0);
					LD1RW_ZXI(z25,p0,bq3,0);

					FMLA_ZSP( z14,p0,z0,z29);
					FMLA_ZSP( z15,p0,z0,z30);
					FMLA_ZSP( z16,p0,z0,z31);

					LD1RW_ZXI(z26,p0,bq4,0);
					LD1W_ZXI(z0,p0,aa,0);

					FMLA_ZSP( z17,p0,z1,z27);
					FMLA_ZSP( z18,p0,z1,z28);
					FMLA_ZSP( z19,p0,z1,z29);
					FMLA_ZSP( z20,p0,z1,z30);
					FMLA_ZSP( z21,p0,z1,z31);


					LD1W_ZXI(z1,p0,aa,1);
				}
				else {
					FMLA_ZSP( z12,p0,z0,z27);
					FMLA_ZSP( z13,p0,z0,z28);
					FMLA_ZSP( z14,p0,z0,z29);
					FMLA_ZSP( z15,p0,z0,z30);
					FMLA_ZSP( z16,p0,z0,z31);
					FMLA_ZSP( z17,p0,z1,z27);
					FMLA_ZSP( z18,p0,z1,z28);
					FMLA_ZSP( z19,p0,z1,z29);
					FMLA_ZSP( z20,p0,z1,z30);
					FMLA_ZSP( z21,p0,z1,z31);
				}
			}
			else if(k==1) {
				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);
				LD1RW_ZXI(z24,p0,bq2,0);
				LD1RW_ZXI(z25,p0,bq3,0);
				LD1RW_ZXI(z26,p0,bq4,0);

				LD1W_ZXI(z0,p0,aa,0);
				LD1W_ZXI(z1,p0,aa,1);

			}

			if(k-ik>=1) {

				FMLA_ZSP( z2, p0,z0,z22);
				FMLA_ZSP( z3, p0,z0,z23);
				FMLA_ZSP( z4, p0,z0,z24);
				FMLA_ZSP( z5, p0,z0,z25);
				FMLA_ZSP( z6, p0,z0,z26);

				LD1W_ZXI(z0,p0,aa,2);

				FMLA_ZSP( z7, p0,z1,z22);
				FMLA_ZSP( z8, p0,z1,z23);
				FMLA_ZSP( z9, p0,z1,z24);
				FMLA_ZSP( z10,p0,z1,z25);
				FMLA_ZSP( z11,p0,z1,z26);

				LD1W_ZXI(z1,p0,aa,3);

				FMLA_ZSP( z12,p0,z0,z22);
				FMLA_ZSP( z13,p0,z0,z23);
				FMLA_ZSP( z14,p0,z0,z24);
				FMLA_ZSP( z15,p0,z0,z25);
				FMLA_ZSP( z16,p0,z0,z26);
				FMLA_ZSP( z17,p0,z1,z22);
				FMLA_ZSP( z18,p0,z1,z23);
				FMLA_ZSP( z19,p0,z1,z24);
				FMLA_ZSP( z20,p0,z1,z25);
				FMLA_ZSP( z21,p0,z1,z26);

				aa+=nul;
				bq0+=1;
				bq1+=1;
				bq2+=1;
				bq3+=1;
				bq4+=1;
				ik++;
			}
			LD1RW_ZXI(z0,p0,alphap,0);
			CALPHA_45;
			cq0+=nul;
			cq1+=nul;
			cq2+=nul;
			cq3+=nul;
			cq4+=nul;
		}

		if(m-im>nvl*3) {
			mm=m-im;
			im+=nvl*3;
			WHILELT_PSX(p1,im,m);
			mm=m-im;
			im=m;
			PRFM_XI(PSTL1KEEP,cq0,256*2);
			PRFM_XI(PSTL1KEEP,cq1,256*2);
			PRFM_XI(PSTL1KEEP,cq2,256*2);
			PRFM_XI(PSTL1KEEP,cq3,256*2);
			PRFM_XI(PSTL1KEEP,cq4,256*2);

			CZERO_45;

			PRFM_XI(PLDL1KEEP,aa,0);
			PRFM_XI(PLDL1KEEP,aa,256*1);
			PRFM_XI(PLDL1KEEP,aa,256*2);
			PRFM_XI(PLDL1KEEP,aa,256*3);
			PRFM_XI(PLDL1KEEP,aa,256*4);
			PRFM_XI(PLDL1KEEP,aa,256*5);
			PRFM_XI(PLDL1KEEP,aa,256*6);
			PRFM_XI(PLDL1KEEP,aa,256*7);
			PRFM_XI(PLDL1KEEP,aa,256*8);

			bq0=bs0;
			bq1=bs1;
			bq2=bs2;
			bq3=bs3;
			bq4=bs4;
			ik=0;
			if(k>=2) { 

				if(im==0) {
					LD1RW_ZXI(z22,p0,bq0,0);
					LD1RW_ZXI(z23,p0,bq1,0);
					LD1RW_ZXI(z24,p0,bq2,0);
					LD1RW_ZXI(z25,p0,bq3,0);
					LD1RW_ZXI(z26,p0,bq4,0);

					LD1RW_ZXI(z27,p0,bq0,SIZE);
					LD1RW_ZXI(z28,p0,bq1,SIZE);

					LD1W_ZXI(z0,p0,aa,0);

					ik+=2;

					for(;ik<k-1;) {
						LD1W_ZXI(z1,p0,aa,1);
						PRFM_XI(PLDL1KEEP,aa,256*9);

						FMLA_ZSP( z2,p0,z0,z22);
						FMLA_ZSP( z3,p0,z0,z23);

						LD1RW_ZXI(z29,p0,bq2,SIZE);
						LD1RW_ZXI(z30,p0,bq3,SIZE);

						FMLA_ZSP( z4,p0,z0,z24);
						FMLA_ZSP( z5,p0,z0,z25);
						FMLA_ZSP( z6,p0,z0,z26);
						FMLA_ZSP( z7,p0,z1,z22);

						LD1W_ZXI(z0,p0,aa,2);
						LD1RW_ZXI(z31,p0,bq4,SIZE);

						FMLA_ZSP( z8,p0,z1,z23);
						FMLA_ZSP( z9,p0,z1,z24);

                                                PRFM_XI(PLDL1KEEP,bq0,256);
                                                PRFM_XI(PLDL1KEEP,bq1,256);

						FMLA_ZSP( z10,p0,z1,z25);
						FMLA_ZSP( z11,p0,z1,z26);

						LD1W_ZXI(z1,p1,aa,3);
						aa+=nvl*3+mm;
						bq0+=2;

						FMLA_ZSP( z12,p0,z0,z22);
						FMLA_ZSP( z13,p0,z0,z23);

						bq1+=2;
						bq2+=2;

						FMLA_ZSP( z14,p0,z0,z24);
						FMLA_ZSP( z15,p0,z0,z25);
						FMLA_ZSP( z16,p0,z0,z26);
						FMLA_ZSP( z17,p1,z1,z22);

						LD1W_ZXI(z0,p0,aa,0);
						bq3+=2;

						FMLA_ZSP( z18,p1,z1,z23);
						FMLA_ZSP( z19,p1,z1,z24);

						LD1RW_ZXI(z22,p0,bq0,0);
						LD1RW_ZXI(z23,p0,bq1,0);

						FMLA_ZSP( z20,p1,z1,z25);
						FMLA_ZSP( z21,p1,z1,z26);

						LD1W_ZXI(z1,p0,aa,1);
						PRFM_XI(PLDL1KEEP,aa,256*9);

						FMLA_ZSP( z2,p0,z0,z27);
						FMLA_ZSP( z3,p0,z0,z28);

						LD1RW_ZXI(z24,p0,bq2,0);
						LD1RW_ZXI(z25,p0,bq3,0);

						FMLA_ZSP( z4,p0,z0,z29);
						FMLA_ZSP( z5,p0,z0,z30);

                                                PRFM_XI(PLDL1KEEP,bq2,256);
                                                PRFM_XI(PLDL1KEEP,bq3,256);

						FMLA_ZSP( z6,p0,z0,z31);
						FMLA_ZSP( z7,p0,z1,z27);

						LD1W_ZXI(z0,p0,aa,2);
						LD1RW_ZXI(z26,p0,bq4,SIZE*2);

						FMLA_ZSP( z8,p0,z1,z28);
						FMLA_ZSP( z9,p0,z1,z29);
						FMLA_ZSP( z10,p0,z1,z30);
						FMLA_ZSP( z11,p0,z1,z31);

						LD1W_ZXI(z1,p1,aa,3);
						aa+=nvl*3+mm;

						FMLA_ZSP( z12,p0,z0,z27);
						FMLA_ZSP( z13,p0,z0,z28);

						bq4+=2;
						ik+=2;

						FMLA_ZSP( z14,p0,z0,z29);
						FMLA_ZSP( z15,p0,z0,z30);

                                                PRFM_XI(PLDL1KEEP,bq4,256);

						FMLA_ZSP( z16,p0,z0,z31);
						FMLA_ZSP( z17,p1,z1,z27);

						LD1W_ZXI(z0,p0,aa,0);

						FMLA_ZSP( z18,p1,z1,z28);
						FMLA_ZSP( z19,p1,z1,z29);

						LD1RW_ZXI(z27,p0,bq0,SIZE);
						LD1RW_ZXI(z28,p0,bq1,SIZE);

						FMLA_ZSP( z20,p1,z1,z30);
						FMLA_ZSP( z21,p1,z1,z31);
					}
				}
				else {
					LD1RW_ZXI(z22,p0,bq0,0);
					LD1RW_ZXI(z23,p0,bq1,0);
					LD1RW_ZXI(z24,p0,bq2,0);
					LD1RW_ZXI(z25,p0,bq3,0);
					LD1RW_ZXI(z26,p0,bq4,0);

					LD1RW_ZXI(z27,p0,bq0,SIZE);
					LD1RW_ZXI(z28,p0,bq1,SIZE);

					LD1W_ZXI(z0,p0,aa,0);

					ik+=2;

					for(;ik<k-1;) {
						LD1W_ZXI(z1,p0,aa,1);
						PRFM_XI(PLDL1KEEP,aa,256*9);

						FMLA_ZSP( z2,p0,z0,z22);
						FMLA_ZSP( z3,p0,z0,z23);

						LD1RW_ZXI(z29,p0,bq2,SIZE);
						LD1RW_ZXI(z30,p0,bq3,SIZE);

						FMLA_ZSP( z4,p0,z0,z24);
						FMLA_ZSP( z5,p0,z0,z25);
						FMLA_ZSP( z6,p0,z0,z26);
						FMLA_ZSP( z7,p0,z1,z22);

						LD1W_ZXI(z0,p0,aa,2);
						LD1RW_ZXI(z31,p0,bq4,SIZE);

						FMLA_ZSP( z8,p0,z1,z23);
						FMLA_ZSP( z9,p0,z1,z24);

						FMLA_ZSP( z10,p0,z1,z25);
						FMLA_ZSP( z11,p0,z1,z26);

						LD1W_ZXI(z1,p1,aa,3);
						aa+=nvl*3+mm;
						bq0+=2;

						FMLA_ZSP( z12,p0,z0,z22);
						FMLA_ZSP( z13,p0,z0,z23);

						bq1+=2;
						bq2+=2;

						FMLA_ZSP( z14,p0,z0,z24);
						FMLA_ZSP( z15,p0,z0,z25);
						FMLA_ZSP( z16,p0,z0,z26);
						FMLA_ZSP( z17,p1,z1,z22);

						LD1W_ZXI(z0,p0,aa,0);
						bq3+=2;

						FMLA_ZSP( z18,p1,z1,z23);
						FMLA_ZSP( z19,p1,z1,z24);

						LD1RW_ZXI(z22,p0,bq0,0);
						LD1RW_ZXI(z23,p0,bq1,0);

						FMLA_ZSP( z20,p1,z1,z25);
						FMLA_ZSP( z21,p1,z1,z26);

						LD1W_ZXI(z1,p0,aa,1);
						PRFM_XI(PLDL1KEEP,aa,256*9);

						FMLA_ZSP( z2,p0,z0,z27);
						FMLA_ZSP( z3,p0,z0,z28);

						LD1RW_ZXI(z24,p0,bq2,0);
						LD1RW_ZXI(z25,p0,bq3,0);

						FMLA_ZSP( z4,p0,z0,z29);
						FMLA_ZSP( z5,p0,z0,z30);
						FMLA_ZSP( z6,p0,z0,z31);
						FMLA_ZSP( z7,p0,z1,z27);

						LD1W_ZXI(z0,p0,aa,2);
						LD1RW_ZXI(z26,p0,bq4,SIZE*2);

						FMLA_ZSP( z8,p0,z1,z28);
						FMLA_ZSP( z9,p0,z1,z29);
						FMLA_ZSP( z10,p0,z1,z30);
						FMLA_ZSP( z11,p0,z1,z31);

						LD1W_ZXI(z1,p1,aa,3);
						aa+=nvl*3+mm;

						FMLA_ZSP( z12,p0,z0,z27);
						FMLA_ZSP( z13,p0,z0,z28);

						bq4+=2;
						ik+=2;

						FMLA_ZSP( z14,p0,z0,z29);
						FMLA_ZSP( z15,p0,z0,z30);
						FMLA_ZSP( z16,p0,z0,z31);
						FMLA_ZSP( z17,p1,z1,z27);

						LD1W_ZXI(z0,p0,aa,0);
						NOP;

						FMLA_ZSP( z18,p1,z1,z28);
						FMLA_ZSP( z19,p1,z1,z29);

						LD1RW_ZXI(z27,p0,bq0,SIZE);
						LD1RW_ZXI(z28,p0,bq1,SIZE);

						FMLA_ZSP( z20,p1,z1,z30);
						FMLA_ZSP( z21,p1,z1,z31);
					}


				}

				LD1W_ZXI(z1,p0,aa,1);

				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z3,p0,z0,z23);

				LD1RW_ZXI(z29,p0,bq2,SIZE);
				LD1RW_ZXI(z30,p0,bq3,SIZE);

				FMLA_ZSP( z4,p0,z0,z24);
				FMLA_ZSP( z5,p0,z0,z25);
				FMLA_ZSP( z6,p0,z0,z26);

				LD1W_ZXI(z0,p0,aa,2);
				LD1RW_ZXI(z31,p0,bq4,SIZE);

				FMLA_ZSP( z7,p0,z1,z22);
				FMLA_ZSP( z8,p0,z1,z23);
				FMLA_ZSP( z9,p0,z1,z24);
				FMLA_ZSP( z10,p0,z1,z25);
				FMLA_ZSP( z11,p0,z1,z26);

				LD1W_ZXI(z1,p1,aa,3);
				aa+=nvl*3+mm;

				FMLA_ZSP( z12,p0,z0,z22);
				FMLA_ZSP( z13,p0,z0,z23);
				FMLA_ZSP( z14,p0,z0,z24);
				FMLA_ZSP( z15,p0,z0,z25);
				FMLA_ZSP( z16,p0,z0,z26);

				LD1W_ZXI(z0,p0,aa,0);

				FMLA_ZSP( z17,p1,z1,z22);
				FMLA_ZSP( z18,p1,z1,z23);
				FMLA_ZSP( z19,p1,z1,z24);
				FMLA_ZSP( z20,p1,z1,z25);
				FMLA_ZSP( z21,p1,z1,z26);

				LD1W_ZXI(z1,p0,aa,1);

				FMLA_ZSP( z2,p0,z0,z27);
				FMLA_ZSP( z3,p0,z0,z28);
				FMLA_ZSP( z4,p0,z0,z29);
				FMLA_ZSP( z5,p0,z0,z30);
				FMLA_ZSP( z6,p0,z0,z31);

				LD1W_ZXI(z0,p0,aa,2);

				FMLA_ZSP( z7,p0,z1,z27);
				FMLA_ZSP( z8,p0,z1,z28);
				FMLA_ZSP( z9,p0,z1,z29);
				FMLA_ZSP( z10,p0,z1,z30);
				FMLA_ZSP( z11,p0,z1,z31);

				LD1W_ZXI(z1,p1,aa,3);

				FMLA_ZSP( z12,p0,z0,z27);
				FMLA_ZSP( z13,p0,z0,z28);
				FMLA_ZSP( z14,p0,z0,z29);
				FMLA_ZSP( z15,p0,z0,z30);
				FMLA_ZSP( z16,p0,z0,z31);
				FMLA_ZSP( z17,p1,z1,z27);
				FMLA_ZSP( z18,p1,z1,z28);
				FMLA_ZSP( z19,p1,z1,z29);
				FMLA_ZSP( z20,p1,z1,z30);
				FMLA_ZSP( z21,p1,z1,z31);

				bq0+=2;
				bq1+=2;
				bq2+=2;
				bq3+=2;
				bq4+=2;
				aa+=nvl*3+mm;

			}
			if(k-ik>=1) {
				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);
				LD1RW_ZXI(z24,p0,bq2,0);
				LD1RW_ZXI(z25,p0,bq3,0);
				LD1RW_ZXI(z26,p0,bq4,0);

				LD1W_ZXI(z0,p0,aa,0);
				LD1W_ZXI(z1,p0,aa,1);

				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z3,p0,z0,z23);
				FMLA_ZSP( z4,p0,z0,z24);
				FMLA_ZSP( z5,p0,z0,z25);
				FMLA_ZSP( z6,p0,z0,z26);
				FMLA_ZSP( z7,p0,z1,z22);
				FMLA_ZSP( z8,p0,z1,z23);
				FMLA_ZSP( z9,p0,z1,z24);
				FMLA_ZSP( z10,p0,z1,z25);
				FMLA_ZSP( z11,p0,z1,z26);

				LD1W_ZXI(z0,p0,aa,2);
				LD1W_ZXI(z1,p1,aa,3);

				FMLA_ZSP( z12,p0,z0,z22);
				FMLA_ZSP( z13,p0,z0,z23);
				FMLA_ZSP( z14,p0,z0,z24);
				FMLA_ZSP( z15,p0,z0,z25);
				FMLA_ZSP( z16,p0,z0,z26);
				FMLA_ZSP( z17,p1,z1,z22);
				FMLA_ZSP( z18,p1,z1,z23);
				FMLA_ZSP( z19,p1,z1,z24);
				FMLA_ZSP( z20,p1,z1,z25);
				FMLA_ZSP( z21,p1,z1,z26);

				aa+=nvl*3+mm;
				bq0+=1;
				bq1+=1;
				bq2+=1;
				bq3+=1;
				bq4+=1;

			}
			LD1RW_ZXI(z0,p0,alphap,0);
			CALPHA_P_45;
			cq0+=nvl*3+mm;
			cq1+=nvl*3+mm;
			cq2+=nvl*3+mm;
			cq3+=nvl*3+mm;
			cq4+=nvl*3+mm;
		}

		else if(m-im>nvl*2) {
			mm=m-im;
			im+=nvl*2;
			WHILELT_PSX(p1,im,m);
			mm=m-im;
			im=m;
			PRFM_XI(PSTL1KEEP,cq0,256*2);
			PRFM_XI(PSTL1KEEP,cq1,256*2);
			PRFM_XI(PSTL1KEEP,cq2,256*2);
			PRFM_XI(PSTL1KEEP,cq3,256*2);
			PRFM_XI(PSTL1KEEP,cq4,256*2);

			CZERO_35;

			PRFM_XI(PLDL1KEEP,aa,0);
			PRFM_XI(PLDL1KEEP,aa,256*1);
			PRFM_XI(PLDL1KEEP,aa,256*2);
			PRFM_XI(PLDL1KEEP,aa,256*3);
			PRFM_XI(PLDL1KEEP,aa,256*4);
			PRFM_XI(PLDL1KEEP,aa,256*5);
			PRFM_XI(PLDL1KEEP,aa,256*6);
			PRFM_XI(PLDL1KEEP,aa,256*7);
			PRFM_XI(PLDL1KEEP,aa,256*8);

			bq0=bs0;
			bq1=bs1;
			bq2=bs2;
			bq3=bs3;
			bq4=bs4;
			ik=0;
			if(k>=2) { 

				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);
				LD1RW_ZXI(z24,p0,bq2,0);
				LD1RW_ZXI(z25,p0,bq3,0);
				LD1RW_ZXI(z26,p0,bq4,0);

				LD1RW_ZXI(z27,p0,bq0,SIZE);
				LD1RW_ZXI(z28,p0,bq1,SIZE);

				LD1W_ZXI(z0,p0,aa,0);

				ik+=2;

				for(;ik<k-1;) {
                                        PRFM_XI(PLDL1KEEP,bq0,256);
                                        PRFM_XI(PLDL1KEEP,bq1,256);
                                        PRFM_XI(PLDL1KEEP,bq2,256);
                                        PRFM_XI(PLDL1KEEP,bq3,256);
                                        PRFM_XI(PLDL1KEEP,bq4,256);

					LD1W_ZXI(z1,p0,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p0,z0,z22);
					FMLA_ZSP( z3,p0,z0,z23);

					LD1RW_ZXI(z29,p0,bq2,SIZE);
					LD1RW_ZXI(z30,p0,bq3,SIZE);

					FMLA_ZSP( z4,p0,z0,z24);
					FMLA_ZSP( z5,p0,z0,z25);
					FMLA_ZSP( z6,p0,z0,z26);
					FMLA_ZSP( z7,p0,z1,z22);

					LD1W_ZXI(z17,p1,aa,2);
					LD1RW_ZXI(z31,p0,bq4,SIZE);

					FMLA_ZSP( z8,p0,z1,z23);
					FMLA_ZSP( z9,p0,z1,z24);

					FMLA_ZSP( z10,p0,z1,z25);
					FMLA_ZSP( z11,p0,z1,z26);

					aa+=nvl*2+mm;
					bq0+=2;

					FMLA_ZSP( z12,p1,z17,z22);
					FMLA_ZSP( z13,p1,z17,z23);

					bq1+=2;
					bq2+=2;

					FMLA_ZSP( z14,p1,z17,z24);
					FMLA_ZSP( z15,p1,z17,z25);
					FMLA_ZSP( z16,p1,z17,z26);

					LD1W_ZXI(z0,p0,aa,0);
					bq3+=2;

					LD1RW_ZXI(z22,p0,bq0,0);
					LD1RW_ZXI(z23,p0,bq1,0);

					LD1W_ZXI(z1,p0,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p0,z0,z27);
					FMLA_ZSP( z3,p0,z0,z28);

					LD1RW_ZXI(z24,p0,bq2,0);
					LD1RW_ZXI(z25,p0,bq3,0);

					FMLA_ZSP( z4,p0,z0,z29);
					FMLA_ZSP( z5,p0,z0,z30);
					FMLA_ZSP( z6,p0,z0,z31);
					FMLA_ZSP( z7,p0,z1,z27);

					LD1W_ZXI(z17,p1,aa,2);
					LD1RW_ZXI(z26,p0,bq4,SIZE*2);

					FMLA_ZSP( z8,p0,z1,z28);
					FMLA_ZSP( z9,p0,z1,z29);
					FMLA_ZSP( z10,p0,z1,z30);
					FMLA_ZSP( z11,p0,z1,z31);

					aa+=nvl*2+mm;

					FMLA_ZSP( z12,p1,z17,z27);
					FMLA_ZSP( z13,p1,z17,z28);

					bq4+=2;
					ik+=2;

					FMLA_ZSP( z14,p1,z17,z29);
					FMLA_ZSP( z15,p1,z17,z30);
					FMLA_ZSP( z16,p1,z17,z31);

					LD1W_ZXI(z0,p0,aa,0);
					NOP;

					LD1RW_ZXI(z27,p0,bq0,SIZE);
					LD1RW_ZXI(z28,p0,bq1,SIZE);

				}
				LD1W_ZXI(z1,p0,aa,1);

				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z3,p0,z0,z23);

				LD1RW_ZXI(z29,p0,bq2,SIZE);
				LD1RW_ZXI(z30,p0,bq3,SIZE);

				FMLA_ZSP( z4,p0,z0,z24);
				FMLA_ZSP( z5,p0,z0,z25);
				FMLA_ZSP( z6,p0,z0,z26);

				LD1W_ZXI(z17,p1,aa,2);
				LD1RW_ZXI(z31,p0,bq4,SIZE);

				FMLA_ZSP( z7,p0,z1,z22);
				FMLA_ZSP( z8,p0,z1,z23);
				FMLA_ZSP( z9,p0,z1,z24);
				FMLA_ZSP( z10,p0,z1,z25);
				FMLA_ZSP( z11,p0,z1,z26);

				aa+=nvl*2+mm;

				FMLA_ZSP( z12,p1,z17,z22);
				FMLA_ZSP( z13,p1,z17,z23);
				FMLA_ZSP( z14,p1,z17,z24);
				FMLA_ZSP( z15,p1,z17,z25);
				FMLA_ZSP( z16,p1,z17,z26);

				LD1W_ZXI(z0,p0,aa,0);

				LD1W_ZXI(z1,p0,aa,1);

				FMLA_ZSP( z2,p0,z0,z27);
				FMLA_ZSP( z3,p0,z0,z28);
				FMLA_ZSP( z4,p0,z0,z29);
				FMLA_ZSP( z5,p0,z0,z30);
				FMLA_ZSP( z6,p0,z0,z31);

				LD1W_ZXI(z17,p1,aa,2);

				FMLA_ZSP( z7,p0,z1,z27);
				FMLA_ZSP( z8,p0,z1,z28);
				FMLA_ZSP( z9,p0,z1,z29);
				FMLA_ZSP( z10,p0,z1,z30);
				FMLA_ZSP( z11,p0,z1,z31);

				FMLA_ZSP( z12,p1,z17,z27);
				FMLA_ZSP( z13,p1,z17,z28);
				FMLA_ZSP( z14,p1,z17,z29);
				FMLA_ZSP( z15,p1,z17,z30);
				FMLA_ZSP( z16,p1,z17,z31);

				bq0+=2;
				bq1+=2;
				bq2+=2;
				bq3+=2;
				bq4+=2;
				aa+=nvl*2+mm;

			}
			if(k-ik>=1) {
				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);
				LD1RW_ZXI(z24,p0,bq2,0);
				LD1RW_ZXI(z25,p0,bq3,0);
				LD1RW_ZXI(z26,p0,bq4,0);

				LD1W_ZXI(z0,p0,aa,0);
				LD1W_ZXI(z1,p0,aa,1);
				LD1W_ZXI(z17,p1,aa,2);

				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z3,p0,z0,z23);
				FMLA_ZSP( z4,p0,z0,z24);
				FMLA_ZSP( z5,p0,z0,z25);
				FMLA_ZSP( z6,p0,z0,z26);

				FMLA_ZSP( z7,p0,z1,z22);
				FMLA_ZSP( z8,p0,z1,z23);
				FMLA_ZSP( z9,p0,z1,z24);
				FMLA_ZSP( z10,p0,z1,z25);
				FMLA_ZSP( z11,p0,z1,z26);

				FMLA_ZSP( z12,p1,z17,z22);
				FMLA_ZSP( z13,p1,z17,z23);
				FMLA_ZSP( z14,p1,z17,z24);
				FMLA_ZSP( z15,p1,z17,z25);
				FMLA_ZSP( z16,p1,z17,z26);

				aa+=nvl*2+mm;
				bq0+=1;
				bq1+=1;
				bq2+=1;
				bq3+=1;
				bq4+=1;
			}
			LD1RW_ZXI(z0,p0,alphap,0);
			CALPHA_P_35;
			cq0+=nvl*2+mm;
			cq1+=nvl*2+mm;
			cq2+=nvl*2+mm;
			cq3+=nvl*2+mm;
			cq4+=nvl*2+mm;
		}

		else if(m-im>nvl*1) {
			mm=m-im;
			im+=nvl*1;
			WHILELT_PSX(p1,im,m);
			mm=m-im;
			im=m;
			PRFM_XI(PSTL1KEEP,cq0,256*2);
			PRFM_XI(PSTL1KEEP,cq1,256*2);
			PRFM_XI(PSTL1KEEP,cq2,256*2);
			PRFM_XI(PSTL1KEEP,cq3,256*2);
			PRFM_XI(PSTL1KEEP,cq4,256*2);

			CZERO_25;

			PRFM_XI(PLDL1KEEP,aa,0);
			PRFM_XI(PLDL1KEEP,aa,256*1);
			PRFM_XI(PLDL1KEEP,aa,256*2);
			PRFM_XI(PLDL1KEEP,aa,256*3);
			PRFM_XI(PLDL1KEEP,aa,256*4);
			PRFM_XI(PLDL1KEEP,aa,256*5);
			PRFM_XI(PLDL1KEEP,aa,256*6);
			PRFM_XI(PLDL1KEEP,aa,256*7);
			PRFM_XI(PLDL1KEEP,aa,256*8);

			bq0=bs0;
			bq1=bs1;
			bq2=bs2;
			bq3=bs3;
			bq4=bs4;
			ik=0;
			if(k>=2) { 

				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);
				LD1RW_ZXI(z24,p0,bq2,0);
				LD1RW_ZXI(z25,p0,bq3,0);
				LD1RW_ZXI(z26,p0,bq4,0);

				LD1RW_ZXI(z27,p0,bq0,SIZE);
				LD1RW_ZXI(z28,p0,bq1,SIZE);

				LD1W_ZXI(z0,p0,aa,0);

				ik+=2;

				for(;ik<k-1;) {
                                        PRFM_XI(PLDL1KEEP,bq0,256);
                                        PRFM_XI(PLDL1KEEP,bq1,256);
                                        PRFM_XI(PLDL1KEEP,bq2,256);
                                        PRFM_XI(PLDL1KEEP,bq3,256);
                                        PRFM_XI(PLDL1KEEP,bq4,256);

					LD1W_ZXI(z1,p1,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p0,z0,z22);
					FMLA_ZSP( z3,p0,z0,z23);

					LD1RW_ZXI(z29,p0,bq2,SIZE);
					LD1RW_ZXI(z30,p0,bq3,SIZE);

					FMLA_ZSP( z4,p0,z0,z24);
					FMLA_ZSP( z5,p0,z0,z25);
					FMLA_ZSP( z6,p0,z0,z26);
					FMLA_ZSP( z7,p1,z1,z22);

					LD1RW_ZXI(z31,p0,bq4,SIZE);

					FMLA_ZSP( z8,p1,z1,z23);
					FMLA_ZSP( z9,p1,z1,z24);

					FMLA_ZSP( z10,p1,z1,z25);
					FMLA_ZSP( z11,p1,z1,z26);

					aa+=nvl*1+mm;
					bq0+=2;
					bq1+=2;
					bq2+=2;

					LD1W_ZXI(z0,p0,aa,0);
					bq3+=2;

					LD1RW_ZXI(z22,p0,bq0,0);
					LD1RW_ZXI(z23,p0,bq1,0);

					LD1W_ZXI(z1,p1,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p0,z0,z27);
					FMLA_ZSP( z3,p0,z0,z28);

					LD1RW_ZXI(z24,p0,bq2,0);
					LD1RW_ZXI(z25,p0,bq3,0);

					FMLA_ZSP( z4,p0,z0,z29);
					FMLA_ZSP( z5,p0,z0,z30);
					FMLA_ZSP( z6,p0,z0,z31);
					FMLA_ZSP( z7,p1,z1,z27);

					LD1RW_ZXI(z26,p0,bq4,SIZE*2);

					FMLA_ZSP( z8,p1,z1,z28);
					FMLA_ZSP( z9,p1,z1,z29);
					FMLA_ZSP( z10,p1,z1,z30);
					FMLA_ZSP( z11,p1,z1,z31);

					aa+=nvl*1+mm;

					bq4+=2;
					ik+=2;

					LD1W_ZXI(z0,p0,aa,0);
					NOP;

					LD1RW_ZXI(z27,p0,bq0,SIZE);
					LD1RW_ZXI(z28,p0,bq1,SIZE);

				}
				LD1W_ZXI(z1,p1,aa,1);

				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z3,p0,z0,z23);

				LD1RW_ZXI(z29,p0,bq2,SIZE);
				LD1RW_ZXI(z30,p0,bq3,SIZE);

				FMLA_ZSP( z4,p0,z0,z24);
				FMLA_ZSP( z5,p0,z0,z25);
				FMLA_ZSP( z6,p0,z0,z26);

				LD1RW_ZXI(z31,p0,bq4,SIZE);

				FMLA_ZSP( z7,p1,z1,z22);
				FMLA_ZSP( z8,p1,z1,z23);
				FMLA_ZSP( z9,p1,z1,z24);
				FMLA_ZSP( z10,p1,z1,z25);
				FMLA_ZSP( z11,p1,z1,z26);

				aa+=nvl*1+mm;

				LD1W_ZXI(z0,p0,aa,0);

				LD1W_ZXI(z1,p1,aa,1);

				FMLA_ZSP( z2,p0,z0,z27);
				FMLA_ZSP( z3,p0,z0,z28);
				FMLA_ZSP( z4,p0,z0,z29);
				FMLA_ZSP( z5,p0,z0,z30);
				FMLA_ZSP( z6,p0,z0,z31);

				FMLA_ZSP( z7,p1,z1,z27);
				FMLA_ZSP( z8,p1,z1,z28);
				FMLA_ZSP( z9,p1,z1,z29);
				FMLA_ZSP( z10,p1,z1,z30);
				FMLA_ZSP( z11,p1,z1,z31);

				bq0+=2;
				bq1+=2;
				bq2+=2;
				bq3+=2;
				bq4+=2;
				aa+=nvl*1+mm;

			}

			if(k-ik>=1) {
				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);
				LD1RW_ZXI(z24,p0,bq2,0);
				LD1RW_ZXI(z25,p0,bq3,0);
				LD1RW_ZXI(z26,p0,bq4,0);

				LD1W_ZXI(z0,p0,aa,0);
				LD1W_ZXI(z1,p1,aa,1);

				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z3,p0,z0,z23);
				FMLA_ZSP( z4,p0,z0,z24);
				FMLA_ZSP( z5,p0,z0,z25);
				FMLA_ZSP( z6,p0,z0,z26);

				FMLA_ZSP( z7,p1,z1,z22);
				FMLA_ZSP( z8,p1,z1,z23);
				FMLA_ZSP( z9,p1,z1,z24);
				FMLA_ZSP( z10,p1,z1,z25);
				FMLA_ZSP( z11,p1,z1,z26);

				aa+=nvl*1+mm;
				bq0+=1;
				bq1+=1;
				bq2+=1;
				bq3+=1;
				bq4+=1;

			}
			LD1RW_ZXI(z0,p0,alphap,0);
			CALPHA_P_25;
			cq0+=nvl*1+mm;
			cq1+=nvl*1+mm;
			cq2+=nvl*1+mm;
			cq3+=nvl*1+mm;
			cq4+=nvl*1+mm;
		}


		else if(m-im>0) {
			mm=m-im;
			WHILELT_PSX(p1,im,m);
			mm=m-im;
			im=m;
			PRFM_XI(PSTL1KEEP,cq0,256*2);
			PRFM_XI(PSTL1KEEP,cq1,256*2);
			PRFM_XI(PSTL1KEEP,cq2,256*2);
			PRFM_XI(PSTL1KEEP,cq3,256*2);
			PRFM_XI(PSTL1KEEP,cq4,256*2);

			CZERO_15;

			PRFM_XI(PLDL1KEEP,aa,0);
			PRFM_XI(PLDL1KEEP,aa,256*1);
			PRFM_XI(PLDL1KEEP,aa,256*2);
			PRFM_XI(PLDL1KEEP,aa,256*3);
			PRFM_XI(PLDL1KEEP,aa,256*4);
			PRFM_XI(PLDL1KEEP,aa,256*5);
			PRFM_XI(PLDL1KEEP,aa,256*6);
			PRFM_XI(PLDL1KEEP,aa,256*7);
			PRFM_XI(PLDL1KEEP,aa,256*8);

			bq0=bs0;
			bq1=bs1;
			bq2=bs2;
			bq3=bs3;
			bq4=bs4;
			ik=0;
			if(k>=2) { 

				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);
				LD1RW_ZXI(z24,p0,bq2,0);
				LD1RW_ZXI(z25,p0,bq3,0);
				LD1RW_ZXI(z26,p0,bq4,0);

				LD1RW_ZXI(z27,p0,bq0,SIZE);
				LD1RW_ZXI(z28,p0,bq1,SIZE);

				LD1W_ZXI(z0,p1,aa,0);

				ik+=2;

				for(;ik<k-1;) {
                                        PRFM_XI(PLDL1KEEP,bq0,256);
                                        PRFM_XI(PLDL1KEEP,bq1,256);
                                        PRFM_XI(PLDL1KEEP,bq2,256);
                                        PRFM_XI(PLDL1KEEP,bq3,256);
                                        PRFM_XI(PLDL1KEEP,bq4,256);

					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p1,z0,z22);
					FMLA_ZSP( z3,p1,z0,z23);

					LD1RW_ZXI(z29,p0,bq2,SIZE);
					LD1RW_ZXI(z30,p0,bq3,SIZE);

					FMLA_ZSP( z4,p1,z0,z24);
					FMLA_ZSP( z5,p1,z0,z25);
					FMLA_ZSP( z6,p1,z0,z26);

					LD1RW_ZXI(z31,p0,bq4,SIZE);

					aa+=mm;
					bq0+=2;
					bq1+=2;
					bq2+=2;

					LD1W_ZXI(z0,p1,aa,0);
					bq3+=2;

					LD1RW_ZXI(z22,p0,bq0,0);
					LD1RW_ZXI(z23,p0,bq1,0);

					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p1,z0,z27);
					FMLA_ZSP( z3,p1,z0,z28);

					LD1RW_ZXI(z24,p0,bq2,0);
					LD1RW_ZXI(z25,p0,bq3,0);

					FMLA_ZSP( z4,p1,z0,z29);
					FMLA_ZSP( z5,p1,z0,z30);
					FMLA_ZSP( z6,p1,z0,z31);

					LD1RW_ZXI(z26,p0,bq4,SIZE*2);

					aa+=mm;

					bq4+=2;
					ik+=2;

					LD1W_ZXI(z0,p1,aa,0);
					NOP;

					LD1RW_ZXI(z27,p0,bq0,SIZE);
					LD1RW_ZXI(z28,p0,bq1,SIZE);

				}

				FMLA_ZSP( z2,p1,z0,z22);
				FMLA_ZSP( z3,p1,z0,z23);

				LD1RW_ZXI(z29,p0,bq2,SIZE);
				LD1RW_ZXI(z30,p0,bq3,SIZE);

				FMLA_ZSP( z4,p1,z0,z24);
				FMLA_ZSP( z5,p1,z0,z25);
				FMLA_ZSP( z6,p1,z0,z26);

				LD1RW_ZXI(z31,p0,bq4,SIZE);

				aa+=mm;

				LD1W_ZXI(z0,p1,aa,0);

				FMLA_ZSP( z2,p1,z0,z27);
				FMLA_ZSP( z3,p1,z0,z28);
				FMLA_ZSP( z4,p1,z0,z29);
				FMLA_ZSP( z5,p1,z0,z30);
				FMLA_ZSP( z6,p1,z0,z31);

				bq0+=2;
				bq1+=2;
				bq2+=2;
				bq3+=2;
				bq4+=2;
				aa+=mm;

			}
			if(k-ik>=1) {
				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);
				LD1RW_ZXI(z24,p0,bq2,0);
				LD1RW_ZXI(z25,p0,bq3,0);
				LD1RW_ZXI(z26,p0,bq4,0);

				LD1W_ZXI(z0,p1,aa,0);

				FMLA_ZSP( z2,p1,z0,z22);
				FMLA_ZSP( z3,p1,z0,z23);
				FMLA_ZSP( z4,p1,z0,z24);
				FMLA_ZSP( z5,p1,z0,z25);
				FMLA_ZSP( z6,p1,z0,z26);

				aa+=mm;
				bq0+=1;
				bq1+=1;
				bq2+=1;
				bq3+=1;
				bq4+=1;

			}
			LD1RW_ZXI(z0,p0,alphap,0);
			CALPHA_P_15;
			cq0+=mm;
			cq1+=mm;
			cq2+=mm;
			cq3+=mm;
			cq4+=mm;
		}

		cq0+=ldc*5-m;
		cq1+=ldc*5-m;
		cq2+=ldc*5-m;
		cq3+=ldc*5-m;
		cq4+=ldc*5-m;
		bs0+=ldb*5;
		bs1+=ldb*5;
		bs2+=ldb*5;
		bs3+=ldb*5;
		bs4+=ldb*5;
	}

	if(n-in==4) {
		in+=4;
		aa=a;
		im=0;
		for(;im<m-nul+1;im+=nul) {
			PRFM_XI(PSTL1KEEP,cq0,256*2);
			PRFM_XI(PSTL1KEEP,cq1,256*2);
			PRFM_XI(PSTL1KEEP,cq2,256*2);
			PRFM_XI(PSTL1KEEP,cq3,256*2);
			CZERO_44;

			PRFM_XI(PLDL1KEEP,aa,0);
			PRFM_XI(PLDL1KEEP,aa,256*1);
			PRFM_XI(PLDL1KEEP,aa,256*2);
			PRFM_XI(PLDL1KEEP,aa,256*3);
			PRFM_XI(PLDL1KEEP,aa,256*4);
			PRFM_XI(PLDL1KEEP,aa,256*5);
			PRFM_XI(PLDL1KEEP,aa,256*6);
			PRFM_XI(PLDL1KEEP,aa,256*7);
			PRFM_XI(PLDL1KEEP,aa,256*8);


			bq0=bs0;
			bq1=bs1;
			bq2=bs2;
			bq3=bs3;
			ik=0;
			if(k>=2) { 

				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);
				LD1RW_ZXI(z24,p0,bq2,0);
				LD1RW_ZXI(z25,p0,bq3,0);

				LD1RW_ZXI(z27,p0,bq0,SIZE);
				LD1RW_ZXI(z28,p0,bq1,SIZE);

				LD1W_ZXI(z0,p0,aa,0);

				ik+=2;

				for(;ik<k-1;) {
                                        PRFM_XI(PLDL1KEEP,bq0,256);
                                        PRFM_XI(PLDL1KEEP,bq1,256);
                                        PRFM_XI(PLDL1KEEP,bq2,256);
                                        PRFM_XI(PLDL1KEEP,bq3,256);

					LD1W_ZXI(z1,p0,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p0,z0,z22);

					FMLA_ZSP( z3,p0,z0,z23);

					LD1RW_ZXI(z29,p0,bq2,SIZE);
					LD1RW_ZXI(z30,p0,bq3,SIZE);

					FMLA_ZSP( z4,p0,z0,z24);
					FMLA_ZSP( z5,p0,z0,z25);
					FMLA_ZSP( z7,p0,z1,z22);

					LD1W_ZXI(z0,p0,aa,2);

					FMLA_ZSP( z8,p0,z1,z23);
					FMLA_ZSP( z9,p0,z1,z24);
					FMLA_ZSP( z10,p0,z1,z25);

					LD1W_ZXI(z1,p0,aa,3);
					bq0+=2;

					FMLA_ZSP( z12,p0,z0,z22);
					FMLA_ZSP( z13,p0,z0,z23);

					bq1+=2;
					bq2+=2;

					FMLA_ZSP( z14,p0,z0,z24);
					FMLA_ZSP( z15,p0,z0,z25);
					FMLA_ZSP( z17,p0,z1,z22);

					LD1W_ZXI(z0,p0,aa,4);
					bq3+=2;

					FMLA_ZSP( z18,p0,z1,z23);
					FMLA_ZSP( z19,p0,z1,z24);

					LD1RW_ZXI(z22,p0,bq0,0);
					LD1RW_ZXI(z23,p0,bq1,0);

					FMLA_ZSP( z20,p0,z1,z25);

					LD1W_ZXI(z1,p0,aa,5);
					PRFM_XI(PLDL1KEEP,aa,256*10);

					FMLA_ZSP( z2,p0,z0,z27);
					FMLA_ZSP( z3,p0,z0,z28);

					LD1RW_ZXI(z24,p0,bq2,0);
					LD1RW_ZXI(z25,p0,bq3,0);

					FMLA_ZSP( z4,p0,z0,z29);
					FMLA_ZSP( z5,p0,z0,z30);
					FMLA_ZSP( z7,p0,z1,z27);

					LD1W_ZXI(z0,p0,aa,6);

					FMLA_ZSP( z8,p0,z1,z28);
					FMLA_ZSP( z9,p0,z1,z29);
					FMLA_ZSP( z10,p0,z1,z30);

					LD1W_ZXI(z1,p0,aa,7);
					aa+=nul*2;

					FMLA_ZSP( z12,p0,z0,z27);
					FMLA_ZSP( z13,p0,z0,z28);

					ik+=2;

					FMLA_ZSP( z14,p0,z0,z29);
					FMLA_ZSP( z15,p0,z0,z30);
					FMLA_ZSP( z17,p0,z1,z27);

					LD1W_ZXI(z0,p0,aa,0);
					NOP;

					FMLA_ZSP( z18,p0,z1,z28);
					FMLA_ZSP( z19,p0,z1,z29);

					LD1RW_ZXI(z27,p0,bq0,SIZE);
					LD1RW_ZXI(z28,p0,bq1,SIZE);

					FMLA_ZSP( z20,p0,z1,z30);
				}
				LD1W_ZXI(z1,p0,aa,1);

				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z3,p0,z0,z23);

				LD1RW_ZXI(z29,p0,bq2,SIZE);
				LD1RW_ZXI(z30,p0,bq3,SIZE);

				FMLA_ZSP( z4,p0,z0,z24);
				FMLA_ZSP( z5,p0,z0,z25);

				LD1W_ZXI(z0,p0,aa,2);

				FMLA_ZSP( z7 ,p0,z1,z22);
				FMLA_ZSP( z8 ,p0,z1,z23);
				FMLA_ZSP( z9 ,p0,z1,z24);
				FMLA_ZSP( z10,p0,z1,z25);

				LD1W_ZXI(z1,p0,aa,3);

				FMLA_ZSP( z12,p0,z0,z22);
				FMLA_ZSP( z13,p0,z0,z23);
				FMLA_ZSP( z14,p0,z0,z24);
				FMLA_ZSP( z15,p0,z0,z25);

				LD1W_ZXI(z0,p0,aa,4);

				FMLA_ZSP( z17,p0,z1,z22);
				FMLA_ZSP( z18,p0,z1,z23);
				FMLA_ZSP( z19,p0,z1,z24);
				FMLA_ZSP( z20,p0,z1,z25);

				LD1W_ZXI(z1,p0,aa,5);

				FMLA_ZSP( z2,p0,z0,z27);
				FMLA_ZSP( z3,p0,z0,z28);
				FMLA_ZSP( z4,p0,z0,z29);
				FMLA_ZSP( z5,p0,z0,z30);


				LD1W_ZXI(z0,p0,aa,6);

				FMLA_ZSP( z7,p0,z1,z27);
				FMLA_ZSP( z8,p0,z1,z28);
				FMLA_ZSP( z9,p0,z1,z29);
				FMLA_ZSP( z10,p0,z1,z30);

				LD1W_ZXI(z1,p0,aa,7);
				bq0+=2;
				bq1+=2;
				bq2+=2;
				bq3+=2;
				aa+=nul*2;

				if(k-ik>=1) {
					LD1RW_ZXI(z22,p0,bq0,0);
					LD1RW_ZXI(z23,p0,bq1,0);

					FMLA_ZSP( z12,p0,z0,z27);
					FMLA_ZSP( z13,p0,z0,z28);

					LD1RW_ZXI(z24,p0,bq2,0);
					LD1RW_ZXI(z25,p0,bq3,0);

					FMLA_ZSP( z14,p0,z0,z29);
					FMLA_ZSP( z15,p0,z0,z30);

					LD1W_ZXI(z0,p0,aa,0);

					FMLA_ZSP( z17,p0,z1,z27);
					FMLA_ZSP( z18,p0,z1,z28);
					FMLA_ZSP( z19,p0,z1,z29);
					FMLA_ZSP( z20,p0,z1,z30);


					LD1W_ZXI(z1,p0,aa,1);
				}
				else {
					FMLA_ZSP( z12,p0,z0,z27);
					FMLA_ZSP( z13,p0,z0,z28);
					FMLA_ZSP( z14,p0,z0,z29);
					FMLA_ZSP( z15,p0,z0,z30);
					FMLA_ZSP( z17,p0,z1,z27);
					FMLA_ZSP( z18,p0,z1,z28);
					FMLA_ZSP( z19,p0,z1,z29);
					FMLA_ZSP( z20,p0,z1,z30);
				}
			}
			else if(k==1) {
				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);
				LD1RW_ZXI(z24,p0,bq2,0);
				LD1RW_ZXI(z25,p0,bq3,0);

				LD1W_ZXI(z0,p0,aa,0);
				LD1W_ZXI(z1,p0,aa,1);

			}

			if(k-ik>=1) {

				FMLA_ZSP( z2, p0,z0,z22);
				FMLA_ZSP( z3, p0,z0,z23);
				FMLA_ZSP( z4, p0,z0,z24);
				FMLA_ZSP( z5, p0,z0,z25);

				LD1W_ZXI(z0,p0,aa,2);

				FMLA_ZSP( z7, p0,z1,z22);
				FMLA_ZSP( z8, p0,z1,z23);
				FMLA_ZSP( z9, p0,z1,z24);
				FMLA_ZSP( z10,p0,z1,z25);

				LD1W_ZXI(z1,p0,aa,3);

				FMLA_ZSP( z12,p0,z0,z22);
				FMLA_ZSP( z13,p0,z0,z23);
				FMLA_ZSP( z14,p0,z0,z24);
				FMLA_ZSP( z15,p0,z0,z25);
				FMLA_ZSP( z17,p0,z1,z22);
				FMLA_ZSP( z18,p0,z1,z23);
				FMLA_ZSP( z19,p0,z1,z24);
				FMLA_ZSP( z20,p0,z1,z25);

				aa+=nul;
				bq0+=1;
				bq1+=1;
				bq2+=1;
				bq3+=1;
				ik++;
			}
			LD1RW_ZXI(z0,p0,alphap,0);
			CALPHA_44;
			cq0+=nul;
			cq1+=nul;
			cq2+=nul;
			cq3+=nul;
		}

		if(m-im>nvl*3) {
			mm=m-im;
			im+=nvl*3;
			WHILELT_PSX(p1,im,m);
			mm=m-im;
			im=m;
			PRFM_XI(PSTL1KEEP,cq0,256*2);
			PRFM_XI(PSTL1KEEP,cq1,256*2);
			PRFM_XI(PSTL1KEEP,cq2,256*2);
			PRFM_XI(PSTL1KEEP,cq3,256*2);

			CZERO_44;

			PRFM_XI(PLDL1KEEP,aa,0);
			PRFM_XI(PLDL1KEEP,aa,256*1);
			PRFM_XI(PLDL1KEEP,aa,256*2);
			PRFM_XI(PLDL1KEEP,aa,256*3);
			PRFM_XI(PLDL1KEEP,aa,256*4);
			PRFM_XI(PLDL1KEEP,aa,256*5);
			PRFM_XI(PLDL1KEEP,aa,256*6);
			PRFM_XI(PLDL1KEEP,aa,256*7);
			PRFM_XI(PLDL1KEEP,aa,256*8);

			bq0=bs0;
			bq1=bs1;
			bq2=bs2;
			bq3=bs3;
			ik=0;
			if(k>=2) { 

				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);
				LD1RW_ZXI(z24,p0,bq2,0);
				LD1RW_ZXI(z25,p0,bq3,0);

				LD1RW_ZXI(z27,p0,bq0,SIZE);
				LD1RW_ZXI(z28,p0,bq1,SIZE);

				LD1W_ZXI(z0,p0,aa,0);

				ik+=2;

				for(;ik<k-1;) {
                                        PRFM_XI(PLDL1KEEP,bq0,256);
                                        PRFM_XI(PLDL1KEEP,bq1,256);
                                        PRFM_XI(PLDL1KEEP,bq2,256);
                                        PRFM_XI(PLDL1KEEP,bq3,256);

					LD1W_ZXI(z1,p0,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p0,z0,z22);
					FMLA_ZSP( z3,p0,z0,z23);

					LD1RW_ZXI(z29,p0,bq2,SIZE);
					LD1RW_ZXI(z30,p0,bq3,SIZE);

					FMLA_ZSP( z4,p0,z0,z24);
					FMLA_ZSP( z5,p0,z0,z25);
					FMLA_ZSP( z7,p0,z1,z22);

					LD1W_ZXI(z0,p0,aa,2);

					FMLA_ZSP( z8,p0,z1,z23);
					FMLA_ZSP( z9,p0,z1,z24);

					FMLA_ZSP( z10,p0,z1,z25);

					LD1W_ZXI(z1,p1,aa,3);
					aa+=nvl*3+mm;
					bq0+=2;

					FMLA_ZSP( z12,p0,z0,z22);
					FMLA_ZSP( z13,p0,z0,z23);

					bq1+=2;
					bq2+=2;

					FMLA_ZSP( z14,p0,z0,z24);
					FMLA_ZSP( z15,p0,z0,z25);
					FMLA_ZSP( z17,p1,z1,z22);

					LD1W_ZXI(z0,p0,aa,0);
					bq3+=2;

					FMLA_ZSP( z18,p1,z1,z23);
					FMLA_ZSP( z19,p1,z1,z24);

					LD1RW_ZXI(z22,p0,bq0,0);
					LD1RW_ZXI(z23,p0,bq1,0);

					FMLA_ZSP( z20,p1,z1,z25);

					LD1W_ZXI(z1,p0,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p0,z0,z27);
					FMLA_ZSP( z3,p0,z0,z28);

					LD1RW_ZXI(z24,p0,bq2,0);
					LD1RW_ZXI(z25,p0,bq3,0);

					FMLA_ZSP( z4,p0,z0,z29);
					FMLA_ZSP( z5,p0,z0,z30);
					FMLA_ZSP( z7,p0,z1,z27);

					LD1W_ZXI(z0,p0,aa,2);

					FMLA_ZSP( z8,p0,z1,z28);
					FMLA_ZSP( z9,p0,z1,z29);
					FMLA_ZSP( z10,p0,z1,z30);

					LD1W_ZXI(z1,p1,aa,3);
					aa+=nvl*3+mm;

					FMLA_ZSP( z12,p0,z0,z27);
					FMLA_ZSP( z13,p0,z0,z28);

					ik+=2;

					FMLA_ZSP( z14,p0,z0,z29);
					FMLA_ZSP( z15,p0,z0,z30);
					FMLA_ZSP( z17,p1,z1,z27);

					LD1W_ZXI(z0,p0,aa,0);
					NOP;

					FMLA_ZSP( z18,p1,z1,z28);
					FMLA_ZSP( z19,p1,z1,z29);

					LD1RW_ZXI(z27,p0,bq0,SIZE);
					LD1RW_ZXI(z28,p0,bq1,SIZE);

					FMLA_ZSP( z20,p1,z1,z30);
				}
				LD1W_ZXI(z1,p0,aa,1);

				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z3,p0,z0,z23);

				LD1RW_ZXI(z29,p0,bq2,SIZE);
				LD1RW_ZXI(z30,p0,bq3,SIZE);

				FMLA_ZSP( z4,p0,z0,z24);
				FMLA_ZSP( z5,p0,z0,z25);

				LD1W_ZXI(z0,p0,aa,2);

				FMLA_ZSP( z7,p0,z1,z22);
				FMLA_ZSP( z8,p0,z1,z23);
				FMLA_ZSP( z9,p0,z1,z24);
				FMLA_ZSP( z10,p0,z1,z25);

				LD1W_ZXI(z1,p1,aa,3);
				aa+=nvl*3+mm;

				FMLA_ZSP( z12,p0,z0,z22);
				FMLA_ZSP( z13,p0,z0,z23);
				FMLA_ZSP( z14,p0,z0,z24);
				FMLA_ZSP( z15,p0,z0,z25);

				LD1W_ZXI(z0,p0,aa,0);

				FMLA_ZSP( z17,p1,z1,z22);
				FMLA_ZSP( z18,p1,z1,z23);
				FMLA_ZSP( z19,p1,z1,z24);
				FMLA_ZSP( z20,p1,z1,z25);

				LD1W_ZXI(z1,p0,aa,1);

				FMLA_ZSP( z2,p0,z0,z27);
				FMLA_ZSP( z3,p0,z0,z28);
				FMLA_ZSP( z4,p0,z0,z29);
				FMLA_ZSP( z5,p0,z0,z30);

				LD1W_ZXI(z0,p0,aa,2);

				FMLA_ZSP( z7,p0,z1,z27);
				FMLA_ZSP( z8,p0,z1,z28);
				FMLA_ZSP( z9,p0,z1,z29);
				FMLA_ZSP( z10,p0,z1,z30);

				LD1W_ZXI(z1,p1,aa,3);

				FMLA_ZSP( z12,p0,z0,z27);
				FMLA_ZSP( z13,p0,z0,z28);
				FMLA_ZSP( z14,p0,z0,z29);
				FMLA_ZSP( z15,p0,z0,z30);
				FMLA_ZSP( z17,p1,z1,z27);
				FMLA_ZSP( z18,p1,z1,z28);
				FMLA_ZSP( z19,p1,z1,z29);
				FMLA_ZSP( z20,p1,z1,z30);

				bq0+=2;
				bq1+=2;
				bq2+=2;
				bq3+=2;
				aa+=nvl*3+mm;

			}
			if(k-ik>=1) {
				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);
				LD1RW_ZXI(z24,p0,bq2,0);
				LD1RW_ZXI(z25,p0,bq3,0);

				LD1W_ZXI(z0,p0,aa,0);
				LD1W_ZXI(z1,p0,aa,1);

				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z3,p0,z0,z23);
				FMLA_ZSP( z4,p0,z0,z24);
				FMLA_ZSP( z5,p0,z0,z25);
				FMLA_ZSP( z7,p0,z1,z22);
				FMLA_ZSP( z8,p0,z1,z23);
				FMLA_ZSP( z9,p0,z1,z24);
				FMLA_ZSP( z10,p0,z1,z25);

				LD1W_ZXI(z0,p0,aa,2);
				LD1W_ZXI(z1,p1,aa,3);

				FMLA_ZSP( z12,p0,z0,z22);
				FMLA_ZSP( z13,p0,z0,z23);
				FMLA_ZSP( z14,p0,z0,z24);
				FMLA_ZSP( z15,p0,z0,z25);
				FMLA_ZSP( z17,p1,z1,z22);
				FMLA_ZSP( z18,p1,z1,z23);
				FMLA_ZSP( z19,p1,z1,z24);
				FMLA_ZSP( z20,p1,z1,z25);

				aa+=nvl*3+mm;
				bq0+=1;
				bq1+=1;
				bq2+=1;
				bq3+=1;
			}
			LD1RW_ZXI(z0,p0,alphap,0);
			CALPHA_P_44;
			cq0+=nvl*3+mm;
			cq1+=nvl*3+mm;
			cq2+=nvl*3+mm;
			cq3+=nvl*3+mm;
		}

		else if(m-im>nvl*2) {
			mm=m-im;
			im+=nvl*2;
			WHILELT_PSX(p1,im,m);
			mm=m-im;
			im=m;
			PRFM_XI(PSTL1KEEP,cq0,256*2);
			PRFM_XI(PSTL1KEEP,cq1,256*2);
			PRFM_XI(PSTL1KEEP,cq2,256*2);
			PRFM_XI(PSTL1KEEP,cq3,256*2);

			CZERO_34;

			PRFM_XI(PLDL1KEEP,aa,0);
			PRFM_XI(PLDL1KEEP,aa,256*1);
			PRFM_XI(PLDL1KEEP,aa,256*2);
			PRFM_XI(PLDL1KEEP,aa,256*3);
			PRFM_XI(PLDL1KEEP,aa,256*4);
			PRFM_XI(PLDL1KEEP,aa,256*5);
			PRFM_XI(PLDL1KEEP,aa,256*6);
			PRFM_XI(PLDL1KEEP,aa,256*7);
			PRFM_XI(PLDL1KEEP,aa,256*8);

			bq0=bs0;
			bq1=bs1;
			bq2=bs2;
			bq3=bs3;
			ik=0;
			if(k>=2) { 

				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);
				LD1RW_ZXI(z24,p0,bq2,0);
				LD1RW_ZXI(z25,p0,bq3,0);

				LD1RW_ZXI(z27,p0,bq0,SIZE);
				LD1RW_ZXI(z28,p0,bq1,SIZE);

				LD1W_ZXI(z0,p0,aa,0);

				ik+=2;

				for(;ik<k-1;) {
                                        PRFM_XI(PLDL1KEEP,bq0,256);
                                        PRFM_XI(PLDL1KEEP,bq1,256);
                                        PRFM_XI(PLDL1KEEP,bq2,256);
                                        PRFM_XI(PLDL1KEEP,bq3,256);

					LD1W_ZXI(z1,p0,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p0,z0,z22);
					FMLA_ZSP( z3,p0,z0,z23);

					LD1RW_ZXI(z29,p0,bq2,SIZE);
					LD1RW_ZXI(z30,p0,bq3,SIZE);

					FMLA_ZSP( z4,p0,z0,z24);
					FMLA_ZSP( z5,p0,z0,z25);
					FMLA_ZSP( z7,p0,z1,z22);

					LD1W_ZXI(z17,p1,aa,2);

					FMLA_ZSP( z8,p0,z1,z23);
					FMLA_ZSP( z9,p0,z1,z24);

					FMLA_ZSP( z10,p0,z1,z25);

					aa+=nvl*2+mm;
					bq0+=2;

					FMLA_ZSP( z12,p1,z17,z22);
					FMLA_ZSP( z13,p1,z17,z23);

					bq1+=2;
					bq2+=2;

					FMLA_ZSP( z14,p1,z17,z24);
					FMLA_ZSP( z15,p1,z17,z25);

					LD1W_ZXI(z0,p0,aa,0);
					bq3+=2;

					LD1RW_ZXI(z22,p0,bq0,0);
					LD1RW_ZXI(z23,p0,bq1,0);

					LD1W_ZXI(z1,p0,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p0,z0,z27);
					FMLA_ZSP( z3,p0,z0,z28);

					LD1RW_ZXI(z24,p0,bq2,0);
					LD1RW_ZXI(z25,p0,bq3,0);

					FMLA_ZSP( z4,p0,z0,z29);
					FMLA_ZSP( z5,p0,z0,z30);
					FMLA_ZSP( z7,p0,z1,z27);

					LD1W_ZXI(z17,p1,aa,2);

					FMLA_ZSP( z8,p0,z1,z28);
					FMLA_ZSP( z9,p0,z1,z29);
					FMLA_ZSP( z10,p0,z1,z30);

					aa+=nvl*2+mm;

					FMLA_ZSP( z12,p1,z17,z27);
					FMLA_ZSP( z13,p1,z17,z28);

					ik+=2;

					FMLA_ZSP( z14,p1,z17,z29);
					FMLA_ZSP( z15,p1,z17,z30);

					LD1W_ZXI(z0,p0,aa,0);
					NOP;

					LD1RW_ZXI(z27,p0,bq0,SIZE);
					LD1RW_ZXI(z28,p0,bq1,SIZE);

				}
				LD1W_ZXI(z1,p0,aa,1);

				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z3,p0,z0,z23);

				LD1RW_ZXI(z29,p0,bq2,SIZE);
				LD1RW_ZXI(z30,p0,bq3,SIZE);

				FMLA_ZSP( z4,p0,z0,z24);
				FMLA_ZSP( z5,p0,z0,z25);

				LD1W_ZXI(z17,p1,aa,2);

				FMLA_ZSP( z7,p0,z1,z22);
				FMLA_ZSP( z8,p0,z1,z23);
				FMLA_ZSP( z9,p0,z1,z24);
				FMLA_ZSP( z10,p0,z1,z25);

				aa+=nvl*2+mm;

				FMLA_ZSP( z12,p1,z17,z22);
				FMLA_ZSP( z13,p1,z17,z23);
				FMLA_ZSP( z14,p1,z17,z24);
				FMLA_ZSP( z15,p1,z17,z25);

				LD1W_ZXI(z0,p0,aa,0);

				LD1W_ZXI(z1,p0,aa,1);

				FMLA_ZSP( z2,p0,z0,z27);
				FMLA_ZSP( z3,p0,z0,z28);
				FMLA_ZSP( z4,p0,z0,z29);
				FMLA_ZSP( z5,p0,z0,z30);

				LD1W_ZXI(z17,p1,aa,2);

				FMLA_ZSP( z7,p0,z1,z27);
				FMLA_ZSP( z8,p0,z1,z28);
				FMLA_ZSP( z9,p0,z1,z29);
				FMLA_ZSP( z10,p0,z1,z30);

				FMLA_ZSP( z12,p1,z17,z27);
				FMLA_ZSP( z13,p1,z17,z28);
				FMLA_ZSP( z14,p1,z17,z29);
				FMLA_ZSP( z15,p1,z17,z30);

				bq0+=2;
				bq1+=2;
				bq2+=2;
				bq3+=2;
				aa+=nvl*2+mm;

			}
			if(k-ik>=1) {
				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);
				LD1RW_ZXI(z24,p0,bq2,0);
				LD1RW_ZXI(z25,p0,bq3,0);

				LD1W_ZXI(z0,p0,aa,0);
				LD1W_ZXI(z1,p0,aa,1);
				LD1W_ZXI(z17,p1,aa,2);

				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z3,p0,z0,z23);
				FMLA_ZSP( z4,p0,z0,z24);
				FMLA_ZSP( z5,p0,z0,z25);

				FMLA_ZSP( z7,p0,z1,z22);
				FMLA_ZSP( z8,p0,z1,z23);
				FMLA_ZSP( z9,p0,z1,z24);
				FMLA_ZSP( z10,p0,z1,z25);

				FMLA_ZSP( z12,p1,z17,z22);
				FMLA_ZSP( z13,p1,z17,z23);
				FMLA_ZSP( z14,p1,z17,z24);
				FMLA_ZSP( z15,p1,z17,z25);

				aa+=nvl*2+mm;
				bq0+=1;
				bq1+=1;
				bq2+=1;
				bq3+=1;
			}
			LD1RW_ZXI(z0,p0,alphap,0);
			CALPHA_P_34;
			cq0+=nvl*2+mm;
			cq1+=nvl*2+mm;
			cq2+=nvl*2+mm;
			cq3+=nvl*2+mm;
		}

		else if(m-im>nvl*1) {	
			mm=m-im;
			im+=nvl*1;
			WHILELT_PSX(p1,im,m);
			mm=m-im;
			im=m;
			PRFM_XI(PSTL1KEEP,cq0,256*2);
			PRFM_XI(PSTL1KEEP,cq1,256*2);
			PRFM_XI(PSTL1KEEP,cq2,256*2);
			PRFM_XI(PSTL1KEEP,cq3,256*2);

			CZERO_24;

			PRFM_XI(PLDL1KEEP,aa,0);
			PRFM_XI(PLDL1KEEP,aa,256*1);
			PRFM_XI(PLDL1KEEP,aa,256*2);
			PRFM_XI(PLDL1KEEP,aa,256*3);
			PRFM_XI(PLDL1KEEP,aa,256*4);
			PRFM_XI(PLDL1KEEP,aa,256*5);
			PRFM_XI(PLDL1KEEP,aa,256*6);
			PRFM_XI(PLDL1KEEP,aa,256*7);
			PRFM_XI(PLDL1KEEP,aa,256*8);

			bq0=bs0;
			bq1=bs1;
			bq2=bs2;
			bq3=bs3;
			ik=0;
			if(k>=2) { 

				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);
				LD1RW_ZXI(z24,p0,bq2,0);
				LD1RW_ZXI(z25,p0,bq3,0);

				LD1RW_ZXI(z27,p0,bq0,SIZE);
				LD1RW_ZXI(z28,p0,bq1,SIZE);

				LD1W_ZXI(z0,p0,aa,0);

				ik+=2;

				for(;ik<k-1;) {
                                        PRFM_XI(PLDL1KEEP,bq0,256);
                                        PRFM_XI(PLDL1KEEP,bq1,256);
                                        PRFM_XI(PLDL1KEEP,bq2,256);
                                        PRFM_XI(PLDL1KEEP,bq3,256);

					LD1W_ZXI(z1,p1,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p0,z0,z22);
					FMLA_ZSP( z3,p0,z0,z23);

					LD1RW_ZXI(z29,p0,bq2,SIZE);
					LD1RW_ZXI(z30,p0,bq3,SIZE);

					FMLA_ZSP( z4,p0,z0,z24);
					FMLA_ZSP( z5,p0,z0,z25);
					FMLA_ZSP( z7,p1,z1,z22);

					FMLA_ZSP( z8,p1,z1,z23);
					FMLA_ZSP( z9,p1,z1,z24);

					FMLA_ZSP( z10,p1,z1,z25);

					aa+=nvl*1+mm;
					bq0+=2;
					bq1+=2;
					bq2+=2;

					LD1W_ZXI(z0,p0,aa,0);
					bq3+=2;

					LD1RW_ZXI(z22,p0,bq0,0);
					LD1RW_ZXI(z23,p0,bq1,0);

					LD1W_ZXI(z1,p1,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p0,z0,z27);
					FMLA_ZSP( z3,p0,z0,z28);

					LD1RW_ZXI(z24,p0,bq2,0);
					LD1RW_ZXI(z25,p0,bq3,0);

					FMLA_ZSP( z4,p0,z0,z29);
					FMLA_ZSP( z5,p0,z0,z30);
					FMLA_ZSP( z7,p1,z1,z27);

					FMLA_ZSP( z8,p1,z1,z28);
					FMLA_ZSP( z9,p1,z1,z29);
					FMLA_ZSP( z10,p1,z1,z30);

					aa+=nvl*1+mm;

					ik+=2;

					LD1W_ZXI(z0,p0,aa,0);
					NOP;

					LD1RW_ZXI(z27,p0,bq0,SIZE);
					LD1RW_ZXI(z28,p0,bq1,SIZE);

				}
				LD1W_ZXI(z1,p1,aa,1);

				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z3,p0,z0,z23);

				LD1RW_ZXI(z29,p0,bq2,SIZE);
				LD1RW_ZXI(z30,p0,bq3,SIZE);

				FMLA_ZSP( z4,p0,z0,z24);
				FMLA_ZSP( z5,p0,z0,z25);

				FMLA_ZSP( z7,p1,z1,z22);
				FMLA_ZSP( z8,p1,z1,z23);
				FMLA_ZSP( z9,p1,z1,z24);
				FMLA_ZSP( z10,p1,z1,z25);

				aa+=nvl*1+mm;

				LD1W_ZXI(z0,p0,aa,0);

				LD1W_ZXI(z1,p1,aa,1);

				FMLA_ZSP( z2,p0,z0,z27);
				FMLA_ZSP( z3,p0,z0,z28);
				FMLA_ZSP( z4,p0,z0,z29);
				FMLA_ZSP( z5,p0,z0,z30);

				FMLA_ZSP( z7,p1,z1,z27);
				FMLA_ZSP( z8,p1,z1,z28);
				FMLA_ZSP( z9,p1,z1,z29);
				FMLA_ZSP( z10,p1,z1,z30);

				bq0+=2;
				bq1+=2;
				bq2+=2;
				bq3+=2;
				aa+=nvl*1+mm;

			}
			if(k-ik>=1) {
				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);
				LD1RW_ZXI(z24,p0,bq2,0);
				LD1RW_ZXI(z25,p0,bq3,0);

				LD1W_ZXI(z0,p0,aa,0);
				LD1W_ZXI(z1,p1,aa,1);

				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z3,p0,z0,z23);
				FMLA_ZSP( z4,p0,z0,z24);
				FMLA_ZSP( z5,p0,z0,z25);

				FMLA_ZSP( z7,p1,z1,z22);
				FMLA_ZSP( z8,p1,z1,z23);
				FMLA_ZSP( z9,p1,z1,z24);
				FMLA_ZSP( z10,p1,z1,z25);

				aa+=nvl*1+mm;
				bq0+=1;
				bq1+=1;
				bq2+=1;
				bq3+=1;
			}
			LD1RW_ZXI(z0,p0,alphap,0);
			CALPHA_P_24;
			cq0+=nvl*1+mm;
			cq1+=nvl*1+mm;
			cq2+=nvl*1+mm;
			cq3+=nvl*1+mm;
		}


		else if(m-im>0) {
			mm=m-im;
			WHILELT_PSX(p1,im,m);
			mm=m-im;
			im=m;
			PRFM_XI(PSTL1KEEP,cq0,256*2);
			PRFM_XI(PSTL1KEEP,cq1,256*2);
			PRFM_XI(PSTL1KEEP,cq2,256*2);
			PRFM_XI(PSTL1KEEP,cq3,256*2);

			CZERO_14;

			PRFM_XI(PLDL1KEEP,aa,0);
			PRFM_XI(PLDL1KEEP,aa,256*1);
			PRFM_XI(PLDL1KEEP,aa,256*2);
			PRFM_XI(PLDL1KEEP,aa,256*3);
			PRFM_XI(PLDL1KEEP,aa,256*4);
			PRFM_XI(PLDL1KEEP,aa,256*5);
			PRFM_XI(PLDL1KEEP,aa,256*6);
			PRFM_XI(PLDL1KEEP,aa,256*7);
			PRFM_XI(PLDL1KEEP,aa,256*8);

			bq0=bs0;
			bq1=bs1;
			bq2=bs2;
			bq3=bs3;
			ik=0;
			if(k>=2) { 

				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);
				LD1RW_ZXI(z24,p0,bq2,0);
				LD1RW_ZXI(z25,p0,bq3,0);

				LD1RW_ZXI(z27,p0,bq0,SIZE);
				LD1RW_ZXI(z28,p0,bq1,SIZE);

				LD1W_ZXI(z0,p1,aa,0);

				ik+=2;

				for(;ik<k-1;) {
                                        PRFM_XI(PLDL1KEEP,bq0,256);
                                        PRFM_XI(PLDL1KEEP,bq1,256);
                                        PRFM_XI(PLDL1KEEP,bq2,256);
                                        PRFM_XI(PLDL1KEEP,bq3,256);

					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p1,z0,z22);
					FMLA_ZSP( z3,p1,z0,z23);

					LD1RW_ZXI(z29,p0,bq2,SIZE);
					LD1RW_ZXI(z30,p0,bq3,SIZE);

					FMLA_ZSP( z4,p1,z0,z24);
					FMLA_ZSP( z5,p1,z0,z25);

					aa+=mm;
					bq0+=2;
					bq1+=2;
					bq2+=2;

					LD1W_ZXI(z0,p1,aa,0);
					bq3+=2;

					LD1RW_ZXI(z22,p0,bq0,0);
					LD1RW_ZXI(z23,p0,bq1,0);

					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p1,z0,z27);
					FMLA_ZSP( z3,p1,z0,z28);

					LD1RW_ZXI(z24,p0,bq2,0);
					LD1RW_ZXI(z25,p0,bq3,0);

					FMLA_ZSP( z4,p1,z0,z29);
					FMLA_ZSP( z5,p1,z0,z30);

					aa+=mm;
					ik+=2;

					LD1W_ZXI(z0,p1,aa,0);
					NOP;

					LD1RW_ZXI(z27,p0,bq0,SIZE);
					LD1RW_ZXI(z28,p0,bq1,SIZE);

				}

				FMLA_ZSP( z2,p1,z0,z22);
				FMLA_ZSP( z3,p1,z0,z23);

				LD1RW_ZXI(z29,p0,bq2,SIZE);
				LD1RW_ZXI(z30,p0,bq3,SIZE);

				FMLA_ZSP( z4,p1,z0,z24);
				FMLA_ZSP( z5,p1,z0,z25);

				aa+=mm;

				LD1W_ZXI(z0,p1,aa,0);

				FMLA_ZSP( z2,p1,z0,z27);
				FMLA_ZSP( z3,p1,z0,z28);
				FMLA_ZSP( z4,p1,z0,z29);
				FMLA_ZSP( z5,p1,z0,z30);

				bq0+=2;
				bq1+=2;
				bq2+=2;
				bq3+=2;
				aa+=mm;

			}
			if(k-ik>=1) {
				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);
				LD1RW_ZXI(z24,p0,bq2,0);
				LD1RW_ZXI(z25,p0,bq3,0);

				LD1W_ZXI(z0,p1,aa,0);

				FMLA_ZSP( z2,p1,z0,z22);
				FMLA_ZSP( z3,p1,z0,z23);
				FMLA_ZSP( z4,p1,z0,z24);
				FMLA_ZSP( z5,p1,z0,z25);

				aa+=mm;
				bq0+=1;
				bq1+=1;
				bq2+=1;
				bq3+=1;

			}
			LD1RW_ZXI(z0,p0,alphap,0);
			CALPHA_P_14;
			cq0+=mm;
			cq1+=mm;
			cq2+=mm;
			cq3+=mm;
		}

		cq0+=ldc*4-m;
		cq1+=ldc*4-m;
		cq2+=ldc*4-m;
		cq3+=ldc*4-m;
		bs0+=ldb*4;
		bs1+=ldb*4;
		bs2+=ldb*4;
		bs3+=ldb*4;
	}

	else if(n-in==3) {
		in+=3;
		aa=a;
		im=0;
		for(;im<m-nul+1;im+=nul) {
			PRFM_XI(PSTL1KEEP,cq0,256*2);
			PRFM_XI(PSTL1KEEP,cq1,256*2);
			PRFM_XI(PSTL1KEEP,cq2,256*2);
			CZERO_43;

			PRFM_XI(PLDL1KEEP,aa,0);
			PRFM_XI(PLDL1KEEP,aa,256*1);
			PRFM_XI(PLDL1KEEP,aa,256*2);
			PRFM_XI(PLDL1KEEP,aa,256*3);
			PRFM_XI(PLDL1KEEP,aa,256*4);
			PRFM_XI(PLDL1KEEP,aa,256*5);
			PRFM_XI(PLDL1KEEP,aa,256*6);
			PRFM_XI(PLDL1KEEP,aa,256*7);
			PRFM_XI(PLDL1KEEP,aa,256*8);


			bq0=bs0;
			bq1=bs1;
			bq2=bs2;
			ik=0;
			if(k>=2) { 

				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);
				LD1RW_ZXI(z24,p0,bq2,0);

				LD1RW_ZXI(z27,p0,bq0,SIZE);
				LD1RW_ZXI(z28,p0,bq1,SIZE);

				LD1W_ZXI(z0,p0,aa,0);

				ik+=2;

				for(;ik<k-1;) {
                                        PRFM_XI(PLDL1KEEP,bq0,256);
                                        PRFM_XI(PLDL1KEEP,bq1,256);
                                        PRFM_XI(PLDL1KEEP,bq2,256);

					LD1W_ZXI(z1,p0,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p0,z0,z22);

					FMLA_ZSP( z3,p0,z0,z23);

					LD1RW_ZXI(z29,p0,bq2,SIZE);

					FMLA_ZSP( z4,p0,z0,z24);
					FMLA_ZSP( z7,p0,z1,z22);

					LD1W_ZXI(z0,p0,aa,2);

					FMLA_ZSP( z8,p0,z1,z23);
					FMLA_ZSP( z9,p0,z1,z24);

					LD1W_ZXI(z1,p0,aa,3);
					bq0+=2;

					FMLA_ZSP( z12,p0,z0,z22);
					FMLA_ZSP( z13,p0,z0,z23);

					bq1+=2;
					bq2+=2;

					FMLA_ZSP( z14,p0,z0,z24);
					FMLA_ZSP( z17,p0,z1,z22);

					LD1W_ZXI(z0,p0,aa,4);

					FMLA_ZSP( z18,p0,z1,z23);
					FMLA_ZSP( z19,p0,z1,z24);

					LD1RW_ZXI(z22,p0,bq0,0);
					LD1RW_ZXI(z23,p0,bq1,0);

					LD1W_ZXI(z1,p0,aa,5);
					PRFM_XI(PLDL1KEEP,aa,256*10);

					FMLA_ZSP( z2,p0,z0,z27);
					FMLA_ZSP( z3,p0,z0,z28);

					LD1RW_ZXI(z24,p0,bq2,0);

					FMLA_ZSP( z4,p0,z0,z29);
					FMLA_ZSP( z7,p0,z1,z27);

					LD1W_ZXI(z0,p0,aa,6);

					FMLA_ZSP( z8,p0,z1,z28);
					FMLA_ZSP( z9,p0,z1,z29);

					LD1W_ZXI(z1,p0,aa,7);
					aa+=nul*2;

					FMLA_ZSP( z12,p0,z0,z27);
					FMLA_ZSP( z13,p0,z0,z28);

					ik+=2;

					FMLA_ZSP( z14,p0,z0,z29);
					FMLA_ZSP( z17,p0,z1,z27);

					LD1W_ZXI(z0,p0,aa,0);
					NOP;

					FMLA_ZSP( z18,p0,z1,z28);
					FMLA_ZSP( z19,p0,z1,z29);

					LD1RW_ZXI(z27,p0,bq0,SIZE);
					LD1RW_ZXI(z28,p0,bq1,SIZE);
				}
				LD1W_ZXI(z1,p0,aa,1);

				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z3,p0,z0,z23);

				LD1RW_ZXI(z29,p0,bq2,SIZE);

				FMLA_ZSP( z4,p0,z0,z24);

				LD1W_ZXI(z0,p0,aa,2);

				FMLA_ZSP( z7 ,p0,z1,z22);
				FMLA_ZSP( z8 ,p0,z1,z23);
				FMLA_ZSP( z9 ,p0,z1,z24);

				LD1W_ZXI(z1,p0,aa,3);

				FMLA_ZSP( z12,p0,z0,z22);
				FMLA_ZSP( z13,p0,z0,z23);
				FMLA_ZSP( z14,p0,z0,z24);

				LD1W_ZXI(z0,p0,aa,4);

				FMLA_ZSP( z17,p0,z1,z22);
				FMLA_ZSP( z18,p0,z1,z23);
				FMLA_ZSP( z19,p0,z1,z24);

				LD1W_ZXI(z1,p0,aa,5);

				FMLA_ZSP( z2,p0,z0,z27);
				FMLA_ZSP( z3,p0,z0,z28);
				FMLA_ZSP( z4,p0,z0,z29);


				LD1W_ZXI(z0,p0,aa,6);

				FMLA_ZSP( z7,p0,z1,z27);
				FMLA_ZSP( z8,p0,z1,z28);
				FMLA_ZSP( z9,p0,z1,z29);

				LD1W_ZXI(z1,p0,aa,7);
				bq0+=2;
				bq1+=2;
				bq2+=2;
				aa+=nul*2;

				if(k-ik>=1) {
					LD1RW_ZXI(z22,p0,bq0,0);
					LD1RW_ZXI(z23,p0,bq1,0);

					FMLA_ZSP( z12,p0,z0,z27);
					FMLA_ZSP( z13,p0,z0,z28);

					LD1RW_ZXI(z24,p0,bq2,0);

					FMLA_ZSP( z14,p0,z0,z29);

					LD1W_ZXI(z0,p0,aa,0);

					FMLA_ZSP( z17,p0,z1,z27);
					FMLA_ZSP( z18,p0,z1,z28);
					FMLA_ZSP( z19,p0,z1,z29);

					LD1W_ZXI(z1,p0,aa,1);
				}
				else {
					FMLA_ZSP( z12,p0,z0,z27);
					FMLA_ZSP( z13,p0,z0,z28);
					FMLA_ZSP( z14,p0,z0,z29);
					FMLA_ZSP( z17,p0,z1,z27);
					FMLA_ZSP( z18,p0,z1,z28);
					FMLA_ZSP( z19,p0,z1,z29);
				}
			}
			else if(k==1) {
				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);
				LD1RW_ZXI(z24,p0,bq2,0);

				LD1W_ZXI(z0,p0,aa,0);
				LD1W_ZXI(z1,p0,aa,1);

			}

			if(k-ik>=1) {

				FMLA_ZSP( z2, p0,z0,z22);
				FMLA_ZSP( z3, p0,z0,z23);
				FMLA_ZSP( z4, p0,z0,z24);

				LD1W_ZXI(z0,p0,aa,2);

				FMLA_ZSP( z7, p0,z1,z22);
				FMLA_ZSP( z8, p0,z1,z23);
				FMLA_ZSP( z9, p0,z1,z24);

				LD1W_ZXI(z1,p0,aa,3);

				FMLA_ZSP( z12,p0,z0,z22);
				FMLA_ZSP( z13,p0,z0,z23);
				FMLA_ZSP( z14,p0,z0,z24);
				FMLA_ZSP( z17,p0,z1,z22);
				FMLA_ZSP( z18,p0,z1,z23);
				FMLA_ZSP( z19,p0,z1,z24);

				aa+=nul;
				bq0+=1;
				bq1+=1;
				bq2+=1;
				ik++;
			}
			LD1RW_ZXI(z0,p0,alphap,0);
			CALPHA_43;
			cq0+=nul;
			cq1+=nul;
			cq2+=nul;
		}

		if(m-im>nvl*3) {
			mm=m-im;
			im+=nvl*3;
			WHILELT_PSX(p1,im,m);
			mm=m-im;
			im=m;
			PRFM_XI(PSTL1KEEP,cq0,256*2);
			PRFM_XI(PSTL1KEEP,cq1,256*2);
			PRFM_XI(PSTL1KEEP,cq2,256*2);

			CZERO_43;

			PRFM_XI(PLDL1KEEP,aa,0);
			PRFM_XI(PLDL1KEEP,aa,256*1);
			PRFM_XI(PLDL1KEEP,aa,256*2);
			PRFM_XI(PLDL1KEEP,aa,256*3);
			PRFM_XI(PLDL1KEEP,aa,256*4);
			PRFM_XI(PLDL1KEEP,aa,256*5);
			PRFM_XI(PLDL1KEEP,aa,256*6);
			PRFM_XI(PLDL1KEEP,aa,256*7);
			PRFM_XI(PLDL1KEEP,aa,256*8);

			bq0=bs0;
			bq1=bs1;
			bq2=bs2;
			ik=0;
			if(k>=2) { 

				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);
				LD1RW_ZXI(z24,p0,bq2,0);

				LD1RW_ZXI(z27,p0,bq0,SIZE);
				LD1RW_ZXI(z28,p0,bq1,SIZE);

				LD1W_ZXI(z0,p0,aa,0);

				ik+=2;

				for(;ik<k-1;) {
                                        PRFM_XI(PLDL1KEEP,bq0,256);
                                        PRFM_XI(PLDL1KEEP,bq1,256);
                                        PRFM_XI(PLDL1KEEP,bq2,256);

					LD1W_ZXI(z1,p0,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p0,z0,z22);
					FMLA_ZSP( z3,p0,z0,z23);

					LD1RW_ZXI(z29,p0,bq2,SIZE);

					FMLA_ZSP( z4,p0,z0,z24);
					FMLA_ZSP( z7,p0,z1,z22);

					LD1W_ZXI(z0,p0,aa,2);

					FMLA_ZSP( z8,p0,z1,z23);
					FMLA_ZSP( z9,p0,z1,z24);

					LD1W_ZXI(z1,p1,aa,3);
					aa+=nvl*3+mm;
					bq0+=2;

					FMLA_ZSP( z12,p0,z0,z22);
					FMLA_ZSP( z13,p0,z0,z23);

					bq1+=2;
					bq2+=2;

					FMLA_ZSP( z14,p0,z0,z24);
					FMLA_ZSP( z17,p1,z1,z22);

					LD1W_ZXI(z0,p0,aa,0);

					FMLA_ZSP( z18,p1,z1,z23);
					FMLA_ZSP( z19,p1,z1,z24);

					LD1RW_ZXI(z22,p0,bq0,0);
					LD1RW_ZXI(z23,p0,bq1,0);

					LD1W_ZXI(z1,p0,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p0,z0,z27);
					FMLA_ZSP( z3,p0,z0,z28);

					LD1RW_ZXI(z24,p0,bq2,0);

					FMLA_ZSP( z4,p0,z0,z29);
					FMLA_ZSP( z7,p0,z1,z27);

					LD1W_ZXI(z0,p0,aa,2);

					FMLA_ZSP( z8,p0,z1,z28);
					FMLA_ZSP( z9,p0,z1,z29);

					LD1W_ZXI(z1,p1,aa,3);
					aa+=nvl*3+mm;

					FMLA_ZSP( z12,p0,z0,z27);
					FMLA_ZSP( z13,p0,z0,z28);

					ik+=2;

					FMLA_ZSP( z14,p0,z0,z29);
					FMLA_ZSP( z17,p1,z1,z27);

					LD1W_ZXI(z0,p0,aa,0);
					NOP;

					FMLA_ZSP( z18,p1,z1,z28);
					FMLA_ZSP( z19,p1,z1,z29);

					LD1RW_ZXI(z27,p0,bq0,SIZE);
					LD1RW_ZXI(z28,p0,bq1,SIZE);
				}
				LD1W_ZXI(z1,p0,aa,1);

				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z3,p0,z0,z23);

				LD1RW_ZXI(z29,p0,bq2,SIZE);

				FMLA_ZSP( z4,p0,z0,z24);

				LD1W_ZXI(z0,p0,aa,2);

				FMLA_ZSP( z7,p0,z1,z22);
				FMLA_ZSP( z8,p0,z1,z23);
				FMLA_ZSP( z9,p0,z1,z24);

				LD1W_ZXI(z1,p1,aa,3);
				aa+=nvl*3+mm;

				FMLA_ZSP( z12,p0,z0,z22);
				FMLA_ZSP( z13,p0,z0,z23);
				FMLA_ZSP( z14,p0,z0,z24);

				LD1W_ZXI(z0,p0,aa,0);

				FMLA_ZSP( z17,p1,z1,z22);
				FMLA_ZSP( z18,p1,z1,z23);
				FMLA_ZSP( z19,p1,z1,z24);

				LD1W_ZXI(z1,p0,aa,1);

				FMLA_ZSP( z2,p0,z0,z27);
				FMLA_ZSP( z3,p0,z0,z28);
				FMLA_ZSP( z4,p0,z0,z29);

				LD1W_ZXI(z0,p0,aa,2);

				FMLA_ZSP( z7,p0,z1,z27);
				FMLA_ZSP( z8,p0,z1,z28);
				FMLA_ZSP( z9,p0,z1,z29);

				LD1W_ZXI(z1,p1,aa,3);

				FMLA_ZSP( z12,p0,z0,z27);
				FMLA_ZSP( z13,p0,z0,z28);
				FMLA_ZSP( z14,p0,z0,z29);
				FMLA_ZSP( z17,p1,z1,z27);
				FMLA_ZSP( z18,p1,z1,z28);
				FMLA_ZSP( z19,p1,z1,z29);

				bq0+=2;
				bq1+=2;
				bq2+=2;
				aa+=nvl*3+mm;

			}
			if(k-ik>=1) {
				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);
				LD1RW_ZXI(z24,p0,bq2,0);

				LD1W_ZXI(z0,p0,aa,0);
				LD1W_ZXI(z1,p0,aa,1);

				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z3,p0,z0,z23);
				FMLA_ZSP( z4,p0,z0,z24);
				FMLA_ZSP( z7,p0,z1,z22);
				FMLA_ZSP( z8,p0,z1,z23);
				FMLA_ZSP( z9,p0,z1,z24);

				LD1W_ZXI(z0,p0,aa,2);
				LD1W_ZXI(z1,p1,aa,3);

				FMLA_ZSP( z12,p0,z0,z22);
				FMLA_ZSP( z13,p0,z0,z23);
				FMLA_ZSP( z14,p0,z0,z24);
				FMLA_ZSP( z17,p1,z1,z22);
				FMLA_ZSP( z18,p1,z1,z23);
				FMLA_ZSP( z19,p1,z1,z24);

				aa+=nvl*3+mm;
				bq0+=1;
				bq1+=1;
				bq2+=1;
			}
			LD1RW_ZXI(z0,p0,alphap,0);
			CALPHA_P_43;
			cq0+=nvl*3+mm;
			cq1+=nvl*3+mm;
			cq2+=nvl*3+mm;
		}

		else if(m-im>nvl*2) {
			mm=m-im;
			im+=nvl*2;
			WHILELT_PSX(p1,im,m);
			mm=m-im;
			im=m;
			PRFM_XI(PSTL1KEEP,cq0,256*2);
			PRFM_XI(PSTL1KEEP,cq1,256*2);
			PRFM_XI(PSTL1KEEP,cq2,256*2);

			CZERO_33;

			PRFM_XI(PLDL1KEEP,aa,0);
			PRFM_XI(PLDL1KEEP,aa,256*1);
			PRFM_XI(PLDL1KEEP,aa,256*2);
			PRFM_XI(PLDL1KEEP,aa,256*3);
			PRFM_XI(PLDL1KEEP,aa,256*4);
			PRFM_XI(PLDL1KEEP,aa,256*5);
			PRFM_XI(PLDL1KEEP,aa,256*6);
			PRFM_XI(PLDL1KEEP,aa,256*7);
			PRFM_XI(PLDL1KEEP,aa,256*8);

			bq0=bs0;
			bq1=bs1;
			bq2=bs2;
			ik=0;
			if(k>=2) { 

				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);
				LD1RW_ZXI(z24,p0,bq2,0);

				LD1RW_ZXI(z27,p0,bq0,SIZE);
				LD1RW_ZXI(z28,p0,bq1,SIZE);

				LD1W_ZXI(z0,p0,aa,0);

				ik+=2;

				for(;ik<k-1;) {
                                        PRFM_XI(PLDL1KEEP,bq0,256);
                                        PRFM_XI(PLDL1KEEP,bq1,256);
                                        PRFM_XI(PLDL1KEEP,bq2,256);

					LD1W_ZXI(z1,p0,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p0,z0,z22);
					FMLA_ZSP( z3,p0,z0,z23);

					LD1RW_ZXI(z29,p0,bq2,SIZE);

					FMLA_ZSP( z4,p0,z0,z24);
					FMLA_ZSP( z7,p0,z1,z22);

					LD1W_ZXI(z17,p1,aa,2);

					FMLA_ZSP( z8,p0,z1,z23);
					FMLA_ZSP( z9,p0,z1,z24);

					aa+=nvl*2+mm;
					bq0+=2;

					FMLA_ZSP( z12,p1,z17,z22);
					FMLA_ZSP( z13,p1,z17,z23);

					bq1+=2;
					bq2+=2;

					FMLA_ZSP( z14,p1,z17,z24);

					LD1W_ZXI(z0,p0,aa,0);

					LD1RW_ZXI(z22,p0,bq0,0);
					LD1RW_ZXI(z23,p0,bq1,0);

					LD1W_ZXI(z1,p0,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p0,z0,z27);
					FMLA_ZSP( z3,p0,z0,z28);

					LD1RW_ZXI(z24,p0,bq2,0);

					FMLA_ZSP( z4,p0,z0,z29);
					FMLA_ZSP( z7,p0,z1,z27);

					LD1W_ZXI(z17,p1,aa,2);

					FMLA_ZSP( z8,p0,z1,z28);
					FMLA_ZSP( z9,p0,z1,z29);

					aa+=nvl*2+mm;

					FMLA_ZSP( z12,p1,z17,z27);
					FMLA_ZSP( z13,p1,z17,z28);

					ik+=2;

					FMLA_ZSP( z14,p1,z17,z29);

					LD1W_ZXI(z0,p0,aa,0);
					NOP;

					LD1RW_ZXI(z27,p0,bq0,SIZE);
					LD1RW_ZXI(z28,p0,bq1,SIZE);

				}
				LD1W_ZXI(z1,p0,aa,1);

				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z3,p0,z0,z23);

				LD1RW_ZXI(z29,p0,bq2,SIZE);

				FMLA_ZSP( z4,p0,z0,z24);

				LD1W_ZXI(z17,p1,aa,2);

				FMLA_ZSP( z7,p0,z1,z22);
				FMLA_ZSP( z8,p0,z1,z23);
				FMLA_ZSP( z9,p0,z1,z24);

				aa+=nvl*2+mm;

				FMLA_ZSP( z12,p1,z17,z22);
				FMLA_ZSP( z13,p1,z17,z23);
				FMLA_ZSP( z14,p1,z17,z24);

				LD1W_ZXI(z0,p0,aa,0);

				LD1W_ZXI(z1,p0,aa,1);

				FMLA_ZSP( z2,p0,z0,z27);
				FMLA_ZSP( z3,p0,z0,z28);
				FMLA_ZSP( z4,p0,z0,z29);

				LD1W_ZXI(z17,p1,aa,2);

				FMLA_ZSP( z7,p0,z1,z27);
				FMLA_ZSP( z8,p0,z1,z28);
				FMLA_ZSP( z9,p0,z1,z29);

				FMLA_ZSP( z12,p1,z17,z27);
				FMLA_ZSP( z13,p1,z17,z28);
				FMLA_ZSP( z14,p1,z17,z29);

				bq0+=2;
				bq1+=2;
				bq2+=2;
				aa+=nvl*2+mm;

			}
			if(k-ik>=1) {
				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);
				LD1RW_ZXI(z24,p0,bq2,0);

				LD1W_ZXI(z0,p0,aa,0);
				LD1W_ZXI(z1,p0,aa,1);
				LD1W_ZXI(z17,p1,aa,2);

				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z3,p0,z0,z23);
				FMLA_ZSP( z4,p0,z0,z24);

				FMLA_ZSP( z7,p0,z1,z22);
				FMLA_ZSP( z8,p0,z1,z23);
				FMLA_ZSP( z9,p0,z1,z24);

				FMLA_ZSP( z12,p1,z17,z22);
				FMLA_ZSP( z13,p1,z17,z23);
				FMLA_ZSP( z14,p1,z17,z24);

				aa+=nvl*2+mm;
				bq0+=1;
				bq1+=1;
				bq2+=1;
			}
			LD1RW_ZXI(z0,p0,alphap,0);
			CALPHA_P_33;
			cq0+=nvl*2+mm;
			cq1+=nvl*2+mm;
			cq2+=nvl*2+mm;
		}

		else if(m-im>nvl*1) {
			mm=m-im;
			im+=nvl*1;
			WHILELT_PSX(p1,im,m);
			mm=m-im;
			im=m;
			PRFM_XI(PSTL1KEEP,cq0,256*2);
			PRFM_XI(PSTL1KEEP,cq1,256*2);
			PRFM_XI(PSTL1KEEP,cq2,256*2);

			CZERO_23;

			PRFM_XI(PLDL1KEEP,aa,0);
			PRFM_XI(PLDL1KEEP,aa,256*1);
			PRFM_XI(PLDL1KEEP,aa,256*2);
			PRFM_XI(PLDL1KEEP,aa,256*3);
			PRFM_XI(PLDL1KEEP,aa,256*4);
			PRFM_XI(PLDL1KEEP,aa,256*5);
			PRFM_XI(PLDL1KEEP,aa,256*6);
			PRFM_XI(PLDL1KEEP,aa,256*7);
			PRFM_XI(PLDL1KEEP,aa,256*8);

			bq0=bs0;
			bq1=bs1;
			bq2=bs2;
			ik=0;
			if(k>=2) { 

				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);
				LD1RW_ZXI(z24,p0,bq2,0);

				LD1RW_ZXI(z27,p0,bq0,SIZE);
				LD1RW_ZXI(z28,p0,bq1,SIZE);

				LD1W_ZXI(z0,p0,aa,0);

				ik+=2;

				for(;ik<k-1;) {
                                        PRFM_XI(PLDL1KEEP,bq0,256);
                                        PRFM_XI(PLDL1KEEP,bq1,256);
                                        PRFM_XI(PLDL1KEEP,bq2,256);

					LD1W_ZXI(z1,p1,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p0,z0,z22);
					FMLA_ZSP( z3,p0,z0,z23);

					LD1RW_ZXI(z29,p0,bq2,SIZE);

					FMLA_ZSP( z4,p0,z0,z24);
					FMLA_ZSP( z7,p1,z1,z22);

					FMLA_ZSP( z8,p1,z1,z23);
					FMLA_ZSP( z9,p1,z1,z24);

					aa+=nvl*1+mm;
					bq0+=2;
					bq1+=2;
					bq2+=2;

					LD1W_ZXI(z0,p0,aa,0);

					LD1RW_ZXI(z22,p0,bq0,0);
					LD1RW_ZXI(z23,p0,bq1,0);

					LD1W_ZXI(z1,p1,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p0,z0,z27);
					FMLA_ZSP( z3,p0,z0,z28);

					LD1RW_ZXI(z24,p0,bq2,0);

					FMLA_ZSP( z4,p0,z0,z29);
					FMLA_ZSP( z7,p1,z1,z27);

					FMLA_ZSP( z8,p1,z1,z28);
					FMLA_ZSP( z9,p1,z1,z29);

					aa+=nvl*1+mm;

					ik+=2;

					LD1W_ZXI(z0,p0,aa,0);
					NOP;

					LD1RW_ZXI(z27,p0,bq0,SIZE);
					LD1RW_ZXI(z28,p0,bq1,SIZE);

				}
				LD1W_ZXI(z1,p1,aa,1);

				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z3,p0,z0,z23);

				LD1RW_ZXI(z29,p0,bq2,SIZE);

				FMLA_ZSP( z4,p0,z0,z24);

				FMLA_ZSP( z7,p1,z1,z22);
				FMLA_ZSP( z8,p1,z1,z23);
				FMLA_ZSP( z9,p1,z1,z24);

				aa+=nvl*1+mm;

				LD1W_ZXI(z0,p0,aa,0);

				LD1W_ZXI(z1,p1,aa,1);

				FMLA_ZSP( z2,p0,z0,z27);
				FMLA_ZSP( z3,p0,z0,z28);
				FMLA_ZSP( z4,p0,z0,z29);

				FMLA_ZSP( z7,p1,z1,z27);
				FMLA_ZSP( z8,p1,z1,z28);
				FMLA_ZSP( z9,p1,z1,z29);

				bq0+=2;
				bq1+=2;
				bq2+=2;
				aa+=nvl*1+mm;

			}
			if(k-ik>=1) {
				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);
				LD1RW_ZXI(z24,p0,bq2,0);

				LD1W_ZXI(z0,p0,aa,0);
				LD1W_ZXI(z1,p1,aa,1);

				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z3,p0,z0,z23);
				FMLA_ZSP( z4,p0,z0,z24);

				FMLA_ZSP( z7,p1,z1,z22);
				FMLA_ZSP( z8,p1,z1,z23);
				FMLA_ZSP( z9,p1,z1,z24);

				aa+=nvl*1+mm;
				bq0+=1;
				bq1+=1;
				bq2+=1;
			}
			LD1RW_ZXI(z0,p0,alphap,0);
			CALPHA_P_23;
			cq0+=nvl*1+mm;
			cq1+=nvl*1+mm;
			cq2+=nvl*1+mm;
		}


		else if(m-im>0) {
			mm=m-im;
			WHILELT_PSX(p1,im,m);
			mm=m-im;
			im=m;
			PRFM_XI(PSTL1KEEP,cq0,256*2);
			PRFM_XI(PSTL1KEEP,cq1,256*2);
			PRFM_XI(PSTL1KEEP,cq2,256*2);

			CZERO_13;

			PRFM_XI(PLDL1KEEP,aa,0);
			PRFM_XI(PLDL1KEEP,aa,256*1);
			PRFM_XI(PLDL1KEEP,aa,256*2);
			PRFM_XI(PLDL1KEEP,aa,256*3);
			PRFM_XI(PLDL1KEEP,aa,256*4);
			PRFM_XI(PLDL1KEEP,aa,256*5);
			PRFM_XI(PLDL1KEEP,aa,256*6);
			PRFM_XI(PLDL1KEEP,aa,256*7);
			PRFM_XI(PLDL1KEEP,aa,256*8);

			bq0=bs0;
			bq1=bs1;
			bq2=bs2;
			ik=0;
			if(k>=2) { 

				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);
				LD1RW_ZXI(z24,p0,bq2,0);

				LD1RW_ZXI(z27,p0,bq0,SIZE);
				LD1RW_ZXI(z28,p0,bq1,SIZE);

				LD1W_ZXI(z0,p1,aa,0);

				ik+=2;

				for(;ik<k-1;) {
                                        PRFM_XI(PLDL1KEEP,bq0,256);
                                        PRFM_XI(PLDL1KEEP,bq1,256);
                                        PRFM_XI(PLDL1KEEP,bq2,256);

					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p1,z0,z22);
					FMLA_ZSP( z3,p1,z0,z23);

					LD1RW_ZXI(z29,p0,bq2,SIZE);

					FMLA_ZSP( z4,p1,z0,z24);

					aa+=mm;
					bq0+=2;
					bq1+=2;
					bq2+=2;

					LD1W_ZXI(z0,p1,aa,0);

					LD1RW_ZXI(z22,p0,bq0,0);
					LD1RW_ZXI(z23,p0,bq1,0);

					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p1,z0,z27);
					FMLA_ZSP( z3,p1,z0,z28);

					LD1RW_ZXI(z24,p0,bq2,0);

					FMLA_ZSP( z4,p1,z0,z29);

					aa+=mm;
					ik+=2;

					LD1W_ZXI(z0,p1,aa,0);
					NOP;

					LD1RW_ZXI(z27,p0,bq0,SIZE);
					LD1RW_ZXI(z28,p0,bq1,SIZE);

				}

				FMLA_ZSP( z2,p1,z0,z22);
				FMLA_ZSP( z3,p1,z0,z23);

				LD1RW_ZXI(z29,p0,bq2,SIZE);

				FMLA_ZSP( z4,p1,z0,z24);

				aa+=mm;

				LD1W_ZXI(z0,p1,aa,0);

				FMLA_ZSP( z2,p1,z0,z27);
				FMLA_ZSP( z3,p1,z0,z28);
				FMLA_ZSP( z4,p1,z0,z29);

				bq0+=2;
				bq1+=2;
				bq2+=2;
				aa+=mm;

			}
			if(k-ik>=1) {
				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);
				LD1RW_ZXI(z24,p0,bq2,0);

				LD1W_ZXI(z0,p1,aa,0);

				FMLA_ZSP( z2,p1,z0,z22);
				FMLA_ZSP( z3,p1,z0,z23);
				FMLA_ZSP( z4,p1,z0,z24);

				aa+=mm;
				bq0+=1;
				bq1+=1;
				bq2+=1;
			}
			LD1RW_ZXI(z0,p0,alphap,0);
			CALPHA_P_13;
			cq0+=mm;
			cq1+=mm;
			cq2+=mm;
		}

		cq0+=ldc*3-m;
		cq1+=ldc*3-m;
		cq2+=ldc*3-m;
		bs0+=ldb*3;
		bs1+=ldb*3;
		bs2+=ldb*3;
	}
	else if(n-in==2) {
		in+=2;
		aa=a;
		im=0;
		for(;im<m-nul+1;im+=nul) {
			PRFM_XI(PSTL1KEEP,cq0,256*2);
			PRFM_XI(PSTL1KEEP,cq1,256*2);
			CZERO_42;

			PRFM_XI(PLDL1KEEP,aa,0);
			PRFM_XI(PLDL1KEEP,aa,256*1);
			PRFM_XI(PLDL1KEEP,aa,256*2);
			PRFM_XI(PLDL1KEEP,aa,256*3);
			PRFM_XI(PLDL1KEEP,aa,256*4);
			PRFM_XI(PLDL1KEEP,aa,256*5);
			PRFM_XI(PLDL1KEEP,aa,256*6);
			PRFM_XI(PLDL1KEEP,aa,256*7);
			PRFM_XI(PLDL1KEEP,aa,256*8);

			bq0=bs0;
			bq1=bs1;
			ik=0;
			if(k>=2) { 

				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);

				LD1RW_ZXI(z27,p0,bq0,SIZE);
				LD1RW_ZXI(z28,p0,bq1,SIZE);

				LD1W_ZXI(z0,p0,aa,0);

				ik+=2;

				for(;ik<k-1;) {
                                        PRFM_XI(PLDL1KEEP,bq0,256);
                                        PRFM_XI(PLDL1KEEP,bq1,256);

					LD1W_ZXI(z1,p0,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p0,z0,z22);

					FMLA_ZSP( z3,p0,z0,z23);

					FMLA_ZSP( z7,p0,z1,z22);

					LD1W_ZXI(z0,p0,aa,2);

					FMLA_ZSP( z8,p0,z1,z23);

					LD1W_ZXI(z1,p0,aa,3);
					bq0+=2;

					FMLA_ZSP( z12,p0,z0,z22);
					FMLA_ZSP( z13,p0,z0,z23);

					bq1+=2;

					FMLA_ZSP( z17,p0,z1,z22);

					LD1W_ZXI(z0,p0,aa,4);

					FMLA_ZSP( z18,p0,z1,z23);

					LD1RW_ZXI(z22,p0,bq0,0);
					LD1RW_ZXI(z23,p0,bq1,0);

					LD1W_ZXI(z1,p0,aa,5);
					PRFM_XI(PLDL1KEEP,aa,256*10);

					FMLA_ZSP( z2,p0,z0,z27);
					FMLA_ZSP( z3,p0,z0,z28);

					FMLA_ZSP( z7,p0,z1,z27);

					LD1W_ZXI(z0,p0,aa,6);

					FMLA_ZSP( z8,p0,z1,z28);

					LD1W_ZXI(z1,p0,aa,7);
					aa+=nul*2;

					FMLA_ZSP( z12,p0,z0,z27);
					FMLA_ZSP( z13,p0,z0,z28);

					ik+=2;

					FMLA_ZSP( z17,p0,z1,z27);

					LD1W_ZXI(z0,p0,aa,0);
					NOP;

					FMLA_ZSP( z18,p0,z1,z28);

					LD1RW_ZXI(z27,p0,bq0,SIZE);
					LD1RW_ZXI(z28,p0,bq1,SIZE);
				}
				LD1W_ZXI(z1,p0,aa,1);

				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z3,p0,z0,z23);

				LD1W_ZXI(z0,p0,aa,2);

				FMLA_ZSP( z7 ,p0,z1,z22);
				FMLA_ZSP( z8 ,p0,z1,z23);

				LD1W_ZXI(z1,p0,aa,3);

				FMLA_ZSP( z12,p0,z0,z22);
				FMLA_ZSP( z13,p0,z0,z23);

				LD1W_ZXI(z0,p0,aa,4);

				FMLA_ZSP( z17,p0,z1,z22);
				FMLA_ZSP( z18,p0,z1,z23);

				LD1W_ZXI(z1,p0,aa,5);

				FMLA_ZSP( z2,p0,z0,z27);
				FMLA_ZSP( z3,p0,z0,z28);

				LD1W_ZXI(z0,p0,aa,6);

				FMLA_ZSP( z7,p0,z1,z27);
				FMLA_ZSP( z8,p0,z1,z28);

				LD1W_ZXI(z1,p0,aa,7);
				bq0+=2;
				bq1+=2;
				aa+=nul*2;

				if(k-ik>=1) {
					LD1RW_ZXI(z22,p0,bq0,0);
					LD1RW_ZXI(z23,p0,bq1,0);

					FMLA_ZSP( z12,p0,z0,z27);
					FMLA_ZSP( z13,p0,z0,z28);

					LD1W_ZXI(z0,p0,aa,0);

					FMLA_ZSP( z17,p0,z1,z27);
					FMLA_ZSP( z18,p0,z1,z28);

					LD1W_ZXI(z1,p0,aa,1);
				}
				else {
					FMLA_ZSP( z12,p0,z0,z27);
					FMLA_ZSP( z13,p0,z0,z28);
					FMLA_ZSP( z17,p0,z1,z27);
					FMLA_ZSP( z18,p0,z1,z28);
				}
			}
			else if(k==1) {
				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);

				LD1W_ZXI(z0,p0,aa,0);
				LD1W_ZXI(z1,p0,aa,1);

			}

			if(k-ik>=1) {

				FMLA_ZSP( z2, p0,z0,z22);
				FMLA_ZSP( z3, p0,z0,z23);

				LD1W_ZXI(z0,p0,aa,2);

				FMLA_ZSP( z7, p0,z1,z22);
				FMLA_ZSP( z8, p0,z1,z23);

				LD1W_ZXI(z1,p0,aa,3);

				FMLA_ZSP( z12,p0,z0,z22);
				FMLA_ZSP( z13,p0,z0,z23);
				FMLA_ZSP( z17,p0,z1,z22);
				FMLA_ZSP( z18,p0,z1,z23);

				aa+=nul;
				bq0+=1;
				bq1+=1;
				ik++;
			}
			LD1RW_ZXI(z0,p0,alphap,0);
			CALPHA_42;
			cq0+=nul;
			cq1+=nul;
		}

		if(m-im>nvl*3) {
			mm=m-im;
			im+=nvl*3;
			WHILELT_PSX(p1,im,m);
			mm=m-im;
			im=m;
			PRFM_XI(PSTL1KEEP,cq0,256*2);
			PRFM_XI(PSTL1KEEP,cq1,256*2);

			CZERO_42;

			PRFM_XI(PLDL1KEEP,aa,0);
			PRFM_XI(PLDL1KEEP,aa,256*1);
			PRFM_XI(PLDL1KEEP,aa,256*2);
			PRFM_XI(PLDL1KEEP,aa,256*3);
			PRFM_XI(PLDL1KEEP,aa,256*4);
			PRFM_XI(PLDL1KEEP,aa,256*5);
			PRFM_XI(PLDL1KEEP,aa,256*6);
			PRFM_XI(PLDL1KEEP,aa,256*7);
			PRFM_XI(PLDL1KEEP,aa,256*8);

			bq0=bs0;
			bq1=bs1;
			ik=0;
			if(k>=2) { 

				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);

				LD1RW_ZXI(z27,p0,bq0,SIZE);
				LD1RW_ZXI(z28,p0,bq1,SIZE);

				LD1W_ZXI(z0,p0,aa,0);

				ik+=2;

				for(;ik<k-1;) {
                                        PRFM_XI(PLDL1KEEP,bq0,256);
                                        PRFM_XI(PLDL1KEEP,bq1,256);

					LD1W_ZXI(z1,p0,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p0,z0,z22);
					FMLA_ZSP( z3,p0,z0,z23);
					FMLA_ZSP( z7,p0,z1,z22);

					LD1W_ZXI(z0,p0,aa,2);

					FMLA_ZSP( z8,p0,z1,z23);

					LD1W_ZXI(z1,p1,aa,3);
					aa+=nvl*3+mm;
					bq0+=2;

					FMLA_ZSP( z12,p0,z0,z22);
					FMLA_ZSP( z13,p0,z0,z23);

					bq1+=2;

					FMLA_ZSP( z17,p1,z1,z22);

					LD1W_ZXI(z0,p0,aa,0);

					FMLA_ZSP( z18,p1,z1,z23);

					LD1RW_ZXI(z22,p0,bq0,0);
					LD1RW_ZXI(z23,p0,bq1,0);

					LD1W_ZXI(z1,p0,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p0,z0,z27);
					FMLA_ZSP( z3,p0,z0,z28);

					FMLA_ZSP( z7,p0,z1,z27);

					LD1W_ZXI(z0,p0,aa,2);

					FMLA_ZSP( z8,p0,z1,z28);

					LD1W_ZXI(z1,p1,aa,3);
					aa+=nvl*3+mm;

					FMLA_ZSP( z12,p0,z0,z27);
					FMLA_ZSP( z13,p0,z0,z28);

					ik+=2;

					FMLA_ZSP( z17,p1,z1,z27);

					LD1W_ZXI(z0,p0,aa,0);
					NOP;

					FMLA_ZSP( z18,p1,z1,z28);

					LD1RW_ZXI(z27,p0,bq0,SIZE);
					LD1RW_ZXI(z28,p0,bq1,SIZE);
				}
				LD1W_ZXI(z1,p0,aa,1);

				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z3,p0,z0,z23);

				LD1W_ZXI(z0,p0,aa,2);

				FMLA_ZSP( z7,p0,z1,z22);
				FMLA_ZSP( z8,p0,z1,z23);

				LD1W_ZXI(z1,p1,aa,3);
				aa+=nvl*3+mm;

				FMLA_ZSP( z12,p0,z0,z22);
				FMLA_ZSP( z13,p0,z0,z23);

				LD1W_ZXI(z0,p0,aa,0);

				FMLA_ZSP( z17,p1,z1,z22);
				FMLA_ZSP( z18,p1,z1,z23);

				LD1W_ZXI(z1,p0,aa,1);

				FMLA_ZSP( z2,p0,z0,z27);
				FMLA_ZSP( z3,p0,z0,z28);

				LD1W_ZXI(z0,p0,aa,2);

				FMLA_ZSP( z7,p0,z1,z27);
				FMLA_ZSP( z8,p0,z1,z28);

				LD1W_ZXI(z1,p1,aa,3);

				FMLA_ZSP( z12,p0,z0,z27);
				FMLA_ZSP( z13,p0,z0,z28);
				FMLA_ZSP( z17,p1,z1,z27);
				FMLA_ZSP( z18,p1,z1,z28);

				bq0+=2;
				bq1+=2;
				aa+=nvl*3+mm;

			}
			if(k-ik>=1) {
				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);

				LD1W_ZXI(z0,p0,aa,0);
				LD1W_ZXI(z1,p0,aa,1);

				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z3,p0,z0,z23);
				FMLA_ZSP( z7,p0,z1,z22);
				FMLA_ZSP( z8,p0,z1,z23);

				LD1W_ZXI(z0,p0,aa,2);
				LD1W_ZXI(z1,p1,aa,3);

				FMLA_ZSP( z12,p0,z0,z22);
				FMLA_ZSP( z13,p0,z0,z23);
				FMLA_ZSP( z17,p1,z1,z22);
				FMLA_ZSP( z18,p1,z1,z23);

				aa+=nvl*3+mm;
				bq0+=1;
				bq1+=1;
			}
			LD1RW_ZXI(z0,p0,alphap,0);
			CALPHA_P_42;
			cq0+=nvl*3+mm;
			cq1+=nvl*3+mm;
		}

		else if(m-im>nvl*2) {
			mm=m-im;
			im+=nvl*2;
			WHILELT_PSX(p1,im,m);
			mm=m-im;
			im=m;
			PRFM_XI(PSTL1KEEP,cq0,256*2);
			PRFM_XI(PSTL1KEEP,cq1,256*2);

			CZERO_32;

			PRFM_XI(PLDL1KEEP,aa,0);
			PRFM_XI(PLDL1KEEP,aa,256*1);
			PRFM_XI(PLDL1KEEP,aa,256*2);
			PRFM_XI(PLDL1KEEP,aa,256*3);
			PRFM_XI(PLDL1KEEP,aa,256*4);
			PRFM_XI(PLDL1KEEP,aa,256*5);
			PRFM_XI(PLDL1KEEP,aa,256*6);
			PRFM_XI(PLDL1KEEP,aa,256*7);
			PRFM_XI(PLDL1KEEP,aa,256*8);

			bq0=bs0;
			bq1=bs1;
			ik=0;
			if(k>=2) { 

				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);

				LD1RW_ZXI(z27,p0,bq0,SIZE);
				LD1RW_ZXI(z28,p0,bq1,SIZE);

				LD1W_ZXI(z0,p0,aa,0);

				ik+=2;

				for(;ik<k-1;) {
                                        PRFM_XI(PLDL1KEEP,bq0,256);
                                        PRFM_XI(PLDL1KEEP,bq1,256);

					LD1W_ZXI(z1,p0,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p0,z0,z22);
					FMLA_ZSP( z3,p0,z0,z23);

					FMLA_ZSP( z7,p0,z1,z22);

					LD1W_ZXI(z17,p1,aa,2);

					FMLA_ZSP( z8,p0,z1,z23);

					aa+=nvl*2+mm;
					bq0+=2;

					FMLA_ZSP( z12,p1,z17,z22);
					FMLA_ZSP( z13,p1,z17,z23);

					bq1+=2;

					LD1W_ZXI(z0,p0,aa,0);

					LD1RW_ZXI(z22,p0,bq0,0);
					LD1RW_ZXI(z23,p0,bq1,0);

					LD1W_ZXI(z1,p0,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p0,z0,z27);
					FMLA_ZSP( z3,p0,z0,z28);

					FMLA_ZSP( z7,p0,z1,z27);

					LD1W_ZXI(z17,p1,aa,2);

					FMLA_ZSP( z8,p0,z1,z28);

					aa+=nvl*2+mm;

					FMLA_ZSP( z12,p1,z17,z27);
					FMLA_ZSP( z13,p1,z17,z28);

					ik+=2;

					LD1W_ZXI(z0,p0,aa,0);

					LD1RW_ZXI(z27,p0,bq0,SIZE);
					LD1RW_ZXI(z28,p0,bq1,SIZE);

				}
				LD1W_ZXI(z1,p0,aa,1);

				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z3,p0,z0,z23);

				LD1W_ZXI(z17,p1,aa,2);

				FMLA_ZSP( z7,p0,z1,z22);
				FMLA_ZSP( z8,p0,z1,z23);

				aa+=nvl*2+mm;

				FMLA_ZSP( z12,p1,z17,z22);
				FMLA_ZSP( z13,p1,z17,z23);

				LD1W_ZXI(z0,p0,aa,0);

				LD1W_ZXI(z1,p0,aa,1);

				FMLA_ZSP( z2,p0,z0,z27);
				FMLA_ZSP( z3,p0,z0,z28);

				LD1W_ZXI(z17,p1,aa,2);

				FMLA_ZSP( z7,p0,z1,z27);
				FMLA_ZSP( z8,p0,z1,z28);

				FMLA_ZSP( z12,p1,z17,z27);
				FMLA_ZSP( z13,p1,z17,z28);

				bq0+=2;
				bq1+=2;
				aa+=nvl*2+mm;
			}
			if(k-ik>=1) {
				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);

				LD1W_ZXI(z0,p0,aa,0);
				LD1W_ZXI(z1,p0,aa,1);
				LD1W_ZXI(z17,p1,aa,2);

				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z3,p0,z0,z23);

				FMLA_ZSP( z7,p0,z1,z22);
				FMLA_ZSP( z8,p0,z1,z23);

				FMLA_ZSP( z12,p1,z17,z22);
				FMLA_ZSP( z13,p1,z17,z23);

				aa+=nvl*2+mm;
				bq0+=1;
				bq1+=1;
			}
			LD1RW_ZXI(z0,p0,alphap,0);
			CALPHA_P_32;
			cq0+=nvl*2+mm;
			cq1+=nvl*2+mm;
		}

		else if(m-im>nvl*1) {
			mm=m-im;
			im+=nvl*1;
			WHILELT_PSX(p1,im,m);
			mm=m-im;
			im=m;
			PRFM_XI(PSTL1KEEP,cq0,256*2);
			PRFM_XI(PSTL1KEEP,cq1,256*2);

			CZERO_22;

			PRFM_XI(PLDL1KEEP,aa,0);
			PRFM_XI(PLDL1KEEP,aa,256*1);
			PRFM_XI(PLDL1KEEP,aa,256*2);
			PRFM_XI(PLDL1KEEP,aa,256*3);
			PRFM_XI(PLDL1KEEP,aa,256*4);
			PRFM_XI(PLDL1KEEP,aa,256*5);
			PRFM_XI(PLDL1KEEP,aa,256*6);
			PRFM_XI(PLDL1KEEP,aa,256*7);
			PRFM_XI(PLDL1KEEP,aa,256*8);

			bq0=bs0;
			bq1=bs1;
			ik=0;
			if(k>=2) { 

				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);

				LD1RW_ZXI(z27,p0,bq0,SIZE);
				LD1RW_ZXI(z28,p0,bq1,SIZE);

				LD1W_ZXI(z0,p0,aa,0);

				ik+=2;

				for(;ik<k-1;) {
                                        PRFM_XI(PLDL1KEEP,bq0,256);
                                        PRFM_XI(PLDL1KEEP,bq1,256);

					LD1W_ZXI(z1,p1,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p0,z0,z22);
					FMLA_ZSP( z3,p0,z0,z23);

					FMLA_ZSP( z7,p1,z1,z22);
					FMLA_ZSP( z8,p1,z1,z23);

					aa+=nvl*1+mm;
					bq0+=2;
					bq1+=2;

					LD1W_ZXI(z0,p0,aa,0);

					LD1RW_ZXI(z22,p0,bq0,0);
					LD1RW_ZXI(z23,p0,bq1,0);

					LD1W_ZXI(z1,p1,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);

					FMLA_ZSP( z2,p0,z0,z27);
					FMLA_ZSP( z3,p0,z0,z28);

					FMLA_ZSP( z7,p1,z1,z27);
					FMLA_ZSP( z8,p1,z1,z28);

					aa+=nvl*1+mm;
					ik+=2;

					LD1W_ZXI(z0,p0,aa,0);
					NOP;

					LD1RW_ZXI(z27,p0,bq0,SIZE);
					LD1RW_ZXI(z28,p0,bq1,SIZE);

				}
				LD1W_ZXI(z1,p1,aa,1);

				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z3,p0,z0,z23);

				FMLA_ZSP( z7,p1,z1,z22);
				FMLA_ZSP( z8,p1,z1,z23);

				aa+=nvl*1+mm;

				LD1W_ZXI(z0,p0,aa,0);
				LD1W_ZXI(z1,p1,aa,1);

				FMLA_ZSP( z2,p0,z0,z27);
				FMLA_ZSP( z3,p0,z0,z28);

				FMLA_ZSP( z7,p1,z1,z27);
				FMLA_ZSP( z8,p1,z1,z28);

				bq0+=2;
				bq1+=2;
				aa+=nvl*1+mm;

			}
			if(k-ik>=1) {
				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);

				LD1W_ZXI(z0,p0,aa,0);
				LD1W_ZXI(z1,p1,aa,1);

				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z3,p0,z0,z23);

				FMLA_ZSP( z7,p1,z1,z22);
				FMLA_ZSP( z8,p1,z1,z23);

				aa+=nvl*1+mm;
				bq0+=1;
				bq1+=1;
			}
			LD1RW_ZXI(z0,p0,alphap,0);
			CALPHA_P_22;
			cq0+=nvl*1+mm;
			cq1+=nvl*1+mm;
		}
		else if(m-im>0) {
			mm=m-im;
			WHILELT_PSX(p1,im,m);
			mm=m-im;
			im=m;
			PRFM_XI(PSTL1KEEP,cq0,256*2);
			PRFM_XI(PSTL1KEEP,cq1,256*2);

			CZERO_12;

			PRFM_XI(PLDL1KEEP,aa,0);
			PRFM_XI(PLDL1KEEP,aa,256*1);
			PRFM_XI(PLDL1KEEP,aa,256*2);
			PRFM_XI(PLDL1KEEP,aa,256*3);
			PRFM_XI(PLDL1KEEP,aa,256*4);
			PRFM_XI(PLDL1KEEP,aa,256*5);
			PRFM_XI(PLDL1KEEP,aa,256*6);
			PRFM_XI(PLDL1KEEP,aa,256*7);
			PRFM_XI(PLDL1KEEP,aa,256*8);

			bq0=bs0;
			bq1=bs1;
			ik=0;
			if(k>=2) { 
				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);
				LD1RW_ZXI(z27,p0,bq0,SIZE);
				LD1RW_ZXI(z28,p0,bq1,SIZE);
				LD1W_ZXI(z0,p1,aa,0);
				ik+=2;
				for(;ik<k-1;) {
                                        PRFM_XI(PLDL1KEEP,bq0,256);
                                        PRFM_XI(PLDL1KEEP,bq1,256);

					PRFM_XI(PLDL1KEEP,aa,256*9);
					FMLA_ZSP( z2,p1,z0,z22);
					FMLA_ZSP( z3,p1,z0,z23);
					aa+=mm;
					bq0+=2;
					bq1+=2;
					LD1W_ZXI(z0,p1,aa,0);
					LD1RW_ZXI(z22,p0,bq0,0);
					LD1RW_ZXI(z23,p0,bq1,0);
					PRFM_XI(PLDL1KEEP,aa,256*9);
					FMLA_ZSP( z2,p1,z0,z27);
					FMLA_ZSP( z3,p1,z0,z28);
					aa+=mm;
					ik+=2;
					LD1W_ZXI(z0,p1,aa,0);
					LD1RW_ZXI(z27,p0,bq0,SIZE);
					LD1RW_ZXI(z28,p0,bq1,SIZE);
				}
				FMLA_ZSP( z2,p1,z0,z22);
				FMLA_ZSP( z3,p1,z0,z23);
				aa+=mm;
				LD1W_ZXI(z0,p1,aa,0);
				FMLA_ZSP( z2,p1,z0,z27);
				FMLA_ZSP( z3,p1,z0,z28);
				bq0+=2;
				bq1+=2;
				aa+=mm;
			}
			if(k-ik>=1) {
				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z23,p0,bq1,0);
				LD1W_ZXI(z0,p1,aa,0);
				FMLA_ZSP( z2,p1,z0,z22);
				FMLA_ZSP( z3,p1,z0,z23);
				aa+=mm;
				bq0+=1;
				bq1+=1;
			}
			LD1RW_ZXI(z0,p0,alphap,0);
			CALPHA_P_12;
			cq0+=mm;
			cq1+=mm;
		}
		cq0+=ldc*2-m;
		cq1+=ldc*2-m;
		bs0+=ldb*2;
		bs1+=ldb*2;
	}
	else if(n-in==1) {
		in+=1;
		aa=a;
		im=0;
		for(;im<m-nul+1;im+=nul) {
			PRFM_XI(PSTL1KEEP,cq0,256*2);
			CZERO_41;

			PRFM_XI(PLDL1KEEP,aa,0);
			PRFM_XI(PLDL1KEEP,aa,256*1);
			PRFM_XI(PLDL1KEEP,aa,256*2);
			PRFM_XI(PLDL1KEEP,aa,256*3);
			PRFM_XI(PLDL1KEEP,aa,256*4);
			PRFM_XI(PLDL1KEEP,aa,256*5);
			PRFM_XI(PLDL1KEEP,aa,256*6);
			PRFM_XI(PLDL1KEEP,aa,256*7);
			PRFM_XI(PLDL1KEEP,aa,256*8);

			bq0=bs0;
			ik=0;
			if(k>=2) { 
				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z27,p0,bq0,SIZE);
				LD1W_ZXI(z0,p0,aa,0);
				ik+=2;
				for(;ik<k-1;) {
                                        PRFM_XI(PLDL1KEEP,bq0,256);

					LD1W_ZXI(z1,p0,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);
					FMLA_ZSP( z2,p0,z0,z22);
					FMLA_ZSP( z7,p0,z1,z22);
					LD1W_ZXI(z0,p0,aa,2);
					LD1W_ZXI(z1,p0,aa,3);
					bq0+=2;
					FMLA_ZSP( z12,p0,z0,z22);
					FMLA_ZSP( z17,p0,z1,z22);
					LD1W_ZXI(z0,p0,aa,4);
					LD1RW_ZXI(z22,p0,bq0,0);
					LD1W_ZXI(z1,p0,aa,5);
					PRFM_XI(PLDL1KEEP,aa,256*10);
					FMLA_ZSP( z2,p0,z0,z27);
					FMLA_ZSP( z7,p0,z1,z27);
					LD1W_ZXI(z0,p0,aa,6);
					LD1W_ZXI(z1,p0,aa,7);
					aa+=nul*2;
					FMLA_ZSP( z12,p0,z0,z27);
					ik+=2;
					FMLA_ZSP( z17,p0,z1,z27);
					LD1W_ZXI(z0,p0,aa,0);
					LD1RW_ZXI(z27,p0,bq0,SIZE);
				}
				LD1W_ZXI(z1,p0,aa,1);
				FMLA_ZSP( z2,p0,z0,z22);
				LD1W_ZXI(z0,p0,aa,2);
				FMLA_ZSP( z7 ,p0,z1,z22);
				LD1W_ZXI(z1,p0,aa,3);
				FMLA_ZSP( z12,p0,z0,z22);
				LD1W_ZXI(z0,p0,aa,4);
				FMLA_ZSP( z17,p0,z1,z22);
				LD1W_ZXI(z1,p0,aa,5);
				FMLA_ZSP( z2,p0,z0,z27);
				LD1W_ZXI(z0,p0,aa,6);
				FMLA_ZSP( z7,p0,z1,z27);
				LD1W_ZXI(z1,p0,aa,7);
				bq0+=2;
				aa+=nul*2;

				if(k-ik>=1) {
					LD1RW_ZXI(z22,p0,bq0,0);
					FMLA_ZSP( z12,p0,z0,z27);
					LD1W_ZXI(z0,p0,aa,0);
					FMLA_ZSP( z17,p0,z1,z27);
					LD1W_ZXI(z1,p0,aa,1);
				}
				else {
					FMLA_ZSP( z12,p0,z0,z27);
					FMLA_ZSP( z17,p0,z1,z27);
				}
			}
			else if(k==1) {
				LD1RW_ZXI(z22,p0,bq0,0);
				LD1W_ZXI(z0,p0,aa,0);
				LD1W_ZXI(z1,p0,aa,1);
			}

			if(k-ik>=1) {
				FMLA_ZSP( z2, p0,z0,z22);
				LD1W_ZXI(z0,p0,aa,2);
				FMLA_ZSP( z7, p0,z1,z22);
				LD1W_ZXI(z1,p0,aa,3);
				FMLA_ZSP( z12,p0,z0,z22);
				FMLA_ZSP( z17,p0,z1,z22);
				aa+=nul;
				bq0+=1;
				ik++;
			}
			LD1RW_ZXI(z0,p0,alphap,0);
			CALPHA_41;

			cq0+=nul;
		}
		if(m-im>nvl*3) {
			mm=m-im;
			im+=nvl*3;
			WHILELT_PSX(p1,im,m);
			mm=m-im;
			im=m;
			PRFM_XI(PSTL1KEEP,cq0,256*2);

			CZERO_41;

			PRFM_XI(PLDL1KEEP,aa,0);
			PRFM_XI(PLDL1KEEP,aa,256*1);
			PRFM_XI(PLDL1KEEP,aa,256*2);
			PRFM_XI(PLDL1KEEP,aa,256*3);
			PRFM_XI(PLDL1KEEP,aa,256*4);
			PRFM_XI(PLDL1KEEP,aa,256*5);
			PRFM_XI(PLDL1KEEP,aa,256*6);
			PRFM_XI(PLDL1KEEP,aa,256*7);
			PRFM_XI(PLDL1KEEP,aa,256*8);

			bq0=bs0;
			ik=0;
			if(k>=2) { 
				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z27,p0,bq0,SIZE);
				LD1W_ZXI(z0,p0,aa,0);
				ik+=2;
				for(;ik<k-1;) {
                                        PRFM_XI(PLDL1KEEP,bq0,256);

					LD1W_ZXI(z1,p0,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);
					FMLA_ZSP( z2,p0,z0,z22);
					FMLA_ZSP( z7,p0,z1,z22);
					LD1W_ZXI(z0,p0,aa,2);
					LD1W_ZXI(z1,p1,aa,3);
					aa+=nvl*3+mm;
					bq0+=2;
					FMLA_ZSP( z12,p0,z0,z22);
					FMLA_ZSP( z17,p1,z1,z22);
					LD1W_ZXI(z0,p0,aa,0);
					LD1RW_ZXI(z22,p0,bq0,0);
					LD1W_ZXI(z1,p0,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);
					FMLA_ZSP( z2,p0,z0,z27);
					FMLA_ZSP( z7,p0,z1,z27);
					LD1W_ZXI(z0,p0,aa,2);
					LD1W_ZXI(z1,p1,aa,3);
					aa+=nvl*3+mm;
					FMLA_ZSP( z12,p0,z0,z27);
					ik+=2;
					FMLA_ZSP( z17,p1,z1,z27);
					LD1W_ZXI(z0,p0,aa,0);
					LD1RW_ZXI(z27,p0,bq0,SIZE);
				}
				LD1W_ZXI(z1,p0,aa,1);
				FMLA_ZSP( z2,p0,z0,z22);
				LD1W_ZXI(z0,p0,aa,2);
				FMLA_ZSP( z7,p0,z1,z22);
				LD1W_ZXI(z1,p1,aa,3);
				aa+=nvl*3+mm;
				FMLA_ZSP( z12,p0,z0,z22);
				LD1W_ZXI(z0,p0,aa,0);
				FMLA_ZSP( z17,p1,z1,z22);
				LD1W_ZXI(z1,p0,aa,1);
				FMLA_ZSP( z2,p0,z0,z27);
				LD1W_ZXI(z0,p0,aa,2);
				FMLA_ZSP( z7,p0,z1,z27);
				LD1W_ZXI(z1,p1,aa,3);
				FMLA_ZSP( z12,p0,z0,z27);
				FMLA_ZSP( z17,p1,z1,z27);
				bq0+=2;
				aa+=nvl*3+mm;
			}
			if(k-ik>=1) {
				LD1RW_ZXI(z22,p0,bq0,0);
				LD1W_ZXI(z0,p0,aa,0);
				LD1W_ZXI(z1,p0,aa,1);
				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z7,p0,z1,z22);
				LD1W_ZXI(z0,p0,aa,2);
				LD1W_ZXI(z1,p1,aa,3);
				FMLA_ZSP( z12,p0,z0,z22);
				FMLA_ZSP( z17,p1,z1,z22);
				aa+=nvl*3+mm;
				bq0+=1;
			}
			LD1RW_ZXI(z0,p0,alphap,0);
			CALPHA_P_41;
			cq0+=nvl*3+mm;
		}

		else if(m-im>nvl*2) {
			mm=m-im;
			im+=nvl*2;
			WHILELT_PSX(p1,im,m);
			mm=m-im;
			im=m;
			PRFM_XI(PSTL1KEEP,cq0,256*2);

			CZERO_31;

			PRFM_XI(PLDL1KEEP,aa,0);
			PRFM_XI(PLDL1KEEP,aa,256*1);
			PRFM_XI(PLDL1KEEP,aa,256*2);
			PRFM_XI(PLDL1KEEP,aa,256*3);
			PRFM_XI(PLDL1KEEP,aa,256*4);
			PRFM_XI(PLDL1KEEP,aa,256*5);
			PRFM_XI(PLDL1KEEP,aa,256*6);
			PRFM_XI(PLDL1KEEP,aa,256*7);
			PRFM_XI(PLDL1KEEP,aa,256*8);

			bq0=bs0;
			ik=0;
			if(k>=2) { 
				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z27,p0,bq0,SIZE);
				LD1W_ZXI(z0,p0,aa,0);
				ik+=2;
				for(;ik<k-1;) {
                                        PRFM_XI(PLDL1KEEP,bq0,256);

					LD1W_ZXI(z1,p0,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);
					FMLA_ZSP( z2,p0,z0,z22);
					FMLA_ZSP( z7,p0,z1,z22);
					LD1W_ZXI(z17,p1,aa,2);
					aa+=nvl*2+mm;
					bq0+=2;
					FMLA_ZSP( z12,p1,z17,z22);
					LD1W_ZXI(z0,p0,aa,0);
					LD1RW_ZXI(z22,p0,bq0,0);
					LD1W_ZXI(z1,p0,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);
					FMLA_ZSP( z2,p0,z0,z27);
					FMLA_ZSP( z7,p0,z1,z27);
					LD1W_ZXI(z17,p1,aa,2);
					aa+=nvl*2+mm;
					FMLA_ZSP( z12,p1,z17,z27);
					ik+=2;
					LD1W_ZXI(z0,p0,aa,0);
					LD1RW_ZXI(z27,p0,bq0,SIZE);
				}
				LD1W_ZXI(z1,p0,aa,1);
				FMLA_ZSP( z2,p0,z0,z22);
				LD1W_ZXI(z17,p1,aa,2);
				FMLA_ZSP( z7,p0,z1,z22);
				aa+=nvl*2+mm;
				FMLA_ZSP( z12,p1,z17,z22);
				LD1W_ZXI(z0,p0,aa,0);
				LD1W_ZXI(z1,p0,aa,1);
				FMLA_ZSP( z2,p0,z0,z27);
				LD1W_ZXI(z17,p1,aa,2);
				FMLA_ZSP( z7,p0,z1,z27);
				FMLA_ZSP( z12,p1,z17,z27);
				bq0+=2;
				aa+=nvl*2+mm;
			}
			if(k-ik>=1) {
				LD1RW_ZXI(z22,p0,bq0,0);

				LD1W_ZXI(z0,p0,aa,0);
				LD1W_ZXI(z1,p0,aa,1);
				LD1W_ZXI(z17,p1,aa,2);

				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z7,p0,z1,z22);
				FMLA_ZSP( z12,p1,z17,z22);

				aa+=nvl*2+mm;
				bq0+=1;
			}
			LD1RW_ZXI(z0,p0,alphap,0);
			CALPHA_P_31;
			cq0+=nvl*2+mm;
		}

		else if(m-im>nvl*1) {
			mm=m-im;
			im+=nvl*1;
			WHILELT_PSX(p1,im,m);
			mm=m-im;
			im=m;
			PRFM_XI(PSTL1KEEP,cq0,256*2);

			CZERO_21;

			PRFM_XI(PLDL1KEEP,aa,0);
			PRFM_XI(PLDL1KEEP,aa,256*1);
			PRFM_XI(PLDL1KEEP,aa,256*2);
			PRFM_XI(PLDL1KEEP,aa,256*3);
			PRFM_XI(PLDL1KEEP,aa,256*4);
			PRFM_XI(PLDL1KEEP,aa,256*5);
			PRFM_XI(PLDL1KEEP,aa,256*6);
			PRFM_XI(PLDL1KEEP,aa,256*7);
			PRFM_XI(PLDL1KEEP,aa,256*8);

			bq0=bs0;
			ik=0;
			if(k>=2) { 
				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z27,p0,bq0,SIZE);
				LD1W_ZXI(z0,p0,aa,0);
				ik+=2;
				for(;ik<k-1;) {
                                        PRFM_XI(PLDL1KEEP,bq0,256);

					LD1W_ZXI(z1,p1,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);
					FMLA_ZSP( z2,p0,z0,z22);
					FMLA_ZSP( z7,p1,z1,z22);
					aa+=nvl*1+mm;
					bq0+=2;
					LD1W_ZXI(z0,p0,aa,0);
					LD1RW_ZXI(z22,p0,bq0,0);
					LD1W_ZXI(z1,p1,aa,1);
					PRFM_XI(PLDL1KEEP,aa,256*9);
					FMLA_ZSP( z2,p0,z0,z27);
					FMLA_ZSP( z7,p1,z1,z27);
					aa+=nvl*1+mm;
					ik+=2;
					LD1W_ZXI(z0,p0,aa,0);
					LD1RW_ZXI(z27,p0,bq0,SIZE);
				}
				LD1W_ZXI(z1,p1,aa,1);
				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z7,p1,z1,z22);
				aa+=nvl*1+mm;
				LD1W_ZXI(z0,p0,aa,0);
				LD1W_ZXI(z1,p1,aa,1);
				FMLA_ZSP( z2,p0,z0,z27);
				FMLA_ZSP( z7,p1,z1,z27);
				bq0+=2;
				aa+=nvl*1+mm;
			}
			if(k-ik>=1) {
				LD1RW_ZXI(z22,p0,bq0,0);
				LD1W_ZXI(z0,p0,aa,0);
				LD1W_ZXI(z1,p1,aa,1);
				FMLA_ZSP( z2,p0,z0,z22);
				FMLA_ZSP( z7,p1,z1,z22);
				aa+=nvl*1+mm;
				bq0+=1;
			}
			LD1RW_ZXI(z0,p0,alphap,0);
			CALPHA_P_21;
			cq0+=nvl*1+mm;
		}
		else if(m-im>0) {
			mm=m-im;
			WHILELT_PSX(p1,im,m);
			mm=m-im;
			im=m;
			PRFM_XI(PSTL1KEEP,cq0,256*2);

			CZERO_11;

			PRFM_XI(PLDL1KEEP,aa,0);
			PRFM_XI(PLDL1KEEP,aa,256*1);
			PRFM_XI(PLDL1KEEP,aa,256*2);
			PRFM_XI(PLDL1KEEP,aa,256*3);
			PRFM_XI(PLDL1KEEP,aa,256*4);
			PRFM_XI(PLDL1KEEP,aa,256*5);
			PRFM_XI(PLDL1KEEP,aa,256*6);
			PRFM_XI(PLDL1KEEP,aa,256*7);
			PRFM_XI(PLDL1KEEP,aa,256*8);
			bq0=bs0;
			ik=0;
			if(k>=2) { 
				LD1RW_ZXI(z22,p0,bq0,0);
				LD1RW_ZXI(z27,p0,bq0,SIZE);
				LD1W_ZXI(z0,p1,aa,0);
				ik+=2;
				for(;ik<k-1;) {
                                        PRFM_XI(PLDL1KEEP,bq0,256);

					PRFM_XI(PLDL1KEEP,aa,256*9);
					FMLA_ZSP( z2,p1,z0,z22);
					aa+=mm;
					bq0+=2;
					LD1W_ZXI(z1,p1,aa,0);
					LD1RW_ZXI(z22,p0,bq0,0);
					PRFM_XI(PLDL1KEEP,aa,256*9);
					FMLA_ZSP( z2,p1,z1,z27);
					aa+=mm;
					ik+=2;
					LD1W_ZXI(z0,p1,aa,0);
					LD1RW_ZXI(z27,p0,bq0,SIZE);
				}
				FMLA_ZSP( z2,p1,z0,z22);
				aa+=mm;
				LD1W_ZXI(z1,p1,aa,0);
				FMLA_ZSP( z2,p1,z1,z27);
				bq0+=2;
				aa+=mm;
			}
			if(k-ik>=1) {
				LD1RW_ZXI(z22,p0,bq0,0);
				LD1W_ZXI(z0,p1,aa,0);
				FMLA_ZSP( z2,p1,z0,z22);
				aa+=mm;
				bq0+=1;
			}
			LD1RW_ZXI(z0,p0,alphap,0);
			CALPHA_P_11;
			cq0+=mm;
		}
		cq0+=ldc*1-m;
		bs0+=ldb*1;
	}

	return 0;
}


