/***************************************************************************
(c) RIKEN 2024, 2024. All rights reserved. sneg_tcopy_sve_v4.c 0.3.26
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

int CNAME(BLASLONG k,BLASLONG m,FLOAT *a,BLASLONG lda,FLOAT *aw)
{
	FLOAT *ap0,*ap1,*awp,*awp2,*apre0,*apre1,*awpre;

	BLASLONG ik,j,im,ii;
	BLASLONG kb,iks,ike,kl,ipre;
	BLASLONG mm;
	BLASLONG nul;
	BLASLONG nvl=svcntw();

	PTRUE_PS(p0);

	nul=4*nvl;

	kb=k;
	iks=0;
	kl=k;
	ike=k-1;

	ipre=lda*2;
	if(lda<=m+32) {
		ipre=32;
	}

	if(kl>0) {
		ap0=a+lda*iks;
		ap1=ap0+lda;
		apre0=ap0+lda*8;
		ik=iks;

		for(;ik<=ike-1;ik+=2) {
			im=0;
			awp=aw+nul*ik;
			if(m-im>=nul) {
				im+=nul;
				apre0=ap0+lda*4;
				apre1=ap1+lda*4;
				PRFM_XI(PLDL2STRM,apre0,0);
				PRFM_XI(PLDL2STRM,apre1,0);
				PRFM_XI(PLDL2STRM,apre0,256);
				PRFM_XI(PLDL2STRM,apre1,256);

				apre0=ap0+lda*2;
				apre1=ap1+lda*2;
				PRFM_XI(PLDL1STRM,apre0,0);
				PRFM_XI(PLDL1STRM,apre1,0);
				PRFM_XI(PLDL1STRM,apre0,256);
				PRFM_XI(PLDL1STRM,apre1,256);

				LD1W_ZXI(z0,p0,ap0,0);
				LD1W_ZXI(z1,p0,ap0,1);
				LD1W_ZXI(z2,p0,ap0,2);
				LD1W_ZXI(z3,p0,ap0,3);

				LD1W_ZXI(z8 ,p0,ap1,0);
				LD1W_ZXI(z9 ,p0,ap1,1);
				LD1W_ZXI(z10,p0,ap1,2);
				LD1W_ZXI(z11,p0,ap1,3);

				FNEG_ZSP(z0,p0,z0);
				FNEG_ZSP(z1,p0,z1);
				FNEG_ZSP(z2,p0,z2);
				FNEG_ZSP(z3,p0,z3);

				FNEG_ZSP(z8 ,p0,z8);
				FNEG_ZSP(z9 ,p0,z9);
				FNEG_ZSP(z10,p0,z10);
				FNEG_ZSP(z11,p0,z11);

				for(;im<m-nul+1;im+=nul) {
					ap0+=nul;
					ap1+=nul;

					apre0=ap0+lda*4;
					apre1=ap1+lda*4;
					PRFM_XI(PLDL2STRM,apre0,0);
					PRFM_XI(PLDL2STRM,apre1,0);

					ST1W_ZXI(z0,p0,awp,0);
					ST1W_ZXI(z1,p0,awp,1);
					ST1W_ZXI(z2,p0,awp,2);
					ST1W_ZXI(z3,p0,awp,3);
					awp2=awp+nul;

					ST1W_ZXI(z8,p0,awp2,0);
					ST1W_ZXI(z9,p0,awp2,1);
					ST1W_ZXI(z10,p0,awp2,2);
					ST1W_ZXI(z11,p0,awp2,3);
					awp+=k*nul;

					LD1W_ZXI(z0,p0,ap0,0);
					LD1W_ZXI(z1,p0,ap0,1);
					LD1W_ZXI(z2,p0,ap0,2);
					LD1W_ZXI(z3,p0,ap0,3);

					LD1W_ZXI(z8 ,p0,ap1,0);
					LD1W_ZXI(z9 ,p0,ap1,1);
					LD1W_ZXI(z10,p0,ap1,2);
					LD1W_ZXI(z11,p0,ap1,3);

					FNEG_ZSP(z0,p0,z0);
					FNEG_ZSP(z1,p0,z1);
					FNEG_ZSP(z2,p0,z2);
					FNEG_ZSP(z3,p0,z3);

					FNEG_ZSP(z8 ,p0,z8);
					FNEG_ZSP(z9 ,p0,z9);
					FNEG_ZSP(z10,p0,z10);
					FNEG_ZSP(z11,p0,z11);

					apre0=ap0+ipre;
					apre1=ap1+ipre;
					PRFM_XI(PLDL1STRM,apre0,0);
					PRFM_XI(PLDL1STRM,apre1,0);

					awp2=awp+64;
					PRFM_XI(PLDL1STRM,awp2,0);
					PRFM_XI(PLDL1STRM,awp2,256);

				}

				ST1W_ZXI(z0,p0,awp,0);
				ST1W_ZXI(z1,p0,awp,1);
				ST1W_ZXI(z2,p0,awp,2);
				ST1W_ZXI(z3,p0,awp,3);

				awp2=awp+nul;
				ST1W_ZXI(z8,p0,awp2,0);
				ST1W_ZXI(z9,p0,awp2,1);
				ST1W_ZXI(z10,p0,awp2,2);
				ST1W_ZXI(z11,p0,awp2,3);

				awp+=k*nul;
				ap0+=nul;
				ap1+=nul;
			}
			mm=m-im;
			awp=aw+k*(m-mm)+ik*mm;
			awp2=awp+mm;

			if(mm>=1) {
				if(mm>nvl*3) {
					j=nvl*3;
					WHILELT_PSX(p1,j,mm);

                                        apre0=ap0+lda*4;
                                        apre1=ap1+lda*4;
                                        PRFM_XI(PLDL2STRM,apre0,0);
                                        PRFM_XI(PLDL2STRM,apre1,0);

                                        apre0=ap0+lda*2;
                                        apre1=ap1+lda*2;
                                        PRFM_XI(PLDL1STRM,apre0,0);
                                        PRFM_XI(PLDL1STRM,apre1,0);

					LD1W_ZXI(z0,p0,ap0,0);
					LD1W_ZXI(z1,p0,ap1,0);
					LD1W_ZXI(z2,p0,ap0,1);
					LD1W_ZXI(z3,p0,ap1,1);
					LD1W_ZXI(z4,p0,ap0,2);
					LD1W_ZXI(z5,p0,ap1,2);
					LD1W_ZXI(z6,p1,ap0,3);
					LD1W_ZXI(z7,p1,ap1,3);

					FNEG_ZSP(z0,p0,z0);
					FNEG_ZSP(z1,p0,z1);
					FNEG_ZSP(z2,p0,z2);
					FNEG_ZSP(z3,p0,z3);
					FNEG_ZSP(z4,p0,z4);
					FNEG_ZSP(z5,p0,z5);
					FNEG_ZSP(z6,p1,z6);
					FNEG_ZSP(z7,p1,z7);

					ST1W_ZXI(z0,p0,awp ,0);
					ST1W_ZXI(z1,p0,awp2,0);
					ST1W_ZXI(z2,p0,awp ,1);
					ST1W_ZXI(z3,p0,awp2,1);
					ST1W_ZXI(z4,p0,awp ,2);
					ST1W_ZXI(z5,p0,awp2,2);
					ST1W_ZXI(z6,p1,awp ,3);
					ST1W_ZXI(z7,p1,awp2,3);

					ap0+=mm;
					ap1+=mm;

				}
				else if(mm>nvl*2) {
					j=nvl*2;
					WHILELT_PSX(p1,j,mm);

                                        apre0=ap0+lda*4;
                                        apre1=ap1+lda*4;
                                        PRFM_XI(PLDL2STRM,apre0,0);
                                        PRFM_XI(PLDL2STRM,apre1,0);

					LD1W_ZXI(z0,p0,ap0,0);
					LD1W_ZXI(z1,p0,ap1,0);
					LD1W_ZXI(z2,p0,ap0,1);
					LD1W_ZXI(z3,p0,ap1,1);
					LD1W_ZXI(z4,p1,ap0,2);
					LD1W_ZXI(z5,p1,ap1,2);

					FNEG_ZSP(z0,p0,z0);
					FNEG_ZSP(z1,p0,z1);
					FNEG_ZSP(z2,p0,z2);
					FNEG_ZSP(z3,p0,z3);
					FNEG_ZSP(z4,p1,z4);
					FNEG_ZSP(z5,p1,z5);

					ST1W_ZXI(z0,p0,awp ,0);
					ST1W_ZXI(z1,p0,awp2,0);
					ST1W_ZXI(z2,p0,awp ,1);
					ST1W_ZXI(z3,p0,awp2,1);
					ST1W_ZXI(z4,p1,awp ,2);
					ST1W_ZXI(z5,p1,awp2,2);

					ap0+=mm;
					ap1+=mm;
				}
				else if(mm>nvl) {
					j=nvl*1;
					WHILELT_PSX(p1,j,mm);

                                        apre0=ap0+lda*4;
                                        apre1=ap1+lda*4;
                                        PRFM_XI(PLDL2STRM,apre0,0);
                                        PRFM_XI(PLDL2STRM,apre1,0);

					LD1W_ZXI(z0,p0,ap0,0);
					LD1W_ZXI(z1,p0,ap1,0);
					LD1W_ZXI(z2,p1,ap0,1);
					LD1W_ZXI(z3,p1,ap1,1);

					FNEG_ZSP(z0,p0,z0);
					FNEG_ZSP(z1,p0,z1);
					FNEG_ZSP(z2,p1,z2);
					FNEG_ZSP(z3,p1,z3);

					ST1W_ZXI(z0,p0,awp ,0);
					ST1W_ZXI(z1,p0,awp2,0);
					ST1W_ZXI(z2,p1,awp ,1);
					ST1W_ZXI(z3,p1,awp2,1);

					ap0+=mm;
					ap1+=mm;
				}
				else if(mm>=1) {
					j=0;
					WHILELT_PSX(p1,j,mm);

                                        apre0=ap0+lda*4;
                                        apre1=ap1+lda*4;
                                        PRFM_XI(PLDL2STRM,apre0,0);
                                        PRFM_XI(PLDL2STRM,apre1,0);

					LD1W_ZXI(z0,p1,ap0,0);
					LD1W_ZXI(z1,p1,ap1,0);

					FNEG_ZSP(z0,p1,z0);
					FNEG_ZSP(z1,p1,z1);

					ST1W_ZXI(z0,p1,awp ,0);
					ST1W_ZXI(z1,p1,awp2,0);
					ap0+=mm;
					ap1+=mm;
				}

			}
			ap0+=lda*2-m;
			ap1+=lda*2-m;
		}
		for(;ik<=ike;ik++) {
			im=0;
			awp=aw+nul*ik;

			if(m-im>=nul) {
				im+=nul;
				apre0=ap0+ipre;
				PRFM_XI(PLDL1STRM,apre0,0);
				PRFM_XI(PLDL1STRM,apre0,256);

				LD1W_ZXI(z0,p0,ap0,0);
				LD1W_ZXI(z1,p0,ap0,1);
				LD1W_ZXI(z2,p0,ap0,2);
				LD1W_ZXI(z3,p0,ap0,3);

				FNEG_ZSP(z0,p0,z0);
				FNEG_ZSP(z1,p0,z1);
				FNEG_ZSP(z2,p0,z2);
				FNEG_ZSP(z3,p0,z3);

				for(;im<m-nul+1;im+=nul) {
					ap0+=nul;
					apre0=ap0+ipre;
					PRFM_XI(PLDL1STRM,apre0,0);
					PRFM_XI(PLDL1STRM,apre0,256);

					ST1W_ZXI(z0,p0,awp,0);
					ST1W_ZXI(z1,p0,awp,1);
					ST1W_ZXI(z2,p0,awp,2);
					ST1W_ZXI(z3,p0,awp,3);

					LD1W_ZXI(z0,p0,ap0,0);
					LD1W_ZXI(z1,p0,ap0,1);
					LD1W_ZXI(z2,p0,ap0,2);
					LD1W_ZXI(z3,p0,ap0,3);

					FNEG_ZSP(z0,p0,z0);
					FNEG_ZSP(z1,p0,z1);
					FNEG_ZSP(z2,p0,z2);
					FNEG_ZSP(z3,p0,z3);

					awp+=k*nul;
				}

				ST1W_ZXI(z0,p0,awp,0);
				ST1W_ZXI(z1,p0,awp,1);
				ST1W_ZXI(z2,p0,awp,2);
				ST1W_ZXI(z3,p0,awp,3);

				awp+=k*nul;
				ap0+=nul;
			}

			mm=m-im;
			awp=aw+k*(m-mm)+ik*mm;

			if(mm>nvl*3) {
				j=nvl*3;
				WHILELT_PSX(p1,j,mm);
				LD1W_ZXI(z0,p0,ap0,0);
				LD1W_ZXI(z2,p0,ap0,1);
				LD1W_ZXI(z4,p0,ap0,2);
				LD1W_ZXI(z6,p1,ap0,3);

				FNEG_ZSP(z0,p0,z0);
				FNEG_ZSP(z2,p0,z2);
				FNEG_ZSP(z4,p0,z4);
				FNEG_ZSP(z6,p1,z6);

				ST1W_ZXI(z0,p0,awp ,0);
				ST1W_ZXI(z2,p0,awp ,1);
				ST1W_ZXI(z4,p0,awp ,2);
				ST1W_ZXI(z6,p1,awp ,3);

				ap0+=mm;
			}
			else if(mm>nvl*2) {
				j=nvl*2;
				WHILELT_PSX(p1,j,mm);

				LD1W_ZXI(z0,p0,ap0,0);
				LD1W_ZXI(z2,p0,ap0,1);
				LD1W_ZXI(z4,p1,ap0,2);

				FNEG_ZSP(z0,p0,z0);
				FNEG_ZSP(z2,p0,z2);
				FNEG_ZSP(z4,p1,z4);

				ST1W_ZXI(z0,p0,awp ,0);
				ST1W_ZXI(z2,p0,awp ,1);
				ST1W_ZXI(z4,p1,awp ,2);

				ap0+=mm;
			}
			else if(mm>nvl) {
				j=nvl*1;
				WHILELT_PSX(p1,j,mm);

				LD1W_ZXI(z0,p0,ap0,0);
				LD1W_ZXI(z1,p1,ap0,1);

				FNEG_ZSP(z0,p0,z0);
				FNEG_ZSP(z1,p1,z1);

				ST1W_ZXI(z0,p0,awp ,0);
				ST1W_ZXI(z1,p1,awp ,1);
				ap0+=mm;
			}
			else if(mm>=1) {
				j=0;
				WHILELT_PSX(p1,j,mm);

				LD1W_ZXI(z0,p1,ap0,0);

				FNEG_ZSP(z0,p1,z0);

				ST1W_ZXI(z0,p1,awp,0);

				ap0+=mm;
			}

			ap0+=lda-m;
		}
	}
	return 0;
}
