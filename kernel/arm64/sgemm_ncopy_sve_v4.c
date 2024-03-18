/***************************************************************************
(c) RIKEN 2024, 2024. All rights reserved. sgemm_ncopy_sve_v4.c 0.3.26
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
#include <arm_sve.h>

int CNAME(BLASLONG k, BLASLONG m, IFLOAT *a, BLASLONG lda, IFLOAT *aw)
{
#define A(ik,im) a[(ik)+(im)*lda]
	BLASLONG im,ik,iaw,ii;
	BLASLONG imb,kk,mm,mr,mbr,kbr;
	BLASLONG mbw;
	BLASLONG nvl=svcntw();
	BLASLONG nsi=4,mb=nvl*nsi;

	for(imb=0;imb<m-mb+1;imb+=mb) {
		for(ik=0;ik<k-mb+1;ik+=mb) {
			iaw=k*imb  +ik*mb;
			for(im=0;im<mb;im++) {
				if(imb+im+mb<m) {
					__builtin_prefetch(&(A(ik+ii,imb+im+mb)),0,2);
					__builtin_prefetch(&(A(ik+ii+31,imb+im+mb)),0,2);
				}
				__builtin_prefetch(&(A(ik+ii,imb+im)),0,0);
				__builtin_prefetch(&(A(ik+ii+31,imb+im)),0,0);
				for(ii=0;ii<mb;ii++) {
					aw[iaw+mb*ii]=A(ik+ii,imb+im);
				}
				iaw=iaw+1;
			}
		}
		kk=ik;
		for(im=0;im<mb;im++) {
			iaw=im+k*imb+kk*mb;
			if(imb+im+mb<m) {
				__builtin_prefetch(&(A(kk,imb+im+mb)),0,2);
			}
			__builtin_prefetch(&(A(kk,imb+im)),0,0);
			for(ik=kk;ik<k;ik++) {
				aw[iaw]=A(ik,imb+im);
				iaw=iaw+mb;
			}
		}
	}
	mbr=m%mb;
	mm=m-mbr;
	mr=mbr;
	if(mr!=0) {
		for(ik=0;ik<k-mb+1;ik+=mb) {
			iaw=k*mm +ik*mr;
			for(im=mm;im<mm+mr;im++) {
				if(imb+im+mb<m) {
					__builtin_prefetch(&(A(ik+ii,imb+mb)),0,2);
					__builtin_prefetch(&(A(ik+ii+31,imb+mb)),0,2);
				}
				__builtin_prefetch(&(A(ik+ii,im-1)),0,0);
				__builtin_prefetch(&(A(ik+ii+31,im-1)),0,0);
				for(ii=0;ii<mb;ii++) {
					aw[iaw+mr*ii]=A(ik+ii,im);
				}
				iaw=iaw+1;
			}
		}
		kbr=k%mb;
		kk=k-kbr;
		iaw=k*mm+kk*mr;
		for(ik=kk;ik<k;ik++) {
			for(im=mm;im<mm+mr;im++) {
				aw[iaw]=A(ik,im);
				iaw=iaw+1;
			}
		}
	}
	return 0;
}
