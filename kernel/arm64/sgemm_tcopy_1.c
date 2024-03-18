/***************************************************************************
(c) RIKEN 2024, 2024. All rights reserved. sgemm_tcopy_1.c 0.3.26
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

int CNAME(BLASLONG n, BLASLONG m, FLOAT *a, BLASLONG lda, FLOAT *b)
{
	BLASLONG ldb=n;
	BLASLONG in,im,iin,iim,nii,mii;
	const static BLASLONG NB=32;
#define A(im,in) a[(im)+(in)*lda]
#define B(in,im) b[(in)+(im)*ldb]

	if(m==5) {
		for(in=0;in<n;in++) {
			__builtin_prefetch(&(A(9,in)),0,2);
			B(in,0)=A(0,in);
			B(in,1)=A(1,in);
			B(in,2)=A(2,in);
			B(in,3)=A(3,in);
			B(in,4)=A(4,in);
		}
		return 0;
	}
	for(in=0;in<n-NB+1;in+=NB) {
		for(im=0;im<m-NB+1;im+=NB) {
			for(iin=0;iin<NB;iin++) {
				__builtin_prefetch(&(A(im+iim   ,in+iin+16)),0,2);
				__builtin_prefetch(&(A(im+iim+31,in+iin+16)),0,2);
				__builtin_prefetch(&(A(im+iim   ,in+iin+2)),0,0);
				__builtin_prefetch(&(A(im+iim+31,in+iin+2)),0,0);
				for(iim=0;iim<NB;iim++) {
					B(in+iin,im+iim)=A(im+iim,in+iin);
				}
			}
		}
		mii=im;
		for(im=mii;im<=m;im++) {
			for(iin=0;iin<NB;iin++) {
				B(in+iin,im)=A(im,in+iin);
			}

		}
	}

	nii=in;
	for(im=0;im<m-NB+1;im+=NB) {
		for(in=nii;in<n;in++) {
			if(n-nii>16) {
				__builtin_prefetch(&(A(im+iim   ,in+iin+16)),0,2);
				__builtin_prefetch(&(A(im+iim+31,in+iin+16)),0,2);
			}
			__builtin_prefetch(&(A(im+iim   ,in+iin+2)),0,0);
			__builtin_prefetch(&(A(im+iim+31,in+iin+2)),0,0);
			for(iim=0;iim<NB;iim++) {
				B(in,im+iim)=A(im+iim,in);
			}
		}
	}
	mii=im;
	for(iim=mii;iim<m;iim++) {
		for(in=nii;in<n;in++) {
			B(in,iim)=A(iim,in);
		}
	}
	return 0;
}
