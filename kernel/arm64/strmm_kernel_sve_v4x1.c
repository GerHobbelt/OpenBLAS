/***************************************************************************
(c) RIKEN 2024, 2024. All rights reserved. strmm_kernel_sve_v4x1.c 0.3.26
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

int CNAME(BLASLONG bm,BLASLONG bn,BLASLONG bk,FLOAT balpha,FLOAT* ba,FLOAT* bb,FLOAT* C,BLASLONG ldc ,BLASLONG offset)
{

   BLASLONG im,j,k;
   FLOAT *C0,*ptrba,*ptrbb;

   svfloat32_t alpha;

   svfloat32_t res0_0;
   svfloat32_t res0_1;
   svfloat32_t res0_2;
   svfloat32_t res0_3;

   svfloat32_t a0;
   svfloat32_t a1;
   svfloat32_t a2;
   svfloat32_t a3;

   svfloat32_t b0;

   BLASLONG off, temp;

   BLASLONG size;
   BLASLONG sve_size = svcntw();
   BLASLONG unroll_size = sve_size*4;

   svbool_t p0,p1;
   p0=svptrue_b32();
   alpha=svdup_f32(balpha);

#if !defined(LEFT)
   off = -offset;
#else
   off = 0;
#endif

   for (j=0; j<bn; j+=1)
   {
        C0 = C;
#if defined(TRMMKERNEL) && defined(LEFT)
		off = offset;
#endif
        ptrba = ba;

        for (im=0; im<bm-unroll_size+1; im+=unroll_size)
        {

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
		ptrbb = bb;
#else
		ptrba += off * unroll_size;
		ptrbb = bb + off;
#endif
                res0_0 = svdup_f32(0);
                res0_1 = svdup_f32(0);
                res0_2 = svdup_f32(0);
                res0_3 = svdup_f32(0);

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
		temp = bk-off;
#elif defined(LEFT)
		temp = off+unroll_size;	// number of values in A
#else
		temp = off+1;	// number of values in B
#endif

		for (k=0; k<temp; k++)
                {
			b0 = svdup_f32_x(p0,*ptrbb);

			a0 = svld1(p0,ptrba);
			a1 = svld1_vnum(p0,ptrba,1);
			a2 = svld1_vnum(p0,ptrba,2);
			a3 = svld1_vnum(p0,ptrba,3);
			res0_0=svmla_x(p0,res0_0,a0,b0);
			res0_1=svmla_x(p0,res0_1,a1,b0);
			res0_2=svmla_x(p0,res0_2,a2,b0);
			res0_3=svmla_x(p0,res0_3,a3,b0);

			ptrba = ptrba+unroll_size;
			ptrbb = ptrbb+1;
                }

		res0_0=svmul_x(p0,res0_0,alpha);
		res0_1=svmul_x(p0,res0_1,alpha);
		res0_2=svmul_x(p0,res0_2,alpha);
		res0_3=svmul_x(p0,res0_3,alpha);

		svst1(p0,C0,res0_0);
		svst1_vnum(p0,C0,1,res0_1);
		svst1_vnum(p0,C0,2,res0_2);
		svst1_vnum(p0,C0,3,res0_3);


#if ( defined(LEFT) && defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
		temp = bk - off;
#ifdef LEFT
		temp -= unroll_size; // number of values in A
#else
		temp -= 1; // number of values in B
#endif
		ptrba += temp*unroll_size;
		ptrbb += temp;
#endif

#ifdef LEFT
		off += unroll_size; // number of values in A
#endif

		C0 = C0+unroll_size;
	}

	if ( bm-im >sve_size*3 )
	{
                size=bm-im;
                p1=svwhilelt_b32(im+sve_size*3,bm);

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
		ptrbb = bb;
#else
		ptrba += off*size;
		ptrbb = bb + off;
#endif

		res0_0 = svdup_f32(0);
		res0_1 = svdup_f32(0);
		res0_2 = svdup_f32(0);
		res0_3 = svdup_f32(0);


#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
		temp = bk-off;
#elif defined(LEFT)
		temp = off+size;	// number of values in A
#else
		temp = off + 1;	// number of values in B
#endif

		for (k=0; k<temp; k++)
                {
                        b0 = svdup_f32_x(p0,*ptrbb);

                        a0 = svld1(p0,ptrba);
                        a1 = svld1_vnum(p0,ptrba,1);
                        a2 = svld1_vnum(p0,ptrba,2);
                        a3 = svld1_vnum(p1,ptrba,3);
                        res0_0=svmla_x(p0,res0_0,a0,b0);
                        res0_1=svmla_x(p0,res0_1,a1,b0);
                        res0_2=svmla_x(p0,res0_2,a2,b0);
                        res0_3=svmla_x(p1,res0_3,a3,b0);

			ptrba = ptrba+size;
			ptrbb = ptrbb+1;
                }

                res0_0=svmul_x(p0,res0_0,alpha);
                res0_1=svmul_x(p0,res0_1,alpha);
                res0_2=svmul_x(p0,res0_2,alpha);
                res0_3=svmul_x(p1,res0_3,alpha);

                svst1(p0,C0,res0_0);
                svst1_vnum(p0,C0,1,res0_1);
                svst1_vnum(p0,C0,2,res0_2);
                svst1_vnum(p1,C0,3,res0_3);

#if ( defined(LEFT) && defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
		temp = bk - off;
#ifdef LEFT
		temp -= size; // number of values in A
#else
		temp -= 1; // number of values in B
#endif
		ptrba += temp*size;
		ptrbb += temp;
#endif

#ifdef LEFT
		off += size; // number of values in A
#endif

		C0 = C0+size;

	}

	else if ( bm - im > sve_size*2 )
	{
                size=bm-im;
                p1=svwhilelt_b32(im+sve_size*2,bm);

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
		ptrbb = bb;
#else
		ptrba += off*size;
		ptrbb = bb + off;
#endif

                res0_0 = svdup_f32(0);
                res0_1 = svdup_f32(0);
                res0_2 = svdup_f32(0);

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
		temp = bk-off;
#elif defined(LEFT)
		temp = off+size; // number of values in A
#else
		temp = off+1;	// number of values in B
#endif

		for (k=0; k<temp; k++)
                {
                        b0 = svdup_f32_x(p0,*ptrbb);

                        a0 = svld1(p0,ptrba);
                        a1 = svld1_vnum(p0,ptrba,1);
                        a2 = svld1_vnum(p1,ptrba,2);
                        res0_0=svmla_x(p0,res0_0,a0,b0);
                        res0_1=svmla_x(p0,res0_1,a1,b0);
                        res0_2=svmla_x(p1,res0_2,a2,b0);

			ptrba = ptrba+size;
                        ptrbb = ptrbb+1;
                }

                res0_0=svmul_x(p0,res0_0,alpha);
                res0_1=svmul_x(p0,res0_1,alpha);
                res0_2=svmul_x(p1,res0_2,alpha);

                svst1(p0,C0,res0_0);
                svst1_vnum(p0,C0,1,res0_1);
                svst1_vnum(p1,C0,2,res0_2);

#if ( defined(LEFT) && defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
		temp = bk - off;
#ifdef LEFT
		temp -= size; // number of values in A
#else
		temp -= 1; // number of values in B
#endif
		ptrba += temp*size;
		ptrbb += temp;
#endif

#ifdef LEFT
		off += size; // number of values in A
#endif

		C0 = C0+size;
 
	}

	else if ( bm - im>sve_size )
	{
                size=bm-im;
                p1=svwhilelt_b32(im+sve_size,bm);

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
		ptrbb = bb;
#else
		ptrba += off*size;
		ptrbb = bb + off;
#endif

                res0_0 = svdup_f32(0);
                res0_1 = svdup_f32(0);

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
		temp = bk-off;
#elif defined(LEFT)
		temp = off+size; // number of values in A
#else
		temp = off+1;	// number of values in B
#endif

		for (k=0; k<temp; k++)
                {
                        b0 = svdup_f32_x(p0,*ptrbb);
                        a0 = svld1(p0,ptrba);
                        a1 = svld1_vnum(p1,ptrba,1);
                        res0_0=svmla_x(p0,res0_0,a0,b0);
                        res0_1=svmla_x(p1,res0_1,a1,b0);

			ptrba = ptrba+size;
                        ptrbb = ptrbb+1;
                }
                res0_0=svmul_x(p0,res0_0,alpha);
                res0_1=svmul_x(p1,res0_1,alpha);

                svst1(p0,C0,res0_0);
                svst1_vnum(p1,C0,1,res0_1);

#if ( defined(LEFT) && defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
		temp = bk - off;
#ifdef LEFT
		temp -= size; // number of values in A
#else
		temp -= 1; // number of values in B
#endif
		ptrba += temp*size;
		ptrbb += temp;
#endif

#ifdef LEFT
		off += size; // number of values in A
#endif

		C0 = C0+size;
	}

	else 
	{
                size=bm-im;
                p1=svwhilelt_b32(im,bm);

#if (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
		ptrbb = bb;
#else
		ptrba += off*size;
		ptrbb = bb + off;
#endif

                res0_0 = svdup_f32(0);

#if (defined(LEFT) && !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
		temp = bk-off;
#elif defined(LEFT)
		temp = off+size;	// number of values in A
#else
		temp = off+1;	// number of values in B
#endif

		for (k=0; k<temp; k++)
                {
                        b0 = svdup_f32_x(p0,*ptrbb);
                        a0 = svld1(p1,ptrba);
                        res0_0=svmla_x(p1,res0_0,a0,b0);
			ptrba = ptrba+size;
                        ptrbb = ptrbb+1;
                }
                res0_0=svmul_x(p1,res0_0,alpha);
                svst1(p1,C0,res0_0);

#if ( defined(LEFT) && defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
		temp = bk - off;
#ifdef LEFT
		temp -= size; // number of values in A
#else
		temp -= 1; // number of values in B
#endif
		ptrba += temp*size;
		ptrbb += temp;
#endif

#ifdef LEFT
		off += size; // number of values in A
#endif

		C0 = C0+size;
	}


#if defined(TRMMKERNEL) && !defined(LEFT)
		off += 1;
#endif

        bb = bb+bk;
        C += ldc;
    }

   return 0;
}
