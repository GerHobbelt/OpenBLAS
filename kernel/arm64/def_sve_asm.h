/***************************************************************************
(c) RIKEN 2025, 2025. All rights reserved. def_sve_asm.h 0.3.26
Copyright 2025 FUJITSU limited
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

#ifndef __DEF_SVE_ASM_INLINE
#define __DEF_SVE_ASM_INLINE

// Argument i must be a multiple of 8 in the range 0 to 32760
#define PRFM_XI(prfop,x,i) PRFM_XI_base(prfop,x,i)
#define PRFM_XI_base(prfop,x,i) \
  __asm__ __volatile__("\t\t\tprfm "#prfop",[%0,#"#i"]"::"r"(x));

// Argument i is in the range -256 to 255
#define PRFUM_XI(prfop,x,i) PRFUM_XI_base(prfop,x,i)
#define PRFUM_XI_base(prfop,x,i) \
  __asm__ __volatile__("\t\t\tprfum "#prfop",[%0,#"#i"]"::"r"(x));

#define PRFM_XXLSL3(prfop,x,i) PRFM_XXLSL3_base(prfop,x,i)
#define PRFM_XXLSL3_base(prfop,x,i) \
  __asm__ __volatile__("\t\t\tprfm "#prfop",[%0,%1,LSL #3]"::"r"(x),"r"(i));

#define CNTW_X(x) CNTW_X_base(x)
#define CNTW_X_base(x) \
  __asm__ __volatile__("\t\t\tcntw %0":"=r"(x));

#define NOP \
  __asm__ __volatile__("\t\t\tnop");

#define DUP_ZSI(Z,X) DUP_ZSI_base(Z,X)
#define DUP_ZSI_base(Z,X) \
	        __asm__ __volatile__("\t\t\tdup\t"#Z".s,#"#X:::#Z);

#define DUP_ZDI(Z,X) DUP_ZDI_base(Z,X)
#define DUP_ZDI_base(Z,X) \
	        __asm__ __volatile__("\t\t\tdup\t"#Z".d,#"#X:::#Z);

#define PLDL1KEEP 0
#define PLDL1STRM 1
#define PLDL2KEEP 2
#define PLDL2STRM 3
#define PLDL3KEEP 4
#define PLDL3STRM 5
#define PSTL1KEEP 8
#define PSTL1STRM 9
#define PSTL2KEEP 10
#define PSTL2STRM 11
#define PSTL3KEEP 12
#define PSTL3STRM 13

#include "def_sve_asm_float.h"
#include "def_sve_asm_double.h"
#include "def_sve_asm_predicate.h"

#endif   /* __DEF_SVE_ASM_INLINE */

