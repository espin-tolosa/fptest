#ifndef FPTEST_INTERNALS_H
#define FPTEST_INTERNALS_H

#include "fptest.h"

/* Floating-Point Layout (valid for arc=x86_64) */

/*
 * f80_t (Extended precision) for arch=x86_64 with x87 FPU
 *
 * In terms of memory alignment, the extender precision type is padded to span 16 bytes
 * But only the lower address corresponding to the first 10 bytes are meaningful.
 */

#define LW0   	4 /* contains MSB */
#define LW1   	3
#define LW2   	2
#define LW3   	1
#define LW4   	0 /* contains LSB */

/* u64_t access to f80_t */

#define LU0     1 /* contains MSB */
#define LU1     0 /* contains LSB */

/*  f64_t (Double precision) */

#define DW0   	3 /* contains MSB */
#define DW1   	2
#define DW2   	1
#define DW3   	0 /* contains LSB */

/*  f32_t (Single precision) */

#define FW0   	1 /* contains MSB */
#define FW1   	0 /* contains LSB */

typedef union
{
	u16_t   w[ 4 ];
	u64_t   u;
	f64_t   f;
} dw_t;

typedef union
{
	u16_t   w[ 2 ];
	u32_t   u;
	f32_t   f;
} fw_t;

/* nbits sign, exponent, digits: 1, 8, 23 (+1 digit bit hidden) */
#define FW0_SMASK 0x8000
#define FW0_EMASK 0x7F80
#define FW0_DMASK 0x007F
#define FWx_DMASK 0xFFFF
#define DOFFS_F32 7         /* offset bits of exponent within FW0 */
#define EMAX_F32  127
#define EMIN_F32  ( 1 - EMAX_F32 )
#define NDIG_F32  24

/* nbits sign, exponent, digits: 1, 11, 52 (+1 digit bit hidden) */
#define DW0_SMASK 0x8000
#define DW0_EMASK 0x7FF0
#define DW0_DMASK 0x000F
#define DWx_DMASK 0xFFFF
#define DOFFS_F64 4         /* offset bits of exponent within DW0 */
#define EMAX_F64  1023
#define EMIN_F64  ( 1 - EMAX_F64 )
#define NDIG_F64  53

/* Constants */
static const fw_t m_inff32 	= { .u = 0x7F800000u, };
static const fw_t m_maxf32 	= { .u = 0x7F7FFFFFu, };
static const fw_t m_zerof32	= { .u = 0x00000000u, };
static const fw_t m_minf32 	= { .u = 0x00800001u, };
static const fw_t m_subnf32	= { .u = 0x00000001u, };
static const fw_t m_qnanf32	= { .u = 0x7FC00000u, };
static const fw_t m_epsf32	= { .u = 0x34000000u, };

static const dw_t m_inff64 	= { .u = 0x7FF0000000000000uL, };
static const dw_t m_maxf64 	= { .u = 0x7FEFFFFFFFFFFFFFuL, };
static const dw_t m_zerof64	= { .u = 0x0000000000000000uL, };
static const dw_t m_minf64 	= { .u = 0x0010000000000001uL, };
static const dw_t m_subnf64	= { .u = 0x0000000000000001uL, };
static const dw_t m_qnanf64	= { .u = 0x7FF8000000000000uL, };
static const fw_t m_epsf64	= { .f = (f64_t) 2.22044604925031308084726333618164062e-16L, };

/* Utilities */
#define F32_UNBIASED_TO_BIASED_EXP(e) ( (e) < EMIN_F32 ? EMIN_F32 : ( (e) > EMAX_F32 ? EMAX_F32 : (e) + EMAX_F32 ) )
#define F64_UNBIASED_TO_BIASED_EXP(e) ( (e) < EMIN_F64 ? EMIN_F64 : ( (e) > EMAX_F64 ? EMAX_F64 : (e) + EMAX_F64 ) )

#define F32_SET_EXP(w0,new_unbiased_exp) do { \
	( w0 ) = ( ( w0 ) & ~( (u16_t) FW0_EMASK ) ) | ( ( F32_UNBIASED_TO_BIASED_EXP(new_unbiased_exp) ) << DOFFS_F32 );\
} while ( 0 )

#define F64_SET_EXP(w0,new_unbiased_exp)\
do { \
	( w0 ) = ( ( w0 ) & ~( (u16_t) DW0_EMASK ) ) | ( ( F64_UNBIASED_TO_BIASED_EXP(new_unbiased_exp) ) << DOFFS_F64 );\
} while ( 0 )

#define F32_GET_EXP(x) ( (i16_t) ( ( ((fw_t){ .f = x }).w[ FW0 ] & FW0_EMASK ) >> DOFFS_F32 ) - EMAX_F32 )
#define F64_GET_EXP(x) ( (i16_t) ( ( ((dw_t){ .f = x }).w[ DW0 ] & DW0_EMASK ) >> DOFFS_F64 ) - EMAX_F64 )

#endif/*FPTEST_INTERNALS_H*/
