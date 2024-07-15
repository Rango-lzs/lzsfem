/**
 * @file matconst.h
 * Defines several material constant (respective their
 * representative number).
 *
 * These constant are used to access material parameters  at material level from
 * property dictionary.
 */

#ifndef matconst_h
#define matconst_h


namespace fem
{
	/* Poissons coeficients */
#define NYxz 300
#define NYyz 301
#define NYxy 302
#define NYzx 303
#define NYzy 304
#define NYyx 305

/* Ortho Young moduli */
#define Ex 400
#define Ey 401
#define Ez 402

/* Thermal dilatation coeffs */
#define tAlphax 403
#define tAlphay 404
#define tAlphaz 405
#define tAlpha 406 // 1d.

/* shear modulus */
#define Gyz 407
#define Gxz 408
#define Gxy 409

/* heat capacity */
#define HeatCapaCoeff 450
#define Mass1CapaCoeff 451

/* viscosity */
#define Viscosity 500
#define YieldStress 501

/* parameters of damage model */
#define e0_ID  800
#define ef_ID  801
#define gf_ID  802
#define ek_ID  803
#define wf_ID  804
#define gft_ID 805

/* strengths */
#define ft_strength 806
#define fc_strength 807

/* nonlocal material parameters */
#define AVERAGING_TYPE 901
#define exponent_ID 902
#define rf_ID 903

/* lattice parameters */
#define eNormal_ID 1101

} // end namespace fem
#endif // matconst_h
