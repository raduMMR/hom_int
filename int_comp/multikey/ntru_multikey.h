#ifndef NTRU_H_MULTIKEY
#define NTRU_H_MULTIKEY

#include <iostream>
#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZXFactoring.h>

#include "../hNTRU/general.h"
#include "../hNTRU/fft_mult.h"
#include "multikey_structs.h"
#include <stdlib.h>
#include <ctime>

using namespace std;
NTL_CLIENT

/*
@brief clasa pentru evaluarea homomorfica ntru multikey
		genereaza cheile si efectueaza operatiile homomorfice
		cheile si ciphertext-urile sunt reprezentate ca polinoame ZZX
		
		pentru evaluarea multikey 
		-------------------------
		1. clasa NtruMultikey retine un vector de pointeri la oate instantele care vor fi create 
		2. metoda compute calculeaza functia homomorfica multikey pe baza indecsilor 
		   corespunzatori partilor si a descrierii polinomului de evaluat
		   
*/
class NtruMultikey{
public:
///////////////////////////////////////////////////////////////////////////////
	NtruMultikey(ZZ my_p, ZZ my_B, int my_N, int my_dm);

	////////////////////////////////
	void	ModulusFind(int num, int max_bit, int diff);
	void	ComputeKeys(int num);
	////////////////////////////////
	void   	ModulusFindRing(int num, int max_bit, int diff, ZZX modu);
	void 	ComputeKeysRingNoRelin(int num);

	////////////////////////////////
	ZZX 	Encrypt(ZZX m, int index);
	ZZX		Decrypt(ZZX c, int index);
	ZZX 	Sample();
	ZZX		Relin(ZZX c, int index);
	ZZX 	ModSwitch(ZZX &x, int index);
	////////////////////////////////
	void 	SetModulus();
	ZZX 	ComputeFastCycModulus(int n);
	void 	SetFastPolyModulus();
	void 	compute_eval(eval_key &ek2, ZZX tppk, ZZX tpek, int index);
	////////////////////////////////
	void	PolyReduce(ZZX &out, ZZX &in);
	void	CoeffReduce(ZZX &out, ZZX &in, int index);
	void 	PolyCoeffReduce(ZZX &out, ZZX &in, int index);
	ZZX		MulModZZX(ZZX &a, ZZX &b, int index);
	ZZX		AddModZZX(ZZX &a, ZZX &b, int index);
	ZZX 	&ReturnPolyMod();
	int		&ReturnPolyModDegree();
	ZZ		&ReturnModQ(int index);
///////////////////////////////////////////////////////////////////////////////
	void 	PrintKeys(int index);
	void 	PrintMods(int index);
	void 	IOReadModulus();
	void 	IOReadModulusExt();
	
	void 	IOReadModulus_Special();
	void 	IOReadModulusExt_Special();
///////////////////////////////////////////////////////////////////////////////
	void	ComputeKeysRingRelin(int num);
	void	compute_eval_onerelin(eval_key &ek2, ZZX tppk, ZZX tpek, int index);
	ZZX		RelinRing(ZZX c, int index);
///////////////////////////////////////////////////////////////////////////////
	void		ComputeKeysRingRelin_FFT(int num, int tblsize);
	void		ConvertKeystoFFT(eval_key &ek2, FFTPolyList &fftkey, int &bitsize, int &num_add);
	ZZX 		RelinRingFFT(ZZX &c, int index);
	ZZX 		RelinRingFFT(ZZX &c, int index, int tableIndx);
	void 	ReduceKeysLevel(int level);
	void 	ReduceKeysLevel(int level, int tblIndex);
	void	GetPolyBitsIndex(ZZX &bits, ZZX &cip, int &bitIndex);
	void 	FFTmult(fftRep *res, fftRep *R1, fftRep *R2, int &priNum);
	void 	FFTadd(fftRep *res, fftRep *R1, fftRep *R2, int &priNum);
	void 	FromCrt2Poly(ZZX &res, zz_pX *polys, int &priNum);
	///////////////////////////////////////////////////////////////////////////////
	//								FOR X^M-1 Poly								 //
	///////////////////////////////////////////////////////////////////////////////
	void	compute_eval_XN1(eval_key &ek2, ZZX tppk, ZZX tpek, int index);
	ZZX		Encrypt_XN1(ZZX m, int i);
	ZZX		MulModZZX_XN1(ZZX &a, ZZX &b, int index);
	ZZX 	Relin_XN1(ZZX c, int index);
	ZZX		ModSwitch_XN1(ZZX &x, int index);
	ZZX		Decrypt_XN1(ZZX c, int index);
	ZZX 	ModSwitchX_M_1(ZZX &x, int index);
	ZZX		PolyModX_N_1(ZZX &in);
	ZZX		PolyModX_Np1(ZZX &in);
	ZZX		MulModPolyX_M_1(ZZX &a, ZZX &b, int index);
	ZZX 	&ReturnPolyMod_M_1();
///////////////////////////////////////////////////////////////////////////////
	ZZX 	CreateMessage(int size);
	ZZX		FFTTestFunc(fftRep *R, int priNum);
	
	/************************* multikey stuff **************************************/
	
	/*
	@brief calculeaza o functie descrisa de un polinom pe mai multe ctxt-uri 
			apartinand mai multor parti 
	@param poly_description reprezinta o descriere a functiei delegate si a partilor
			implicate in computatie
	*/
	static HomEvalResult compute_function(PolyDescription poly_description);
	
	

private:
	int N, degree_m;
	ZZ 	p, q, B;
	ZZX modulus, modulus_m_1;

	ZZ  	*q_list, *pr_list;
	ZZX 	*pk, *sk;
	vec_ZZ 	*q_list_mult;

	struct eval_key *ek, *ek2;
	FFTPolyList *fftKeys;
	int tableSize;

	myReduction myr;
	bool reducsetflag;
	int rand_seed, message_rand, max_bitsize, init_bitsize;

	ZZ mySeed;
	
	/************************* multikey stuff **************************************/
	
};

///////////////////////////////////////////////////////////////////////////////
void 	clear(ZZX &x, int size);
void 	find_inverse(ZZX &f_inv, ZZX &f, ZZ &q, int &degree, bool &isfound, ZZX modulus);
void 	coeff_reduction(ZZX &out, ZZX &in, ZZ &q, int &degree);
void 	coeff_reduction_q_2(ZZX &out, ZZX &in, ZZ &q, int N);
void 	getClosestMod2(ZZ& out, ZZ& in, ZZ& p1, ZZ& p2);
ZZ 		ComputePR(ZZ q, ZZ m);
vec_ZZ 	findFactors(ZZ f);
int 	MobuisFunction(int n);
///////////////////////////////////////////////////////////////////////////////

#endif /* NTRU_H_MULTIKEY */




