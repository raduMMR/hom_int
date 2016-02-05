#ifndef CTXTBIT_
#define CTXTBIT_

#include "Hom_NTRU_library/def.h"
#include "Hom_NTRU_library/ntru.h"
#include "Hom_NTRU_library/crt.h"
#include "Hom_NTRU_library/fft_mult.h"
#include "Hom_NTRU_library/general.h"

/*
clasa CtxtBit este clona nereusita a clasei CtxtPart din
biblioteca HElib pentru contextul NTRU.
In aceasta clasa retin polinomul corespunzator ctxt-ului
indexul corespunzator nivelului din "scara de moduli" pe 
care ma aflu, contextul NTRU(echivalent FHEcontext din HElib)
*/
class CtxtBit
{
	ZZX ctxt;
	int qIndex;
	// int num_relin;
	// int num_mul;
	ntru *n;
	
public:
	
	CtxtBit(ntru *n);
	
	CtxtBit(ntru *n, int bit);
	
	CtxtBit& operator+(CtxtBit& ctxt2);
	
	CtxtBit& operator+=(CtxtBit& ctxt2);
	
	CtxtBit& operator*(CtxtBit& ctxt2);
	
	CtxtBit& operator*=(CtxtBit& ctxt2);
	
	int decryptBit()const;
	
	void setCtxt(ZZX ctxt)
	{
		this->ctxt = ctxt;
	}
	
	ZZX getCtxt()const
	{
		return ctxt;
	}
	
	// void modSwitch();
	// void reliniarizeCtxt();
	// void encryptBit(int bit);
};

#endif