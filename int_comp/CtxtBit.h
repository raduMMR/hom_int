#ifndef CTXTBIT_
#define CTXTBIT_

#include "def.h"
#include "ntru.h"
#include "crt.h"
#include "fft_mult.h"
#include "general.h"

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