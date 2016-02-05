
#include "CtxtBit.h"

CtxtBit::CtxtBit(ntru *n)
{
	this->n = n;
	qIndex = 0;
	ctxt = 0;
}

CtxtBit::CtxtBit(ntru *n, int bit)
{
	this->n = n;
	qIndex = 0;
	
	// num_relin = 0;
	// num_mul = 0;
	
	ZZX message;
	message = bit;
	ctxt = n->Encrypt(message, 0);
}

/*************************************************************************************/
CtxtBit& CtxtBit::operator+(CtxtBit& ctxt2)
{
	// int qIndex_add = qIndex > ctxt2.qIndex ? qIndex : ctxt2.qIndex;
	// int num_mul_add = num_mul > ctxt2.num_mul ? num_mul : ctxt2.num_mul;
	// int num_relin_add = num_relin > ctxt2.num_relin ? num_relin : ctxt2.num_relin;
	// s-ar putea ca inainte sa fac operatiile sa trebuiasca sa aduc ctxt-urile la acelasi
	// nivel, plus inainte trebuie sa verific daca sunt criptate cu aceeasi cheie sk(0).
	// prepareCtxtForOperation(n, (*this), ctxt2);
	// ZZX c_add = n->AddModZZX(this->ctxt, ctxt2.ctxt, qIndex);
	// CtxtBit *ctxt_add = new CtxtBit(n, c_add);
	// ctxt_add -> num_mul = num_mul_add;
	// ctxt_add -> num_relin = num_relin_add;
	// ctxt_add -> qIndex = qIndex_add;
	// ctxt_add -> level = level;
	
	// pentru o singura cheie nu este nevoie de reliniarizare
	ZZX c_add = n->AddModZZX(ctxt, ctxt2.ctxt, 0);
	CtxtBit *ctxt_add = new CtxtBit(n);
	ctxt_add->setCtxt(c_add);
	
	return (*ctxt_add);
}
	
/*************************************************************************************/
CtxtBit& CtxtBit::operator+=(CtxtBit& ctxt2)
{
	(*this) = (*this) + ctxt2;
}
	
	
/*************************************************************************************/	
CtxtBit& CtxtBit::operator*(CtxtBit& ctxt2)
{
	/*int qIndex_mul = qIndex > ctxt2.qIndex ? qIndex : ctxt2.qIndex;
	int num_mul_mul = num_mul > ctxt2.num_mul ? num_mul : ctxt2.num_mul;
	int num_relin_mul = num_relin > ctxt2.num_relin ? num_relin : ctxt2.num_relin;*/
	
	// s-ar putea ca inainte sa fac operatiile sa trebuiasca sa aduc ctxt-urile la acelasi
	// nivel, plus inainte trebuie sa verific daca sunt criptate cu aceeasi cheie sk(0).
	// prepareCtxtForOperation(n, (*this), ctxt2);
	// inmultire modulo phi_m
	// ZZX c_mul = n->MulModZZX_XN1(this->ctxt, ctxt2.ctxt, qIndex);	
	// reliniarizare
	/*c_mul = n->Relin(c_mul, qIndex);
	c_mul = n->ModSwitch(c_mul, qIndex);
	n->CoeffReduce(c_mul, c_mul, qIndex+1);	
	CtxtBit *ctxt_mul = new CtxtBit(n, c_mul);
	ctxt_mul -> num_mul = num_mul_mul;
	ctxt_mul -> num_relin = num_relin_mul;
	ctxt_mul -> qIndex = qIndex_mul ;*/
	// ctxt_mul -> level = level + 1;
	
	ZZX c_mul = n->MulModZZX(ctxt, ctxt2.ctxt, 0);
	c_mul = n->Relin(c_mul, 0);
	n->ModSwitch(c_mul, 0);

	CtxtBit *ctxt_mul = new CtxtBit(n);
	ctxt_mul->setCtxt(c_mul);
	
	return (*ctxt_mul);
}
	
	
/*************************************************************************************/	
CtxtBit& CtxtBit::operator*=(CtxtBit& ctxt2)
{
	(*this) = (*this) * ctxt2;
}

/*************************************************************************************/	
/*void CtxtBit::encryptBit(int bit)
{
	qIndex = 0;
	ZZX b;
	b = bit;
	ctxt = n->Encrypt(b, qIndex);
}*/

/*************************************************************************************/	
int CtxtBit::decryptBit()const
{
	int bit = 0;
	ZZX bit_decrypted = n->Decrypt(ctxt, qIndex);
	
	if( bit_decrypted == 0 )
		bit = 0;
	else
		conv( bit, bit_decrypted[0]);
		
	return bit;
}

/*************************************************************************************/	
/*void CtxtBit::modSwitch()
{
	ctxt = n->ModSwitch(ctxt, qIndex);
	n->CoeffReduce(ctxt, ctxt, qIndex+1);
	qIndex++;
}*/
	
/*************************************************************************************/
/*void CtxtBit::reliniarizeCtxt()
{
	n->PolyCoeffReduce(ctxt, ctxt, qIndex);	
	ctxt = n->RelinRingFFT(ctxt, qIndex, 0);
	n->PolyCoeffReduce(ctxt, ctxt, qIndex);
	num_relin++;
}*/

/*void prepareCtxtForOperation(ntru *n, CtxtBit& ctxt1, CtxtBit& ctxt2)
{	
	for(int i=ctxt1.num_relin; i<ctxt2.num_relin; i++)
	{
		n->PolyCoeffReduce(ctxt1.ctxt, ctxt1.ctxt, ctxt1.qIndex);
	
		ctxt1.ctxt = n->RelinRingFFT(ctxt1.ctxt, ctxt1.qIndex, 0);
		
		n->PolyCoeffReduce(ctxt1.ctxt, ctxt1.ctxt, ctxt1.qIndex);

	}
	
	for(int i=ctxt2.num_relin; i<ctxt1.num_relin; i++)
	{
		n->PolyCoeffReduce(ctxt2.ctxt, ctxt2.ctxt, ctxt2.qIndex);
	
		ctxt2.ctxt = n->RelinRingFFT(ctxt2.ctxt, ctxt2.qIndex, 0);
		
		n->PolyCoeffReduce(ctxt2.ctxt, ctxt2.ctxt, ctxt2.qIndex);

	}
}*/