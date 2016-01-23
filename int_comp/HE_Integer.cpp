#include "HE_Integer.h"
#include <NTL/ZZX.h>
#include <NTL/ZZ.h>

/*************************************************************************************/
HE_Integer::HE_Integer(int t_bits)
	: m_t_bits(t_bits)
{
	
}

HE_Integer::~HE_Integer()
{
	
}

/*************************************************************************************/
void HE_Integer::encryptIntValue (ntru *n, vector<CtxtBit>& vct, int val)
{
	unsigned int mask = 0x00000001;
    int bit;

	for (int i = 0; i < m_t_bits; i++)
    {
        bit = (val >> i) & mask;        	
		CtxtBit ct(n, bit);
        vct.push_back(ct);
    }  
}

/*************************************************************************************/
int HE_Integer::decryptIntValue (ntru *n, const vector<CtxtBit>& vct)
{
	int val = 0;

 	for (int i = 0; i < m_t_bits /* vct.size()*/; i++)
	{
		int bit = vct[i].decryptBit();
		val += bit * pow (2, i);
	}
	
	return val;
}

/*************************************************************************************/	  
void HE_Integer::encryptIntVector(ntru *n, vector<vector<CtxtBit> >& vvct, const vector<int>& values)
{
	vector<CtxtBit> vct;
	
    for (int i = 0; i < values.size(); i++)
    {
        vct.clear();
        encryptIntValue(n, vct, values[i]);		
        vvct.push_back(vct);
    }    
}

/*************************************************************************************/	  
 void HE_Integer::decryptIntVector(ntru *n, vector<vector<CtxtBit> >& vvct, vector<int>& values)
 {
	 int intreg = 0;
	 
	 for (int i = 0; i < vvct.size(); i++)
     {
        intreg = decryptIntValue(n, vvct[i]);		
        values.push_back(intreg);
     } 
 }

/*************************************************************************************/	  
/*int	HE_Integer::decryptBit(ntru *n, const CtxtBit &ct, int qIndex)
{
	ZZX bit_decrypted = m_n->Decrypt_XN1(ct, qIndex);
	
	int bit;

	if( bit_decrypted == 0 )
		bit = 0;
	else
		conv( bit, bit_decrypted[0]);
	
	return bit;
}*/

/*************************************************************************************/ 
/*int	HE_Integer::decryptBit (ntru *n, const CtxtBit &ct)
{	
	ZZX bit_decrypted = m_n->Decrypt_XN1(ct, 0);
	
	int bit;
	
	// daca polinomul rezultat in urma decriptarii este 0
	// clar bit-ul este zero. Functia conv(poly_dest, poly_src) apelata cu un polinom 0
	// da eroare la rulare pentru ca bit_decrypted este egal cu [].
	if( bit_decrypted == 0 )
		bit = 0;
	else
		conv( bit, bit_decrypted[0]);
	
	return bit;
}*/
