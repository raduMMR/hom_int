#ifndef HE_INTEGER_H
#define HE_INTEGER_H

#include <iostream>
#include <vector>
#include "Hom_NTRU_library/def.h"
#include "Hom_NTRU_library/ntru.h"
#include "Hom_NTRU_library/crt.h"
#include "Hom_NTRU_library/fft_mult.h"
#include "Hom_NTRU_library/general.h"
#include "CtxtBit.h"

using namespace std;

/* integer bitwise HE */
class HE_Integer
{
    int m_t_bits;
public:    
    HE_Integer(int t_bits);
    ~HE_Integer();

     void	encryptIntValue (ntru *n, vector<CtxtBit>& vct, int val);
     int 	decryptIntValue (ntru *n, const vector<CtxtBit>& vct);
	  
     void	encryptIntVector(ntru *n, vector<vector<CtxtBit> >& vvct, const vector<int>& values);
     void   decryptIntVector(ntru *n, vector<vector<CtxtBit> >& vvct, vector<int>& values);
     
     ///
     // void	encryptBit (ntru *n, CtxtBit& ct, int bit);
     // int	    decryptBit (ntru *n, const CtxtBit &ct);
     // int	    decryptBit (ntru *n, const CtxtBit &ct, int qIndex);
};

#endif // HE_INTEGER_H
