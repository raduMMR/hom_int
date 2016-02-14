#include "test_multikey.h"

/*
run: build 
	./main < hNTRU/key.txt

build: 
	g++ -g hNTRU/crt.cpp hNTRU/fft_mult.cpp hNTRU/general.cpp hNTRU/ntru_multikey.cpp test_multikey.cpp main.cpp -o main -lntl
clean:
	rm -rf main
*/


void Test_NTRUMultikey_Comparison()
{
	ZZ q 	= to_ZZ("2");
	ZZ B 	= to_ZZ("1");
	
	GlobalParam gp;
	Set(gp, q, B, Modulus_M);

	NtruMultikey *ntru_context = new NtruMultikey(gp.q, gp.B, gp.N_PolyDegree, gp.M_CycDegree);
	
	ntru_context->IOReadModulus();
	ntru_context->IOReadModulusExt();

	ntru_context->ModulusFindRing(Num_Primes, Max_Prime, Dif_Prime, ntru_context->ReturnPolyMod());
	ntru_context->ComputeKeysRingRelin_FFT(Num_Primes, 1);
	
	ZZX m0;
	ZZX m1;
	ZZX result;
	
	int clear_result = 0;
	
	int bit0;
	int bit1;
	
	ZZX c0;
	ZZX c1;
	
	// operatia pe care o voi executa in simularea polinomului.
	int operation = 0;
	
	srand(time(NULL));
	
	int nr_polinoame = 2;
	int nr_operatii_polinom = 50;
	
	for(int i=0; i<nr_polinoame; i++)
	{
		clear_result = 0;
		m0 = 0;
		result = ntru_context->Encrypt(m0, 0);   // result = ENC(0)
		// cout<<"result = (" << ntru_context->Decrypt(result, 0);
		
		for(int j=0; j<nr_operatii_polinom; j++)
		{
			bit0 = rand() % 2;
			bit1 = rand() % 2;
			
			m0 = bit0;
			m1 = bit1;
			
			c0 = ntru_context->Encrypt(m0, 0);
			c1 = ntru_context->Encrypt(m1, 0);
			
			operation = ( operation + rand() ) % 2;
			
			if( operation == 0 )
			{
				// result = result + c0 + c1;
				//cout<< "+" << bit0 << "+" << bit1;
				
				result += ntru_context->AddModZZX(c0, c1, 0);
				
				clear_result = (clear_result + bit0 + bit1) % 2;
			}
			else
			{
				// result = result * c0 * c1
				// cout<< ")*" << bit0 << "*" << bit1;
				
				result= ntru_context->MulModZZX(result, c0, 0);
				result = ntru_context->Relin(result, 0);
				ntru_context->ModSwitch(result, 0);
				
				result = ntru_context->MulModZZX(result, c1, 0);
				result = ntru_context->Relin(result, 0);
				ntru_context->ModSwitch(result, 0);
				
				clear_result = clear_result * bit0 * bit1;
			} 
			
		}
		// cout<<" = " << ntru_context->Decrypt(result, 0) << " = " << clear_result << endl;
		
		if( ntru_context->Decrypt(result, 0) != clear_result )
			cout<< "Eroare la calculul polinomului de test nr. " << i << endl;
	
	}
	
	// @TODO : imbunatatirea testarii
	// se poate imbunatati testarea prin alegerea unui grad aleator cuprins intr-un
	// anumit interval si construiesc polinomul ca o suma de polinoame de un 
	// anumit grad
	
	delete ntru_context;
	ntru_context = NULL;
	
	std::cout << "\nTest incheiat\n";
}