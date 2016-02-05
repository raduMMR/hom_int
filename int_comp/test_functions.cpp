#include "test_functions.h"

// ctxt_of_1 = ENC(1)
CtxtBit *ctxt_of_1;			

// @TODO: shared pointer !!!
ntru*    ntru_context;

void Test_NTRU_Comparison(){
	
	ZZ q 	= to_ZZ("2");
	ZZ B 	= to_ZZ("1");
	
	GlobalParam gp;
	Set(gp, q, B, Modulus_M);

	ntru_context = new ntru(gp.q, gp.B, gp.N_PolyDegree, gp.M_CycDegree);
	
	ntru_context->IOReadModulus();
	ntru_context->IOReadModulusExt();

	ntru_context->ModulusFindRing(Num_Primes, Max_Prime, Dif_Prime, ntru_context->ReturnPolyMod());
	ntru_context->ComputeKeysRingRelin_FFT(Num_Primes, 1);
	
	ctxt_of_1 = new CtxtBit(ntru_context, 1);
	
	vector<CtxtBit> ctxt_x;
	vector<CtxtBit> ctxt_y;
	
	
	HE_Integer he_integer(3);
		
	he_integer.encryptIntValue(ntru_context, ctxt_x, 2);
	he_integer.encryptIntValue(ntru_context, ctxt_y, 5);
		
		
	evaluate_X_eq_Y(ctxt_x, ctxt_y, 3);
	
	
	
// cleanup
	delete ctxt_of_1;
	ctxt_of_1 = NULL;
	delete ntru_context;
	ntru_context = NULL;
}





/*******************   Testarea unor functionalitati ale bibliotecii NTRU *********************/	  

/*
in aceasta functie construiesc un polinom cu operatii aleatoare si
valori ale variabilelor alese aleator la fiecare pas,
fiecare polinom cu un grad ales la testare
in timp ce evaluez homomorfic polinomul, calculez
polinom si in clar iar la afarsitul unui ciclu compar 
cele doua valori pentru a testa daca evaluarea a fost
efectuata corect.
*/
void Test_NTRU(int nr_polinoame, int nr_operatii_polinom)
{
	ZZ q 	= to_ZZ("2");
	ZZ B 	= to_ZZ("1");
	
	GlobalParam gp;
	Set(gp, q, B, Modulus_M);

	ntru_context = new ntru(gp.q, gp.B, gp.N_PolyDegree, gp.M_CycDegree);
	
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
}


void test_ntru_multiplication_depth()
{
	
	ZZ q 	= to_ZZ("2");
	ZZ B 	= to_ZZ("1");
	
	GlobalParam gp;
	Set(gp, q, B, Modulus_M);

	ntru *ntru_context = new ntru(gp.q, gp.B, gp.N_PolyDegree, gp.M_CycDegree);
	
	ntru_context->IOReadModulus();
	ntru_context->IOReadModulusExt();

	ntru_context->ModulusFindRing(Num_Primes, Max_Prime, Dif_Prime, ntru_context->ReturnPolyMod());
	ntru_context->ComputeKeysRingRelin_FFT(Num_Primes, 1);
	
	ZZX m0;
	ZZX m1;
	
	ZZX result;
	m0 = 1;
	result = ntru_context->Encrypt(m0, 0); // result = ENC(1)
	int clear_result = 1;
	
	int bit0 = 1;
	int bit1 = 0;
	ZZX c0;
	ZZX c1;
	
	m0 = bit0;
	m1 = bit1;
	c0 = ntru_context->Encrypt(m0, 0);
	c1 = ntru_context->Encrypt(m1, 0);
	
	int level = 0;
	
	int max_mult_depth = 0;
	
	for(int i=0; i<1000; i++)
	{
		result= ntru_context->MulModZZX(result, c0, level);
		result = ntru_context->Relin(result, level);
		ntru_context->ModSwitch(result, level);
		
		// level++;
		
		result = ntru_context->MulModZZX(result, c1, level);
		result = ntru_context->Relin(result, level);
		ntru_context->ModSwitch(result, level);
		
		// level++;
				
		clear_result = clear_result * bit0 * bit1;
		
		if( ntru_context->Decrypt(result, 1) != clear_result )
			break;
		max_mult_depth++;
	}

	cout<<"Schema accepta o adancime multiplicativa de " << max_mult_depth << " niveluri.\n";

	delete ntru_context;
}