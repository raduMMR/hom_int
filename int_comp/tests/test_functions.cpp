#include "test_functions.h"

/*
run: build 
	./main < hNTRU/key.txt

build: 
	g++ -g hNTRU/crt.cpp hNTRU/fft_mult.cpp hNTRU/general.cpp hNTRU/ntru_multikey.cpp hNTRU/fheaes.cpp hcomp/CtxtBit.cpp hcomp/HE_Integer.cpp hcomp/HE_Utils.cpp test_functions.cpp main.cpp -o main -lntl
clean:
	rm -rf main
*/


// ctxt_of_1 = ENC(1)
CtxtBit *ctxt_of_1;			

// @TODO: shared pointer !!!
ntru*    ntru_context;

void ntru_tutorial()
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
	
	ZZX c0;
	ZZX c1;
	ZZX r;
	
	int bit0 = 0;
	int bit1 = 0;
	m0 = bit0;
	m1 = bit1;	
		
	c0 = ntru_context->Encrypt(m0, 0);
	c1 = ntru_context->Encrypt(m1, 0);
	
	result= ntru_context->MulModZZX(c0, c1, 0);
	result = ntru_context->Relin(result, 0);
	result = ntru_context->ModSwitch(result, 0);
	
	c1 = ntru_context->ModSwitch(c1, 0);

	result= ntru_context->MulModZZX(result, c1, 1);
	result = ntru_context->Relin(result, 0);
	result = ntru_context->ModSwitch(result, 1);
	
	cout << "Decrypt = " << ntru_context->Decrypt(result, 0);
	
	delete ntru_context;
}


//////////////////////////////////////////////////////////////////////

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


/****************     Functiile de test ale lui Yarkin       **********************************/

//////////////////////////////KEY GENERATION//////////////////////////////////
void KeyGenOut(){
	ZZ q 	= to_ZZ("2");
	ZZ B 	= to_ZZ("1");

	GlobalParam gp;
	Set(gp, q, B, Modulus_M);

	myCRT c;

	ntru n(gp.q, gp.B, gp.N_PolyDegree, gp.M_CycDegree);
	n.SetModulus();

	c.SetModulus(n.ReturnPolyMod());
	c.ComputeFactors(gp.D_FactDegree, gp.FactSize);
	c.CalculateMs();
	c.CalculateNs();
	c.CalculateMxNs();

	cout << n.ReturnPolyMod() << endl;
	cout << n.ReturnPolyMod_M_1() << endl;
/*
	cout << c.Returnmodulus() << endl;
	cout << c.Returnsize() << endl;
	vec_ZZ_pX temp;
	temp = c.Returnfactors();

	for(int i=0; i<c.Returnsize(); i++)
		cout << temp[i] << endl;

	temp = c.ReturnM();
	for(int i=0; i<c.Returnsize(); i++)
			cout << temp[i] << endl;

	temp = c.ReturnN();
	for(int i=0; i<c.Returnsize(); i++)
			cout << temp[i] << endl;

	temp = c.ReturnMxN();
	for(int i=0; i<c.Returnsize(); i++)
			cout << temp[i] << endl;
*/

}

void FHETry(){
	myTimer t;
	FheAES f;

	t.Start();
	f.LTVKeySetUp();
	t.Stop();
	t.ShowTime("Set Time:\t");

	unsigned char m[16] = {
			0x00, 0x44, 0x88, 0xCC,
			0x11, 0x55, 0x99, 0xDD,
			0x22, 0x66, 0xAA, 0xEE,
			0x33, 0x77, 0xBB, 0xFF
	};
	f.SetMessageBits(m);
			cout << "Encrypt Message Start" << endl;



	t.Start();
	f.EncryptMessage();
	t.Stop();
	t.ShowTime("Enc Mess Time:\t");

	cout << "Encrypt Message Done" << endl;

	t.Start();
	f.SetKeys();
	t.Stop();
	t.ShowTime("Key Set Time:\t");
	cout << "Set Keys Done" << endl;

	f.AESEncryption();

	cout << dec << f.num_mul << endl;
	cout << dec << f.num_relin << endl;

}


void Test(){

	ZZ q 	= to_ZZ("2");
	ZZ B 	= to_ZZ("1");
	GlobalParam gp;
	Set(gp, q, B, Modulus_M);

	ntru *n = new ntru(gp.q, gp.B, gp.N_PolyDegree, gp.M_CycDegree);
	myCRT c;
	n->IOReadModulus();
	n->IOReadModulusExt();
	c.IOReadAll();

	n->ModulusFindRing(Num_Primes, Max_Prime, Dif_Prime, n->ReturnPolyMod());
	n->ComputeKeysRingRelin_FFT(Num_Primes, 1);


	myTimer t;
	ZZX m0, m1;
	m0 = 1;
	m1 = 1;
	ZZX a = n->Encrypt(m0, 0);
	ZZX b = n->Encrypt(m1, 0);
	ZZX r;

	t.Start();
		for(int i=0; i<10; i++)
			r = n->MulModZZX(a, b, 0);
	t.Stop();
	t.ShowTime("Time is:\t");

}



ZZX CreateMess(int size, int q){
	ZZX x;
	ZZ r = to_ZZ(time(NULL));
	SetSeed(r);

	ZZ l;
	ZZ one = to_ZZ("1");
	for(int i=0; i<size; i++){
		l = RandomBnd(one<<q);
		SetCoeff(x, i, l);
	}
	return x;
}


void readTime(){

	int N 		= 32768;
	int log_q	= 1271;
	ZZX t;

	ZZX a[log_q];

	for(int i=0; i<log_q; i++)
		 a[i] = CreateMess(N, log_q);

	myTimer mt;
	mt.Start();
		for(int i=0; i<log_q; i++)
			cout << a[i] << endl;
	mt.Stop();
	mt.ShowTime("Write:\t");

}

void MemTest(){

	int N = 32768;
//	int q = 1024;
	int l = 256;

	ZZX t[l];
	for(int i=0; i<l; i++){
		for(int j=0; j<N; j++){
			SetCoeff(t[i], j, 1);
		}
	}
	cout << "Sleep" << endl;
	sleep(10);
}

void TT(){

	MemTest();
	cout << "Out Func" << endl;
	sleep(10);
	MemTest();
	cout << "Out Func" << endl;
	sleep(10);
	MemTest();
	cout << "Out Func" << endl;
	sleep(10);


}

void DisplayRealTime(){
	time_t e = time(0);   // get time now
	struct tm * now = localtime( & e );
	cout << now->tm_hour << ":" << now->tm_min << ":" <<  now->tm_sec << endl;
}

void mm(ZZX &out, ZZX &in0, ZZX &in1){
	out = in0*in1;
}

void TimeTest(){
	int NN = 32768/2/2;
	int NN2 = 32768;

	ZZ q = GenPrime_ZZ(1271, 80);
	ZZX a,b;

	for(int i=0; i<NN; i++)
		SetCoeff(a, i, RandomBnd(q));

	for(int i=0; i<NN2; i++)
		SetCoeff(b, i, RandomBnd(q));


	ZZX c;


	myTimer t;
	t.Start();
	for(int i=0; i<10; i++)
		mm(c, a,b);
	t.Stop();
	t.ShowTime("Time:\t");

}
