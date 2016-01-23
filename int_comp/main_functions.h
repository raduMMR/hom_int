#include "def.h"
#include "ntru.h"
#include "crt.h"
#include "fft_mult.h"
#include "general.h"
#include "fheaes.h"
using namespace std;

#include <assert.h>

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
	
	n->IOReadModulus();
	n->IOReadModulusExt();

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
		{
			r = n->MulModZZX(a, b, 0);
			a = r;
		}
	t.Stop();
	t.ShowTime("Time is:\t");
}



/*ZZX CreateMess(int size, int q){
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
*/





