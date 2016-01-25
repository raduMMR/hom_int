//#include<iostream>
//#include <NTL\ZZ.h>
//
//using namespace std;
//using namespace NTL;
//
//class Elgamal
//{
//	ZZ q;
//	ZZ n;
//	ZZ alpha;
//	ZZ y_a;
//	ZZ x_a;
//
//	void ElgamalSetup()
//	{
//		q = 283;
//		n = 47;
//		alpha = 60;
//		x_a = 7;
//		PowerMod(y_a, alpha, x_a, q);
//	}
//
//public:
//
//	Elgamal()
//	{
//		ElgamalSetup();
//	}
//
//	void Encrypt(ZZ M, ZZ &c1, ZZ &c2)
//	{
//		long len = 16;
//		ZZ k;
//		RandomBits(k, len);
//		PowerMod(c1, alpha, k, q);
//		ZZ K;
//		PowerMod(K, y_a, k, q);
//		MulMod(c2, K, M, q);
//	}
//
//	void Decrypt(ZZ &M, ZZ c1, ZZ c2)
//	{
//		ZZ K;
//		PowerMod(K, c1, x_a, q);
//		ZZ K_1;
//		InvMod(K_1, K, q);
//		MulMod(M, c2, K_1, q);
//	}
//
//};
//
//int main()
//{
//	Elgamal elgamal;
//	ZZ m, c1, c2, mp;
//
//	elgamal.Encrypt(m, c1, c2);
//	elgamal.Decrypt(mp, c1, c2);
//
//	if( m == mp )
//	{
//		cout<<"Cryptosystemul functioneaza, pana acum ;)";
//	}
//	else
//	{
//		cout<<"Eroare in cryptosystem.\n";
//	}
//
//
//	cout<<"Test Elgamal terminat.\n";
//	return 0;
//}


