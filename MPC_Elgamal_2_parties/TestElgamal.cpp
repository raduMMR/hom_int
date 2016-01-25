#include "TestElgamal.h"

using namespace std;

Party *TestElgamal::parties = NULL;
int TestElgamal::nb_of_parties = 0;

void TestElgamal::create_context(int nb)
{
	nb_of_parties = nb;

	Party::MPC_setup(nb);

	parties = new Party[nb];

	for(int i=0; i<nb; i++)
	{
		for(int j=0; j<nb; j++)
		{
			parties[i].broadcast_g_s_ij(parties[j]);
		}
	}
}

void TestElgamal::test_2_parties()
{
	try
	{
		Party::MPC_setup(2);

		Party p[2];

		ZZ m1 = ZZ(2);
		ZZ m2 = ZZ(2000);
		ZZ c1, c2;

		p[0].broadcast_g_s_ij(p[0]);
		p[0].broadcast_g_s_ij(p[1]);

		p[1].broadcast_g_s_ij(p[0]);
		p[1].broadcast_g_s_ij(p[1]);

		// Party::print_public_key();

		for(int i=0; i<1000; i++)
		{
			m1 = ( i % 11 );

			p[0].encrypt(m1, c1, c2);

			// cout<< "c1 = " << c1 << " c2 = " << c2 << endl;

			p[0].decrypt(m2, c1, c2, p);

			if( m1 == m2 )
			{
				// cout<<"Cryptosistemul functioneaza.\n";
			}
			else
			{
				cout<<"Eroare la decriptare la iteratia " << i << "\n";
			}
			// cout<<"m2 = "<< m2 << endl;
		}
	}
	catch(Party::PartyException *pe)
	{
		cout<<"!!!!!!!!!!!Exceptie!!!!!!!!!!!!!!\n";
		pe->printExceptie();
		delete pe;
	}
}

void TestElgamal::test_encryption_decryption()
{
	ZZ m1 = ZZ(2);
	ZZ m2 = ZZ(3);
	ZZ c1, c2;

	try
	{
		for(int i=0; i<1000; i++)
		{
			m1 = ( i % 11 );

			parties[0].encrypt(m1, c1, c2);

			// cout<< "c1 = " << c1 << " c2 = " << c2 << endl;

			parties[0].decrypt(m2, c1, c2, parties);

			if( m1 == m2 )
			{
				// cout<<"Cryptosistemul functioneaza.\n";
			}
			else
			{
				cout<<"Eroare la decriptare la iteratia " << i << "\n";
			}
			// cout<<"m2 = "<< m2 << endl;
		}
	}
	catch(Party::PartyException *pe)
	{
		cout<<"!!!!!!!!!!!Exceptie!!!!!!!!!!!!!!\n";
		pe->printExceptie();
		delete pe;
	}

	cout<<"Final test_encryption_decryption\n";
}
