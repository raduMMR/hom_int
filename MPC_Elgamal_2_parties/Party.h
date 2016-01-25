#pragma once
#include <NTL/ZZ_pX.h>

using namespace NTL;

/*
@brief 1 obiect Party = 1 Participant la protocolul MPC bazat pe Threshold Elgamal.
*/
class Party
{
	int party_number;

	ZZ_p x_i;

	ZZ_pX a_i;

	ZZ_p *di;

	///

	static int N;

	static int count_party;

	static ZZ q;

	static ZZ n;

	static ZZ_p alpha;
	static ZZ alpha_;
	
	static ZZ_p y;

	// static ZZ_p A_x;

	///

	void sample_polynomial();

	void broadcast_g_s_i()const;

	ZZ_p compute_di(ZZ_p A)const;

	// TBD : void compute_sigma_proof()const;

public:
	Party(void);

	~Party(void);

	void broadcast_g_s_ij(Party &party_j)const;
	
	static void MPC_setup(int nb_of_parties);

	/// 

	void encrypt(ZZ m, ZZ &c1, ZZ &c2)const;

	void decrypt(ZZ& m, ZZ A, ZZ B, Party *otherParties)const;

	///

	// Exceptii aruncate de obiectele clasei Party.
	class PartyException
	{
		char exceptie[300];
	public:
		PartyException(char *ex)
		{
			strcpy(exceptie, ex);
		}

		void printExceptie()const
		{
			std::cout<<"Exceptie : clasa Party "<<exceptie<<std::endl;
		}
	};

	///

	/// test purpose only
	void static print_public_key()
	{
		std::cout<<"Cheia publica = " << y << std::endl;
	}
	///
};