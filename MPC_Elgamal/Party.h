#pragma once
#include <NTL/ZZ_pX.h>

using namespace NTL;

class Party
{
	int party_number;

	ZZ x_i;

	ZZ_pX a_i;


	// contextul MPC in care vom lucra.
	static int N;

	static int count_party;

	static ZZ_p y;

	static ZZ_p alpha;

	static ZZ q;

	static ZZ_p A_x;

	///

	void sample_polinom();

	void compute_di(ZZ_p A)const;

	// TBD : void compute_sigma_proof()const;

public:
	Party(void);

	~Party(void);

	void broadcast_g_s_i()const;

	void broadcast_g_s_ij(Party &party_j)const;
	
	static void MPC_setup(int nb_of_parties);

	// 

	void encrypt(int message)const;

	void decrypt(ZZ_p A, ZZ_p B)const;

	///

	// Exceptii aruncate de obiectele clasei Party.
	class PartyException{};

	class NumberOfPartiesExceeded:PartyException{};

	///
};

int Party::N = 0;
int Party::count_party = 0;
ZZ_p Party::y = ZZ_p(1);
ZZ_p Party::alpha = ZZ_p(0);
ZZ Party::q = ZZ(0);
ZZ_p Party::A_x = ZZ_p(1);

