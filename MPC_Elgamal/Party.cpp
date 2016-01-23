#include "Party.h"
#include <NTL/ZZ.h>

Party::Party(void)
{
	if( Party::count_party == Party::N )
		throw new Party::NumberOfPartiesExceeded();

	this->party_number = Party::count_party;
	Party::count_party++;
}


Party::~Party(void)
{
	if( Party::count_party > 0 )
		Party::count_party--;
}

void Party::sample_polinom()
{
	ZZ_p::init(Party::q);
	NTL::random(a_i, N);
}

void Party::broadcast_g_s_i()const
{
	NTL::ZZ x;
	x = NTL::rep(a_i[0]);

	ZZ_p r;
	r = NTL::power(Party::alpha, x);

	Party::y = Party::y * r;
}

void Party::broadcast_g_s_ij(Party &party_j)const
{

}

void Party::compute_di(ZZ_p A)const
{
	A_x = A_x * power(A, x_i);
}
	
void Party::MPC_setup(int nb_of_parties)
{
	Party::N = nb_of_parties;

	// valori ales doar pentru test.
	Party::alpha = 3;
	Party::q = 7;
	Party::y = 1;
}


/*****************************************************************************/

void Party::encrypt(int message)const
{

}

void Party::decrypt(ZZ_p A, ZZ_p B)const
{

}

/*****************************************************************************/