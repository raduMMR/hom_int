#include "Party.h"
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>

// random_zz_pX
// numarul de parti implicate in computatie.
int Party::N = 0;

// un indice pentru partea care este creata la un moment dat.
int Party::count_party = 1;

// cheia publica Elgamal comuna.
ZZ_p Party::y = ZZ_p(1);

// generatorul grupului folosit pentru schema de criptare.
ZZ_p Party::alpha = ZZ_p(0);
ZZ Party::alpha_ = ZZ(0);

// modulul grupului, alpha^n = 1 mod q.
ZZ Party::q = ZZ(11);

// ordinul grupului 
ZZ Party::n = ZZ(47);

// cheia secreta comuna, care va fi calculata la afarsit din
// cheile secrete de verificare ale partilor, adica din x_i.
// ZZ_p Party::A_x = ZZ_p(1);



Party::Party(void)
{
	if( Party::N == 0 )
		throw new Party::PartyException("Constructorul clasei Party\nNumarul de parti este 0.");

	if( Party::count_party > Party::N )
		throw new Party::
		PartyException("Constructorul clasei Party\nS-a atins limita maxima de parti admise in context");

	this->party_number = Party::count_party;
	Party::count_party++;

	sample_polynomial();

	// ZZ_p i = ZZ_p(this->party_number);
	// x_i = eval(a_i, i);
	x_i = 0;

	broadcast_g_s_i();
}

void Party::MPC_setup(int nb_of_parties)
{
	Party::N = nb_of_parties;
	// Party::q = 283;
	Party::q = 11;
	ZZ_p::init(q);
	// valori alese doar pentru test.
	// Party::alpha = ZZ_p(60);
	// Party::alpha_ = ZZ(60);
	Party::alpha = ZZ_p(2);
	Party::alpha_ = ZZ(2);
	Party::n = 47;
}



Party::~Party(void)
{
}

/***************************************************************************/
/*
@esantonioneaza un polinom cu grad < N, nr. de parti
*/
void Party::sample_polynomial()
{
	// ZZ_p::init(Party::q);

	NTL::random(a_i, N);

	// std::cout<< "a_" << this->party_number << " = " << a_i << std::endl;
}



/*
@brief folosita pentru calculul cheii secrete comune
@param1 A ciphertextul de decriptat este (A,B)
@return calculeaza d_i = A ^ x_i 
*/
ZZ_p Party::compute_di(ZZ_p A)const
{
	ZZ_p lambda = ZZ_p(1);

	for(int j=1; j<=Party::N; j++)
	{
		if( j == this->party_number )
			continue;
		lambda *= ( j /( j - this->party_number ) );
	}

	/*ZZ e;
	conv(e, x_i);
	ZZ_p di = power(A, e);*/

	lambda *= x_i;

	ZZ l;
	conv(l, lambda);

	ZZ_p di = power(A, l);

	// std::cout<< "d_" << this->party_number << " = " << di << std::endl;

	return di;
}


/***************************************************************************/

/*
prin publicarea g^(s_i) se calculeaza cheia publica comuna
y = Produs( g^(s_i) ).
*/
void Party::broadcast_g_s_i()const
{
	NTL::ZZ s_i;
	s_i = NTL::rep(a_i[0]);

	ZZ_p r;
	r = NTL::power(Party::alpha, s_i);

	Party::y = Party::y * r;
}




/*
@brief fiecare Parte imparte celorlalte parti polinomul sau secret
	printr-un protocol Verifiable Secret Sharing.
	Aceste "share-uri" vor fi folosite pentru constructia cheii secrete
	a fiecarei parti, cheie secreta care va combinata cu celelalte chei 
	secrete pentru a calcula cheia secreta comuna.
@param party_j Partea catre care este distribuit share-ul
*/
void Party::broadcast_g_s_ij(Party &party_j)const
{
	ZZ_p j = ZZ_p(party_j.party_number);

	ZZ_p s_ij;

	s_ij = eval(a_i, j);

	party_j.x_i += s_ij;

	// std::cout<< "x_" << party_j.party_number << " = " << party_j.x_i << std::endl;
}



/*****************************************************************************/
/*
@brief criptare Elgamal
*/
void Party::encrypt(ZZ m, ZZ &c1, ZZ &c2)const
{
	int len = 16;
	ZZ k;
	RandomBits(k, len);
	k = k % Party::q;

	c1 = PowerMod(alpha_, k, q);

	ZZ K;
	ZZ y_;

	conv(y_, y);

	PowerMod(K, y_, k, q);

	MulMod(c2, K, m, q);
	 
}

/*
decriptarea "threshold" a ciphertextului (A,B) unde 
			m  = B / produs( A ^ x_i ).
*/
void Party::decrypt(ZZ &m, ZZ A, ZZ B, Party *otherParties)const
{
	ZZ_p A_x;
	ZZ_p A_p = ZZ_p(1);
	ZZ pi_A_x(1);
	ZZ inv_A;

	conv(A_p, A);

	A_x = compute_di(A_p);

	for(int i=1; i<=Party::N; i++)
	{
		if( i == this->party_number )
		{
			continue;
		}

		A_x *= otherParties[i-1].compute_di(A_p);
	}

	conv(pi_A_x, A_x);
	inv_A = InvMod(pi_A_x, Party::q);

	// conv(joint_A, A_x);
	// std::cout<< "joint secret key = " << joint_A << std::endl;

	m = MulMod(inv_A, B, Party::q);
}

/*****************************************************************************/