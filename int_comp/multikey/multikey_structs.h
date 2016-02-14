#include "NTL/ZZX.h"

/*
@brief fiecare key este o cheie de evaluare pentru un level ek_t
*/
struct eval_key{
	ZZX *key;
};


/*
@brief structura pentru descrierea unui polinom de catre party ( client ) 
	un polinom este o reprezentare a unei functiei delegate pentru calcul cloud-ului
														---
	P(x_i) = Suma( coeff * Produs( x_i^(exp_i) ) ), i = 1,N.
	
*/
struct PolyDescription
{
	int nb_of_levels;		// pentru suma
	
	struct PolyLevel
	{
		int party_index;		// x_i
		ZZ party_exponent;		// exp_i
	};
																			//		    ---
	// reprezentarea polinoamelor de forma coeff * produs ( x_i ^ ( exp_i ) ),  pt. i = 1,N.
	PolyLevel *level_poly;		// vector de polinoame pt. fiecare produs intermediar
	ZZ *level_coeff;			// coeficientul corespunzator produsului intermediar
};


/*
@brief structura care reprezinta un bit criptat dat ca intrare
	   cloud-ului intr-un context multiparti
*/
struct PartyEncBit
{
	ZZX ctxt;
	int *party_indexes;			// retine indecsii partilor implicate in calculul ctxt-ului pt. dec
	int nb_of_parties;
};

/*
@brief structura care reprezinta un intreg criptat dat ca intrare
		cloud-ului intr-un context multiparti
*/
struct PartyEncInt
{
	PartyEncBit *int_bits; 
	int nb_of_bits;	
};

/*
@brief structura prin care o parte isi transmite cloud-ului intrarea
*/
struct HomEvalInput
{
	PartyEncBit *enc_bits;
	ZZ nb_of_bits;
	int *party_indexes;
	int nb_of_parties;
};

/*
@brief structura prin care se reprezinta rezultatul functiei de evaluare homomorfica
@member PartyEncBit * bitii rezultati in urma operatiei de evaluare
@member nb_of_bits numarul de biti al rezultatului
*/
struct HomEvalOutput
{
	PartyEncBit *enc_bits;
	ZZ nb_of_bits;
	int *party_indexes;	
	int nb_of_parties;
};