#pragma once
#include "Party.h"

class TestElgamal
{
	static int nb_of_parties;

	static Party *parties;

	TestElgamal(void);

public:

	static void create_context(int nb);

	static void test_2_parties();

	static void test_encryption_decryption();
};

