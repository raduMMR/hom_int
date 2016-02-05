#include <time.h>
#include <stdlib.h>
#include <vector>
#include <assert.h>
#include <time.h>
#include "hNTRU/def.h"
#include "hNTRU/ntru.h"
#include "hNTRU/crt.h"
#include "hNTRU/fft_mult.h"
#include "hNTRU/general.h"
#include "hcomp/HE_Utils.h"
#include "hcomp/HE_Integer.h"

using namespace std;

/*
in aceasta functie sunt testate primitivele ce sunt 
folosite pentru comparatia intregilor si gasirea maximului
dintr-un vector de numere intregi criptate
*/
void Test_NTRU_Comparison();

/*
in aceasta functie construiesc un polinom cu operatii aleatoare si
valori ale variabilelor alese aleator la fiecare pas,
fiecare polinom cu un grad ales la testare
in timp ce evaluez homomorfic polinomul, calculez
polinom si in clar iar la afarsitul unui ciclu compar 
cele doua valori pentru a testa daca evaluarea a fost
efectuata corect.
*/
void Test_NTRU(int nr_polinoame, int nr_operatii_polinom);

void test_ntru_multiplication_depth();