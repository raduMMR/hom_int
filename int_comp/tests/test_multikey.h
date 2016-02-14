#include <time.h>
#include <stdlib.h>
#include <vector>
#include <assert.h>
#include <time.h>
#include <unistd.h>
#include "../hNTRU/def.h"
#include "../hNTRU/ntru_multikey.h"
#include "../hNTRU/crt.h"
#include "../hNTRU/fft_mult.h"
#include "../hNTRU/general.h"

/*
in aceasta functie sunt testate primitivele ce sunt 
folosite pentru comparatia intregilor si gasirea maximului
dintr-un vector de numere intregi criptate
*/
void Test_NTRUMultikey_Comparison();