#include <unistd.h>
#include <sys/sysinfo.h>

#include "main_functions.h"
#include "test_comparison.h"

int main()
{
	// Need to compute and extract the key first. Run only KeyGenOut as: ./main > key.txt
	// KeyGenOut();
	
	// Later command KeyGenOut and command out FHETry. Compile and run as: ./main < key.txt
	// FHETry();
	
	
	// testez functiile de baza pentru comparatie
	// In lucru
	Test_NTRU_Comparison();
	
	
	
	// aici am testat multiplicarea si diverse polinoame pentru a 
	// verifica corectitudinea bibliotecii.
	
	// test_ntru_multiplication_depth();
	// Test_NTRU(10, 80);
	
	
	
	return 0;
}











