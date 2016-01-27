#ifndef DEF_H_
#define DEF_H_

//////FOR PHI OF M
#define M_SIZE 51

#if M_SIZE == 51			//Toy setting
	#define Modulus_M	51
	#define Factor_Deg	8
#elif M_SIZE == 8191
	#define Modulus_M	8191
	#define Factor_Deg	13
#elif M_SIZE == 21845
	#define Modulus_M	21845
	#define Factor_Deg	16
#elif M_SIZE == 32767
	#define Modulus_M	32767
	#define Factor_Deg	16
#elif M_SIZE == 65535
	#define Modulus_M	65535
	#define Factor_Deg	16
#endif

////FOR MODULUS
#define Num_Primes  41
#define Max_Prime	(1271)
#define Dif_Prime   31

#endif /* DEF_H_ */
