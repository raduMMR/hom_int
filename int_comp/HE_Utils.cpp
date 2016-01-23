#include "HE_Utils.h"
#include "HE_Integer.h"

extern CtxtBit*		ctxt_of_1;                  // encryption of 1 (constant)
extern int		    t_bits;
extern ntru*        ntru_context;

static clock_t  	t1, t2, t3;

using namespace std;

/*************************************************************************************/
float clock_diff(clock_t &t1, clock_t &t2)
{
	return (((float)t2 - (float)t1) / CLOCKS_PER_SEC); // - seconds
}
/*************************************************************************************/


/*************************************************************************************/
CtxtBit compute_z (int i, int j, vector<CtxtBit>& ct_x, vector<CtxtBit>& ct_y)
{
	assert (ct_x.size() > 0 && ct_y.size() > 0);
			
	if (j == 1)
	{
		CtxtBit ret ( (*ctxt_of_1) );				
		ret += ct_x[i];
		ret += ct_y[i];
		
		return ret;
	}

	int l;
	l = (j%2 == 0) ? j/2: j/2 + 1;

	CtxtBit ret = compute_z(i+l, j-l, ct_x, ct_y);
	CtxtBit ct = compute_z (i, l, ct_x, ct_y);
    
	// ret.multiplyBy(ct);                // ret *= ct;	
	ret *= ct;
	
	return ret;	
}

/*************************************************************************************/
CtxtBit compute_t (int i, int j, vector<CtxtBit>& ct_x, vector<CtxtBit>& ct_y)
{
	assert (ct_x.size() > 0 && ct_y.size() > 0);	
	if (j == 1)
	{	
		CtxtBit ret (ct_x[i]);
		ret *= ct_y[i];
		ret += ct_x[i];

		return ret;
	}
	
	int l;
	l = (j%2 == 0) ? j/2: j/2 + 1;

	CtxtBit ret = compute_t(i+l, j-l, ct_x, ct_y);
	CtxtBit ct_z = compute_z (i+l, j-l, ct_x, ct_y);
	CtxtBit ct_t = compute_t (i, l, ct_x, ct_y);

	ct_z *= ct_t;
	ret += ct_z;
	
	return ret;	
}

/*************************************************************************************/
CtxtBit compute_s (int i, int j, vector<CtxtBit>& ct_x, vector<CtxtBit>& ct_y)
{
	assert (ct_x.size() > 0 && ct_y.size() > 0);	
	if (j == 1)
	{				
		CtxtBit ret (ct_x[i]);
		ret *= ct_y[i];
		ret += ct_y[i];
		ret += *ctxt_of_1;
		
		return ret;
	}

	int l;
	l = (j%2 == 0) ? j/2: j/2 + 1;

	CtxtBit ret = compute_t(i+l, j-l, ct_x, ct_y);
	CtxtBit ct_z = compute_z (i+l, j-l, ct_x, ct_y);
	CtxtBit ct_s = compute_s (i, l, ct_x, ct_y);

	ct_z *= ct_s;
	ret += ct_z;

	return ret;
}

/*************************************************************************************/
void evaluate_X_gt_Y (vector<CtxtBit>& ct_x, vector<CtxtBit>& ct_y, int t_bits)
{
	// clock_t t1, t2;	
	HE_Integer heInt (t_bits);
	
	cout << endl << "Evaluating (X>Y)...." << std::flush;
 	cout << "Evaluating " << heInt.decryptIntValue(ntru_context, ct_x) << " > " << heInt.decryptIntValue(ntru_context, ct_y) << " compute_t: " << std::flush;
	// t1 = clock();
	CtxtBit ct_t = compute_t (0, t_bits, ct_x, ct_y);
	// t2 = clock();		
	
	int dec_t = ct_t.decryptBit();
	
	if( dec_t == 1 )
	{
		cout<<" X gt Y\n";
	}
	else
	{
		if(dec_t == 0)
		{
			cout<<"X not(gt) Y\n";
		}
		else
		{
			cout<<"Eroare: valoarea lui dec_t este " << dec_t <<endl;
		}
	}
	
 	// cout << endl<< dec_t << " (" << clock_diff(t1,t2)<< "ms)" << endl << std::flush;
}

/*************************************************************************************/
void evaluate_X_ge_Y (vector<CtxtBit>& ct_x, vector<CtxtBit>& ct_y, int t_bits)
{
	// clock_t t1, t2;
	HE_Integer heInt (t_bits);
	
	cout << endl << "Evaluating (X>=Y)..." << std::flush;
	cout << "Evaluating " << heInt.decryptIntValue(ntru_context, ct_x) << " >= " << heInt.decryptIntValue(ntru_context, ct_y) << " compute_s: " << std::flush;
	// t1 = clock();	
	CtxtBit ct_s = compute_s (0, t_bits, ct_x, ct_y);
	// t2 = clock();
	
	int dec_s = ct_s.decryptBit();
	
	if( dec_s == 1 )
	{
		cout<<" X ge Y\n";
	}
	else
	{
		if(dec_s == 0)
		{
			cout<<"X lt Y\n";
		}
		else
		{
			cout<<"Eroare: valoarea lui dec_s este " << dec_s <<endl;
		}
	}
	
	// cout << endl<< dec_s << " (" << clock_diff(t1,t2)<< "ms)" << endl << std::flush;
}

/*************************************************************************************/
void evaluate_X_eq_Y(vector<CtxtBit>& ct_x, vector<CtxtBit>& ct_y, int t_bits)
{
	// clock_t t1, t2;
	HE_Integer heInt (t_bits);
	
	cout << endl << "Evaluating (X=Y)...." << std::flush;
	cout << "Evaluating " << heInt.decryptIntValue(ntru_context, ct_x) << " = " << heInt.decryptIntValue(ntru_context, ct_y) << " compute_z: " << std::flush;
	// t1 = clock();
	CtxtBit ct_z = compute_z (0, t_bits, ct_x, ct_y);
	// t2 = clock();	
		
	int dec_z = ct_z.decryptBit();
	
	if( dec_z == 1 )
	{
		cout<<" X eq Y\n";
	}
	else
	{
		if(dec_z == 0)
		{
			cout<<"X eq Y\n";
		}
		else
		{
			cout<<"Eroare: valoarea lui dec_z este " << dec_z <<endl;
		}
	}
	
	// cout << endl<< dec_z << " (" << clock_diff(t1,t2)<< "ms)" << endl << std::flush;
}

/*************************************************************************************/
vector<CtxtBit> select (CtxtBit& c, vector<CtxtBit>& a, vector<CtxtBit>& b)
{
    vector<CtxtBit> ret;
	
	vector<CtxtBit> vt1, vt2;		
	for (int i = 0; i < a.size(); i++)
	{
		vt1.push_back (CtxtBit(c));
		vt2.push_back (CtxtBit(c));		
	}
		
	for  (int i = 0; i < a.size(); i++)
	{
        vt1[i] *= a[i];                      // vt1[i] *= a[i];
		vt2[i] += (*ctxt_of_1);           // vt2[i] += ENC(1);

		vt2[i] *= b[i];              // vt2[i] *= b[i];
        vt1[i] += vt2[i];                       // vt1[i] += vt2[i];
				
		ret.push_back(vt1[i]);		
	}

// 	if (ret[0].findBaseLevel() < 2) 
// 	{
// 		cout << endl << "select: recryption need =>findBaseLevel: " << ret[0].findBaseLevel() << std::flush;
// 		batchRecrypt(ret);
// 	}
	
	return ret;
}

/*************************************************************************************/
vector<CtxtBit> getmax (vector<vector<CtxtBit> >& vvct)
{	
	/*HE_Integer heInt(t_bits);	
	vector<CtxtBit>  ct_max = vvct[0];	
	
 	for (int i = 1; i < vvct.size(); i++)
 	{
		cout << endl << endl << "Step# " << i << "...";		
 		cout << endl << "Evaluation of " << heInt.decryptIntValue(ct_max) << " > " << heInt.decryptIntValue(vvct[i], *secretKey) << std::flush;
		cout << endl << "getmax (before enter to round) =>findBaseLevel: " << ct_max[0].findBaseLevel() << std::flush;		
		
		t1 = clock();
		if (ct_max[0].findBaseLevel() < 5)
		{
			cout << endl << "before compute_t: recryption need =>findBaseLevel: " << ct_max[0].findBaseLevel() << std::flush;
			batchRecrypt(ct_max);
			cout << endl << "before compute_t: recryption done =>findBaseLevel: " << ct_max[0].findBaseLevel() << std::flush;			
		}
		
 		CtxtBit ct_t = compute_t (0, t_bits, ct_max, vvct[i]);
		t2 = clock();
		cout << endl << "compute_t: " << heInt.decryptBit(ct_t, *secretKey) << " (" << clock_diff(t1,t2) << " ms)" 
		<< "; (baseLevel of ct_t:" << ct_t.findBaseLevel() << ")" << std::flush;
		
		if (ct_t.findBaseLevel() < 2) 
		{
			cout << endl << "after compute_t: recryption need =>findBaseLevel: " << ct_t.findBaseLevel() << std::flush;
			c txtRecrypt(ct_t);
			cout << endl << "after compute_t: recryption done =>findBaseLevel: " << ct_t.findBaseLevel() << std::flush;
		}
		
		t1 = clock();
 		ct_max = select (ct_t, ct_max, vvct[i]);
 		t2 = clock();
 		cout << endl << "select_ct_max: " << heInt.decryptIntValue(ct_max, *secretKey) << " (" << clock_diff(t1,t2) << " ms)" 
		<< "; (baseLevel of ct_max:" << ct_max[0].findBaseLevel() << ")" << std::flush;
 	}

	return ct_max; */
}

/*************************************************************************************
void batchRecrypt (vector<CtxtBit>& vct)
{
	if (vct.size() == 0) return;
	
	FHEPubKey& pk = (FHEPubKey &)vct[0].getPubKey();
	if (!pk.isBootstrappable()) 
	  return;
  		
	for (int i=0; i < vct.size(); i++)
	{
		cout << endl << "recrypt ctxt.......1" << std::flush;
		clock_t t1 = clock();
		pk.reCrypt(vct[i]);
		clock_t t2 = clock();
		cout << endl << "recrypt ctxt.......2" << " (" << clock_diff(t1,t2) << "ms)" << std::flush;
	}
}

*************************************************************************************
void ctxtRecrypt (CtxtBit& ct)
{
	FHEPubKey& pk = (FHEPubKey&) ct.getPubKey();
	if (!pk.isBootstrappable()) 
		return;
  
	pk.reCrypt(ct);  	
}

*************************************************************************************/

vector<CtxtBit> getmax (vector<vector<CtxtBit> >& vvct, int start, int n)
{	
 	/*assert (n >= 1);
	
	HE_Integer heInt(t_bits);	
 	if (n == 1) return vvct[start];
 		
	vector<CtxtBit> ct_max_1 = getmax (vvct, start, n/2);
 	vector<CtxtBit> ct_max_2 = getmax (vvct, start+n/2, n%2 == 0? n/2: n/2 +1);

	cout << "\n\nEvaluation of " << heInt.decryptIntValue(ct_max_1, *secretKey) << " > " << heInt.decryptIntValue(ct_max_2, *secretKey)
		 << "; baseLevel: " << ct_max_1[0].findBaseLevel() << ", " << ct_max_2[0].findBaseLevel() << std::flush;
	*/
	
	// FHE_NTIMER_START(getmax);	
// 	if (ct_max_1[0].findBaseLevel() < 6)
// 	{
// 		cout << endl << "recryption need =>findBaseLevel(ct_max_1): " << ct_max_1[0].findBaseLevel() << std::flush;
// 		batchRecrypt(ct_max_1);
// 		cout << endl << "recryption done =>findBaseLevel (ct_max_1): " << ct_max_1[0].findBaseLevel() << std::flush;			
// 	}
// 	if (ct_max_2[0].findBaseLevel() < 6)
// 	{
// 		cout << endl << "recryption need =>findBaseLevel(ct_max_2): " << ct_max_2[0].findBaseLevel() << std::flush;
// 		batchRecrypt(ct_max_2);
// 		cout << endl << "recryption done =>findBaseLevel (ct_max_2): " << ct_max_2[0].findBaseLevel() << std::flush;			
// 	}
		
 	/*FHE_NTIMER_START(compute_t);
	t1 = clock();
	CtxtBit ct_t = compute_t (0, t_bits, ct_max_1, ct_max_2);
	t2 = clock();
	FHE_NTIMER_STOP(compute_t);
	cout << endl << "ct_t: " << heInt.decryptBit(ct_t, *secretKey)<<" (" << clock_diff(t1,t2)	<< " sec); (baseLevel ct_t:" << ct_t.findBaseLevel() << ")" << std::flush;	
	
		
	FHE_NTIMER_START(select);
	t1 = clock();
	vector<CtxtBit> ct_max = select (ct_t, ct_max_1, ct_max_2);
	t2 = clock();
	FHE_NTIMER_STOP(select);
	FHE_NTIMER_STOP(getmax);
	cout << endl << "ct_max: " << heInt.decryptIntValue(ct_max, *secretKey) << " (" << clock_diff(t1,t2) << " sec); (baseLevel ct_max: " << ct_max[0].findBaseLevel() <<  ")" << std::flush;
		
	return ct_max;*/
	
//	cout << endl << "ct_t: " << heInt.decryptBit(ct_t, *secretKey) << "; (baseLevel ct_t:" << ct_t.findBaseLevel() << ")" << std::flush;		
//	cout << "\nct_max: " << heInt.decryptIntValue(ct_max, *secretKey) << "; (baseLevel ct_max: " << ct_max[0].findBaseLevel() <<  ")" << std::flush;
//	cout << endl << "ct_t: " << heInt.decryptBit(ct_t, *secretKey) << " (" << tm << " sec); (baseLevel ct_t:" << ct_t.findBaseLevel() << ")" << std::flush;	
//	cout << endl << "ct_max: " << heInt.decryptIntValue(ct_max, *secretKey) << " (" << tm << " sec); (baseLevel ct_max: " << ct_max[0].findBaseLevel() <<  ")" << std::flush;
}

 int gmax (vector<int>& v, int start, int n)
 {
 	assert (n >= 1);	
 	if (n == 1) return v[start];
 	
 	int max_1 = gmax (v, start, n/2);
 	int max_2 = gmax (v, start+n/2, n%2 == 0? n/2: n/2 +1);
 	
 	return max_1 > max_2? max_1: max_2;	
 }
