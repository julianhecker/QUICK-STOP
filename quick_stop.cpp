#include <iostream>
#include <stdlib.h>
using namespace std;

#include <malloc.h>
#include <stdio.h>


#include <fstream>

#include <cstring>

#include <errno.h>

#include <limits.h>

#include <float.h>

#include <math.h>


#include <stdio.h>
#include <time.h>

#ifndef MERSENNETWISTER_H
#define MERSENNETWISTER_H

/*
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.
   Before using, initialize the state by using init_genrand(seed)
   or init_by_array(init_key, key_length).
   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:
   1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
   2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
   3. The names of its contributors may not be used to endorse or promote
   products derived from this software without specific prior written
   permission.
   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
   */

struct MersenneTwister
{
  /* Period parameters */
  const static int N = 624;
  const static int M = 397;
  const static unsigned long long MATRIX_A = 0x9908b0dfUL;  /* constant vector a */
  const static unsigned long long UPPER_MASK = 0x80000000UL;  /* most significant w-r bits */
  const static unsigned long long LOWER_MASK = 0x7fffffffUL;  /* least significant r bits */

  unsigned long mt[N]; /* the array for the state vector  */
  int mti; /* mti==N+1 means mt[N] is not initialized */

  MersenneTwister()
  {
    mti = N + 1;
    unsigned long init[4] = { 0x123, 0x234, 0x345, 0x456 }, length = 4;
    init_by_array( init, length );
  }

  MersenneTwister( unsigned int seed )
  {
    mti = N + 1;
    init_genrand( seed );
  }

private:
  /* initializes mt[N] with a seed */
  void init_genrand( unsigned long s )
  {
    mt[0] = s & 0xffffffffUL;
    for( mti = 1; mti < N; mti++ ) {
      mt[mti] =
        (1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
      /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
      /* In the previous versions, MSBs of the seed affect   */
      /* only MSBs of the array mt[].                        */
      /* 2002/01/09 modified by Makoto Matsumoto             */
      mt[mti] &= 0xffffffffUL;
      /* for >32 bit machines */
    }
  }

  /* initialize by an array with array-length */
  /* init_key is the array for initializing keys */
  /* key_length is its length */
  /* slight change for C++, 2004/2/26 */
  void init_by_array( unsigned long init_key[], int key_length )
  {
    int i, j, k;
    init_genrand( 19650218UL );
    i = 1; j = 0;
    k = (N > key_length ? N : key_length);
    for( ; k; k-- ) {
      mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1664525UL))
        + init_key[j] + j; /* non linear */
      mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
      i++; j++;
      if( i >= N ) { mt[0] = mt[N - 1]; i = 1; }
      if( j >= key_length ) j = 0;
    }
    for( k = N - 1; k; k-- ) {
      mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1566083941UL))
        - i; /* non linear */
      mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
      i++;
      if( i >= N ) { mt[0] = mt[N - 1]; i = 1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
  }

  /* generates a random number on [0,0xffffffff]-interval */
  unsigned long randUInt()
  {
    unsigned long y;
    static unsigned long mag01[2] = { 0x0UL, MATRIX_A };
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if( mti >= N ) { /* generate N words at one time */
      int kk;

      if( mti == N + 1 )   /* if init_genrand() has not been called, */
        init_genrand( 5489UL ); /* a default initial seed is used */

      for( kk = 0; kk < N - M; kk++ ) {//227 iterations
        y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
        mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1UL];
      }
      for( ; kk < N - 1; kk++ ) {//goes the 397 remaining iterations (624 iterations total)
        y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
        mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
      }
      y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
      mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];

      mti = 0;
    }

    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
  }

public:
  /* generates a random number on [0,0x7fffffff]-interval */
  // (0 and 2147483647)
  inline int randInt() { return randUInt() >> 1; }
  inline float randFloat() { return randInt()*(1.f / 4294967295.f); }
  inline double randDouble() { return randInt()*(1.0 / 4294967295.0); }
};

#endif

inline double bin_log(unsigned long ctr,unsigned long n, double p)
{
   double tmp=0;
   if(fabs(p)>pow(10,-16))
   {
	  tmp+=(double)(ctr)*log(p);
   }
   if(fabs(1.0-p)>pow(10,-16))
   {
	  tmp+=((double)(n)-(double)(ctr))*log(1.0-p);
   }
   return tmp;

}


int main(int argc, char * argv[])
{

        if(argc!=7){ cout << "incorrect number of parameters."<<endl; exit(1);}
		double p_true =atof(argv[1]);
		double p_cut1  =atof(argv[2]);
		double p_cut2  =atof(argv[3]);
		double alpha  =atof(argv[4]);
		int seed=atoi(argv[5]);
		unsigned long repsp=atoi(argv[6]);
		
		
		unsigned long n;
		unsigned long ctr;
		unsigned long reps=repsp;
        unsigned int i;
		
		
	    double lw_bd;
		double lw_bd1,lw_bd2;
		double av_ps;
		
		unsigned int stop;
		unsigned int tmpx;
		double mle;
		double mle_adaptive;
		double d1,d2,d3;
		double stat1,stat2;
		double adaptive_stat;
		double logalpha=log(1/alpha);
		MersenneTwister mt(seed);
		
	    double overallperms=0.0;
		
		for(i=1;i<=reps;i++)
		{
				ctr=0;
				n=0;
				stop=0;
				
				stat1=0;
				stat2=0;
				while(stop==0)
				{
				   tmpx=0;
				   if(2.0*mt.randDouble()<=p_true)
				   {
					 tmpx=1;
				   }
				   n++;
				   if(n==1)
				   {
					 adaptive_stat=log(0.5);
				   }
				   if(n>1)
				   {
					  mle_adaptive=((double)(ctr)+0.5)/((double)(n));
					 
					  adaptive_stat+=bin_log(tmpx,(unsigned long)1,mle_adaptive);
				   }
				   ctr+=tmpx;
				   mle=((double)(ctr))/((double)n);
				   
				   d1=bin_log(ctr,n,mle);
				   d2=bin_log(ctr,n,p_cut1);
				   d3=bin_log(ctr,n,p_cut2);
				   if(mle<p_cut1)
				   {
					  stat1=adaptive_stat-d1;
				   }
				   else
				   {
				      stat1=adaptive_stat-d2;
				   }
				   if(mle>p_cut2)
				   {
					  stat2=adaptive_stat-d1;
				   }
				   else
				   {
				      stat2=adaptive_stat-d3;
				   }
				   if(stat1>=logalpha)
				   {
					 stop=2;
				   }
				   if(stat2>=logalpha)
				   {
					 stop=1;
				   }
				   
				  
				}
				overallperms+=(double)n;
				
		}
		cout << "QUICK-STOP: "<<endl;
		cout << "p: " << p_true<<endl;
		cout << "p_1: " << p_cut1<<endl;
		cout << "p_2: " << p_cut2<<endl;
		cout << "alpha: " << alpha<<endl;
		cout << "average number of permutations/simulations: "<< (double)overallperms/(double)reps<<endl;
		av_ps=(double)overallperms/(double)reps;
		if(p_true>p_cut2 && p_true!=0.0 && p_true!=1.0)
		{ 
	       cout << "lower bound: " << log(1/alpha)/(p_true*log(p_true/p_cut1)+(1-p_true)*log((1-p_true)/(1-p_cut1)))<<endl;
		   lw_bd=log(1/alpha)/(p_true*log(p_true/p_cut1)+(1-p_true)*log((1-p_true)/(1-p_cut1)));
		}
		if(p_true<p_cut1 && p_true!=0.0 && p_true!=1.0)
		{ 
	       cout << "lower bound: " << log(1/alpha)/(p_true*log(p_true/p_cut2)+(1-p_true)*log((1-p_true)/(1-p_cut2)))<<endl;
		   lw_bd=log(1/alpha)/(p_true*log(p_true/p_cut2)+(1-p_true)*log((1-p_true)/(1-p_cut2)));
		}
		if(p_true>p_cut1 && p_true<p_cut2)
		{ 
	       lw_bd1=log(1/alpha)/(p_true*log(p_true/p_cut1)+(1-p_true)*log((1-p_true)/(1-p_cut1)));
		   lw_bd2=log(1/alpha)/(p_true*log(p_true/p_cut2)+(1-p_true)*log((1-p_true)/(1-p_cut2)));
		   lw_bd=lw_bd2;
	       if(lw_bd1<=lw_bd2)
		   {
			   lw_bd=lw_bd1;
			   
		   }
		   cout << "lower bound: " << lw_bd <<endl;
		}
		if(p_true==0.0)
		{
			lw_bd=log(1/alpha)/(log((1-p_true)/(1-p_cut2)));
			cout << "lower bound: " << lw_bd <<endl;
		}
		if(p_true==1.0)
		{
			lw_bd=log(1/alpha)/(log(p_true/p_cut1));
			cout << "lower bound: " << lw_bd <<endl;
		}
		cout << "ratio QUICK-STOP/lower bound: " << av_ps/lw_bd<<endl;	
		
	
        
}