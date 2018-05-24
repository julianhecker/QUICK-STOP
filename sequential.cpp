unsigned long ctr=0; // counter # successes
unsigned long n=0; // counter permutations

unsigned long tmpx;	// temporary variable, specifies if current permutation resulted in more extreme
// 					   test statistic or not.			
double mle,d1,d2,d3,mle_adaptive; // temporary variables


double tau_1=0; // tau_1 as described in reference
double tau_2=0; // tau_2 as described in reference
double p_cut1,p_cut2; // p_1 and p_2 as defined in reference => d is determined by p_cut2-p_cut1
double alpha1; //alpha1 as defined in reference
double alpha2; //alpha2 as defined in reference
double logalpha1=log(1/alpha1); 
double logalpha2=log(1/alpha2); 
int delta=0; // decision variable as defined in reference

while(delta==0)
{
	               #############################################
				   tmpx=new_permuted_stat_greater_or_equal(...);  // arbitrary function that performs permutation and compares the 
				   // 												 permuted test statistic with the observed test statistic. Returns 1
				   //												 if the permuted test statistic is greater or equal and 0 otherwise
				   #############################################
				   n++;
				   if(n==1){ adaptive_stat=log(0.5); } // initial value for p_1 in reference
				   if(n>1)
				   {  
					 mle_adaptive=((double)(ctr)+0.5)/((double)(n)); // update of p_i, see reference
					 adaptive_stat+=bin_log(tmpx,(unsigned long)1,mle_adaptive); 
				   }
				   ctr+=tmpx;
				   mle=((double)(ctr))/((double)n); // current MLE estimate := #successes/#permutations
				   
				   // Update of objects tau_1 and tau_2
				   d1=bin_log(ctr,n,mle);
				   d2=bin_log(ctr,n,p_cut1);
				   d3=bin_log(ctr,n,p_cut2);
				   if(mle<p_cut1){ tau_1=adaptive_stat-d1; }
				   else{  tau_1=adaptive_stat-d2;  }
				   if(mle>p_cut2) {  tau_2=adaptive_stat-d1; }
				   else{ tau_2=adaptive_stat-d3; }
				   // make decision if threshold is reached
				   if(tau_1>=logalpha1){ delta=2; }
				   if(tau_2>=logalpha2){ delta=1; }
}

//....// help function to evaluate binomial distribution
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