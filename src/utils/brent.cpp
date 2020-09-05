#include <stdio.h>
#include <algorithm>
#include <cmath>
#include "brent.h"

// This will work generally for any class where the constraint function being solved is set up to take the variable it is solving for as its first argument, and as its second argument, a struct containing any additional arugments that might change between instances of the solver, such as a zone index.
// That member function is passed here as func

#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember)) // This would be better done with std::invoke, but that requires c++17

// implementation of Brent solver, by Nathaniel Roth and Paul Duffell

template <class cl,class args>
double brent_solver<cl,args>::solve(cl& cl_instance, constraintMemFn func, args* passed_args, double aa, double bb, double eps, int max_iters, int* n) 
{
  double a = aa;
  double b = bb;

  double fa = CALL_MEMBER_FN(cl_instance,func)(a,passed_args);
  double fb = CALL_MEMBER_FN(cl_instance,func)(b,passed_args);
  *n = 0;
  if( fa*fb >= 0 ){
    printf("Brent failed; not bracketed properly\n");
    *n = -1;
    return(0.0);
  }
  if( fabs(fa) < fabs(fb) ){
    double temp = a;
    a = b;
    b = temp;
    temp = fa;
    fa = fb;
    fb = temp;
  }
     double c = a;
   double fc = fa;
   int stillrunning = 1;
   bool bisectflag = true;
   double bdelta = 5.e-16 * std::min(fabs(a),fabs(b)); // maybe could be further optimized?
   double s;
   double d;
   while( stillrunning ){
     //      if (*n > 0)  printf("%e %e %e %e %e %e %d\n",b,a,c,func(b),func(a),func(c),*n); // debugging output
      if( fa != fc && fb != fc ){
	s = a*fb*fc/(fa-fb)/(fa-fc) + b*fc*fa/(fb-fc)/(fb-fa) + c*fa*fb/(fc-fa)/(fc-fb); // inverse quadratic interpolation
      }else{
	s = b - fb*(b-a)/(fb-fa); // secant method
      }

      bool cond1 = (s < (3. * a + b)/4. && s < b) || (s > (3. * a + b)/4. && s > b);
      bool cond2 = bisectflag && (fabs(s - b) >= fabs(b - c)/2.) ;
      bool cond3 = !bisectflag && (fabs(s - b) >= fabs(c - d)/2.) ;
      bool cond4 = bisectflag && (fabs(b - c) < bdelta) ;
      bool cond5 = !bisectflag && (fabs(c - d) < bdelta);

      if (cond1 || cond2 || cond3 || cond4 || cond5)
	{
	  s = (a + b)/2.;
	  bisectflag = true;
	}
      else
	{
	  bisectflag = false;
	}

      double fs = CALL_MEMBER_FN(cl_instance,func)(s,passed_args);
      d = c;
      c = b;
      fc = fb;

      if( fa*fs < 0 ){ 
         b = s;
	 fb = fs;
      }else{
         a = s;
         fa = fs; 
      }

      if( fabs(fa) < fabs(fb) ){
         double temp = a; 
         a = b; 
         b = temp;
         temp = fa;
         fa = fb;
         fb = temp;      
      }
      
      *n = *n+1;
      if( fb == 0. || fs == 0. || fabs(a-b) < eps * std::min(fabs(a),fabs(b)) ) stillrunning = 0;

   }
   return( b );

}

