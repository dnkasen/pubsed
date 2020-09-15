#ifndef _BRENT_H
#define _BRENT_H

// This will work generally for any class that has a member function with a constraint function being solved, provided that the constrant function takes the variable it is solving for as its first argument, and as its second argument, a struct containing any additional arugments that might change between instances of the solver, such as a zone index.
// That member function is passed here through the solve() method, and typecast as constraintMemFn

#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember)) // This would be better done with std::invoke, but that requires c++17

template <class cl, class args>
class brent_solver
{

public:

  // Right now only set up for solving equations that work on double precision
  // Could possibly be modified to use SedonaReal?
  typedef double (cl::*constraintMemFn)(double, args*);


  double solve(cl& cl_instance, constraintMemFn func, args* passed_args, double aa, double bb, double eps, int max_iters, int* n) 
{

  // This is the interface to the outside world
  // Takes the class containing the constraint equation, the member function of that class which expresses the constraint equation, a struct with any additional arguments, the lower bracket for the root, the upper bracket for the root, the desired relative precision for the solution, the maximum number of allowed iterations, and an int that will store the number of iterations required by the solver

  // The class member function that expresses the constraint being solved is passed here as func
  
  // implementation of Brent solver, by Nathaniel Roth and Paul Duffell, largely based on the wikipedia description as of September 5, 2020
  
  double a = aa;
  double b = bb;

  double fa = CALL_MEMBER_FN(cl_instance,func)(a,passed_args);
  double fb = CALL_MEMBER_FN(cl_instance,func)(b,passed_args);
  *n = 0;
  if( fa*fb >= 0 ){
    //printf("Brent failed; not bracketed properly\n");
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
   double bdelta = 1.e-15 * std::min(fabs(a),fabs(b)); // maybe could be further optimized?
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
      if( *n > max_iters)
      {
	printf("Brent failed; %d iterations out of %d\n",*n,max_iters);
	return(b);	
      }
      if( fb == 0. || fs == 0. || fabs(a-b) < eps * std::min(fabs(a),fabs(b)) ) stillrunning = 0;

   }
   return( b );

}

  
};

#endif
