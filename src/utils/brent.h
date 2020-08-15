#ifndef _BRENT_H
#define _BRENT_H

// This will work generally for any class that has a member function with a constraint function being solved, provided that the constrant function takes the variable it is solving for as its first argument, and as its second argument, a struct containing any additional arugments that might change between instances of the solver, such as a zone index.
// That member function is passed here through the solve() method, and typecast as constraintMemFn

template <class cl, class args>
class brent_solver
{

public:

  // Right now only set up for solving equations that work on double precision
  // Could possibly be modified to use SedonaReal?
  typedef double (cl::*constraintMemFn)(double, args*);

  // This is the interface to the outside world
  // Takes the class containing the constraint equation, the member function of that class which expresses the constraint equation, a struct with any additional arguments, the lower bracket for the root, the upper bracket for the root, the desired relative precision for the solution, the maximum number of allowed iterations, and an int that will store the number of iterations required by the solver
  double solve(cl&, constraintMemFn, args*, double, double, double, int, int*);
};

#endif
