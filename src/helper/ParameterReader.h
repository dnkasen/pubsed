#ifndef _PARAMETER_READER_H
#define _PARAMETER_READER_H 1

#include "Lua.h"
#include <iostream>
#include <tuple>

//-------------------------------------------------------
// parameter_reader
// This is a simple wrapper class that uses the
// Lua interface to read parameter file.
// If it can't find an entry in the parameter file
// it looks in the default file 
//-------------------------------------------------------
class ParameterReader
{
 private:
  
  Lua params_;
  Lua defaults_;
  int verbose_;

 public:

  // simple constructor
  ParameterReader()
  { 
    verbose_ = 0;
  }

  // constructor with init
  ParameterReader(std::string pfile,int v)
  {
    initialize(pfile,v);
  }

  
  // initialize the parameter file and defaults file
  void initialize(std::string pfile, int v)
  {
    verbose_ = v;

    // open parameter file
    params_.init(pfile);

    // open defaults file
    std::pair<std::string,bool> dfile;
    dfile = params_.scalar_pair<std::string>("defaults_file");
    if (!dfile.second) 
    {
      std::cout << "Must list a defaults_file in input param file!, exiting\n";
      exit(1); 
    }
    defaults_.init(dfile.first);
  }


  // get a scalar quantity
  template< typename T > T getScalar( const char* label )
  {
    // look for  value in the parameter file
    std::pair< T, bool > value = params_.scalar_pair< T >( label );
    // if there, return it
    if (value.second) return value.first;
    // otherwise look in the defaults file
    value = defaults_.scalar_pair< T >( label );
    // if there, return it
    if (value.second) return value.first;
    // otherwise we have a problem
    if (verbose_) std::cout << "Can't find parameter \"" << label << 
		   "\" in either param file or defaults; exiting\n";
    exit(1);
  }


  // get a vector quantity
  template< typename T > std::vector< T > getVector( const char* label )
  {
    // look for  value in the parameter file
    std::pair< std::vector<T>, bool > value = params_.vector_pair< T >( label );
    //if there, return it
    if (value.second) return value.first;
    // otherwise look in the defaults file
    value = defaults_.vector_pair< T >( label );
    // if there, return it
    if (value.second) return value.first;
    // otherwise we have a problem
    if (verbose_) std::cout << "Can't find parameter \"" << label << 
		   "\" in either param file or defaults; exiting\n";
    exit(1);
  }

  template< typename T > std::vector< T > getArray( const char* label )
  {   
    std::vector<T> list = getVector<double>(label);
    std::vector<T> arr;

    double start,stop,del;
    if (list.size() == 1)
    {
      start = list[0];
      stop  = list[0];
      del   = 1;
    }
    else if (list.size() == 2)
    {
      std::cout << "List " << label << " cannot have 2 numbers\n";
      exit(1);
    }
    else 
    {
      start = list[0];
      stop  = list[1];
      del   = list[2];
    }
  
    for (double v = start; v <= stop;)
    {
      arr.push_back(v);
      if (list.size() == 4) v = v*(1.0+del);
      else v += del;
    }
    return arr;
  }
};


#endif
