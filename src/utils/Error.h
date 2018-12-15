// 
// File    : Error.hh
// ------------------
// Created : Tue May 27 10:52:24 2008
// Authors : Rollin C. Thomas (RCT) - rcthomas@lbl.gov
// Purpose : Tired of not having errors show up right...
//  
// $Header: /cvs/mp22/ddmc_sn/Error.hh,v 1.1 2008/05/27 18:06:30 rthomas Exp $ 
// 

#ifndef __ERROR__
#define __ERROR__

#include <stdexcept>

class Error : public std::runtime_error
{

   public :

      Error( const std::string& message ) : std::runtime_error( message ) {}

};

#endif
