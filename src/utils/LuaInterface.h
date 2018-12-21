//
// File    : Lua.hh
// ----------------
// Created : Mon Mar 17 14:51:58 2008
// Authors : Rollin C. Thomas (RCT) - rcthomas@lbl.gov
// Purpose : Lua interface.
//
// $Header: /cvs/mp22/ddmc_sn/Lua.hh,v 1.4 2008/05/09 00:46:07 rthomas Exp $
//

#ifndef __LUA__
#define __LUA__

#include <vector>
#include <string>

extern "C"
{
#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
}

/// @class Lua
/// @brief Interface to Lua via a wrapped Lua state.

class LuaInterface
{

   private :

      bool        _status;       ///< Status of the Lua state.
      lua_State*  _lua_state;    ///< Lua state itself.

      bool        _has_add_path; ///< Whether or not to add path information.
      std::string _add_path;     ///< Additional path information.

   public :

      /// Constructor.
      LuaInterface();

      /// Destructor.
      ~LuaInterface();

      /// Constructor, initializes state from a script.
      LuaInterface( std::string script );

      /// Constructor, initializes state from a script with additional path information.
      LuaInterface( std::string script, std::string path );

      /// Access to the lua state.
      lua_State* state() const;

      /// Initialize existing object from a script.
      void init( std::string script );

      /// Close the Lua state.
      void close();

      /// Report an error.
      void error( const char *format, ... );

      /// Obtain scalar data from the Lua state, forcing an error if it is nil.
      template< typename T > T scalar( const char* param );

      /// Obtain vector data from the Lua state, forcing an error if it is nil.
      template< typename T > std::vector< T > vector( const char* param );

      /// Obtain scalar data from the Lua state, without error even if nil.  Second value is true if not nil.
      template< typename T > std::pair< T, bool > scalar_pair( const char* param );

      /// Obtain vector data from the Lua state, without error even if nil.  Second value is true if not nil.
      template< typename T > std::pair< std::vector< T >, bool > vector_pair( const char* param );

      /// Obtain a scalar element from a table in the Lua state, forcing an error if nil.
      template< typename T, typename U > T scalar( const char* param, U const key );

      /// Obtain a vector element from a table in the Lua state, forcing an error if it is nil.
      template< typename T, typename U > std::vector< T > vector( const char* param, U const key );

      /// Obtain a scalar element from a table in the Lua state, without error even if nil.  Second value is true if not nil.
      template< typename T, typename U > std::pair< T, bool > scalar_pair( const char* param, U const key );

      /// Obtain a vector element from a table in the Lua state, without error even if nil.  Second value is true if not nil.
      template< typename T, typename U > std::pair< std::vector< T >, bool > vector_pair( const char* param, U const key );

      std::pair<double,bool> get_function_pair(const char* param, double x);


   private :

      /// Puts a key on the end of the Lua state.
      template< typename T > void put_key( T const key );
      void put_key( const char* key );

      /// Obtains a single value from the Lua state.
      template< typename T > void get_value( const char* param, T& value );
      void get_value( const char* param, std::string& value );
      void get_value( const char* param, bool& value );

};

template< typename T >
T LuaInterface::scalar( const char* param )
{
   std::pair< T, bool > value = scalar_pair< T >( param );
   if( ! value.second ) error( "'%s' must not be nil\n", param );
   return value.first;
}

template< typename T >
std::vector< T > LuaInterface::vector( const char* param )
{
   std::pair< std::vector< T >, bool > value = vector_pair< T >( param );
   if( ! value.second ) error( "'%s' must not be nil\n", param );
   return value.first;
}

template< typename T >
std::pair< T, bool > LuaInterface::scalar_pair( const char* param )
{
   std::pair< T, bool > value;
   lua_getglobal( _lua_state, param );
   if( lua_isnil( _lua_state, -1 ) )
   {
      value.second = false;
   }
   else
   {
      get_value( param, value.first );
      value.second = true;
   }
   lua_pop( _lua_state, -1 );
   return value;
}

template< typename T >
std::pair< std::vector< T >, bool > LuaInterface::vector_pair( const char* param )
{
   T buffer;
   std::pair< std::vector< T >, bool > value;
   lua_getglobal( _lua_state, param );
   if( lua_isnil( _lua_state, -1 ) )
   {
      value.second = false;
   }
   else
   {
      if( ! lua_istable( _lua_state, -1 ) ) error( "'%s' should be a table\n", param );
      lua_pushnil( _lua_state );
      while( lua_next( _lua_state, -2 ) != 0 )
      {
         get_value( param, buffer );
         (value.first).push_back( buffer );
         lua_pop( _lua_state, 1 );
      }
      value.second = true;
   }
   lua_pop( _lua_state, 1 );
   return value;
}

template< typename T, typename U >
T LuaInterface::scalar( const char* param, U const key )
{
   std::pair< T, bool > value = scalar_pair< T >( param, key );
   if( ! value.second ) error( "'%s' table lookup is nil\n", param );
   return value.first;
}

template< typename T, typename U >
std::vector< T > LuaInterface::vector( const char* param, U const key )
{
   std::pair< std::vector< T >, bool > value = vector_pair< T >( param, key );
   if( ! value.second ) error( "'%s' table lookup is nil\n", param );
   return value.first;
}

template< typename T, typename U >
std::pair< T, bool > LuaInterface::scalar_pair( const char* param, U const key )
{
   std::pair< T, bool > value;
   lua_getglobal( _lua_state, param );
   if( lua_isnil( _lua_state, -1 ) )
   {
      value.second = false;
   }
   else
   {
      if( ! lua_istable( _lua_state, -1 ) ) error( "'%s' should be a table\n", param );
      put_key( key );
      lua_gettable( _lua_state, -2 );
      get_value( param, value.first );
      lua_pop( _lua_state, 1 );
      value.second = true;
   }
   lua_pop( _lua_state, 1 );
   return value;
}

template< typename T, typename U >
std::pair< std::vector< T >, bool > LuaInterface::vector_pair( const char* param, U const key )
{
   T buffer;
   std::pair< std::vector< T >, bool > value;
   lua_getglobal( _lua_state, param );
   if( lua_isnil( _lua_state, -1 ) )
   {
      value.second = false;
   }
   else
   {
      if( ! lua_istable( _lua_state, -1 ) ) error( "'%s' should be a table\n", param );
      put_key( key );
      lua_gettable( _lua_state, -2 );
      if( lua_isnil( _lua_state, -1 ) )
      {
         value.second = false;
      }
      else
      {
         if( ! lua_istable( _lua_state, -1 ) ) error( "'%s' should be a table\n", param );
         lua_pushnil( _lua_state );
         while( lua_next( _lua_state, -2 ) != 0 )
         {
            get_value( param, buffer );
            (value.first).push_back( buffer );
            lua_pop( _lua_state, 1 );
         }
         value.second = true;
      }
      lua_pop( _lua_state, 1 );
   }
   lua_pop( _lua_state, 1 );
   return value;
}

template< typename T >
void LuaInterface::put_key( T const key )
{
   lua_pushnumber( _lua_state, key );
}

template< typename T >
void LuaInterface::get_value( const char* param, T& value )
{
   if( ! lua_isnumber( _lua_state, -1 ) ) error( "'%s' should be number\n", param );
   value = (T)lua_tonumber( _lua_state, -1 );
}





#endif
