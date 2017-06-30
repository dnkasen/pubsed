// 
// File    : Lua.hh
// ----------------
// Created : Mon Mar 17 14:51:58 2008
// Authors : Rollin C. Thomas (RCT) - rcthomas@lbl.gov
// Purpose : Lua interface.
//  
// $Header: /cvs/mp22/ddmc_sn/Lua.cc,v 1.4 2008/05/27 18:06:54 rthomas Exp $
// 

#include "Error.h"
#include "Lua.h"

Lua::Lua()
{
   _status = false;
}

Lua::~Lua() 
{ 
   if( _status ) close(); 
}

Lua::Lua( std::string script )
{ 
   _has_add_path = false; 
   init( script ); 
}

Lua::Lua( std::string script, std::string path ) 
{ 
   _has_add_path = true; 
   _add_path = path; 
   init( script ); 
}

lua_State* Lua::state() const 
{ 
   return _lua_state; 
}

void Lua::init( std::string script )
{
   _lua_state = luaL_newstate();
   luaL_openlibs( _lua_state );

   if( _has_add_path )
   {
      lua_getglobal( _lua_state, "package" );
      lua_pushstring( _lua_state, "path" );
      lua_gettable( _lua_state, -2 );
      std::string path = (std::string) lua_tostring( _lua_state, -1 );
      lua_pop( _lua_state, 1 );
      path += ";" + _add_path;
      
      lua_pushstring( _lua_state, "path" );
      lua_pushstring( _lua_state, path.c_str() );
      lua_settable( _lua_state, -3 );
      lua_pop( _lua_state, 1 );
   }

   if( luaL_loadfile( _lua_state, script.c_str() ) || lua_pcall( _lua_state, 0, 0, 0 ) )
      error( "Configuration file error: %s\n", lua_tostring( _lua_state, -1 ) );
   _status = true;
}

void Lua::close()
{
   lua_close( _lua_state );
   _status = false;
   _has_add_path = false;
}

void Lua::error( const char *format, ... )
{
   va_list argp;
   va_start( argp, format );
   vfprintf( stderr, format, argp );
   va_end( argp );
   lua_close( _lua_state );
   throw Error( "Lua error" );
}

void Lua::put_key( const char* key )
{
   lua_pushstring( _lua_state, key );
}

void Lua::get_value( const char* param, std::string& value )
{
   if( ! lua_isstring( _lua_state, -1 ) ) error( "'%s' should be a string\n", param );
   value = std::string( lua_tostring( _lua_state, -1 ) );
}

void Lua::get_value( const char* param, bool& value )
{
   if( ! lua_isboolean( _lua_state, -1 ) ) error( "'%s' should be a boolean\n", param );
   value = (bool)lua_toboolean( _lua_state, -1 );
}

std::pair< double, bool > Lua::get_function_pair(const char* param, double x)
{
   std::pair< double, bool > value;
   value.second = true;

   // set variable name
   lua_getglobal(_lua_state, param);


   // check if it exists
   if( lua_isnil( _lua_state, -1 ) ) 
   {
      value.second = false;
   }
   else
   {
         printf("HI\n");
      // check if this variable is a function
      if (lua_isfunction(_lua_state, -1))
      {
         // pass the first argument 
         lua_pushnumber(_lua_state, x);
         // call the function with 1 arguments, return 1 result 
         lua_call(_lua_state, 1, 1);
         // get the result
         value.first = (double)lua_tonumber(_lua_state, -1);
         printf("YEAH\n");
       }
      // if not function try to return scalar
      else
      {
         printf("GOG\n");
         get_value( param, value.first );
         value.second = true;
      }
    }
   lua_pop(_lua_state, 1);
   return value;
}
