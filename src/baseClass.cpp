/**
    @file baseClass.cpp
    @author Mateusz Janiak
*/

#include "baseClass.hpp"

baseClass::baseClass( scfgp* _cfg, jetGeometry* _r, std::string _id ) : id(_id), cfg(_cfg), r(_r) { }

baseClass::~baseClass() {
  cfg = NULL;
  r = NULL; }

std::string baseClass::whoAmI() { return id; }

double baseClass::beta( double x ){ return sqrt(1.0-1.0/(x*x)); }
