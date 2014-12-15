/**
    @file    blazar.hpp
    Main blazar++'s header file

    @author Mateusz Janiak
*/

#ifndef _blazar_hpp
#define _blazar_hpp 1

#include <scfgp.hpp>
#include <bazinga.hpp>
#include "electrons.hpp"

/**
   Read main config switches from configuration file
   @param scfgp 
*/

void getConfigSwitches( scfgp* _cfg ) {
  _cfg -> add<int>("lumSyn",0);
  _cfg -> add<int>("lumSsc",0);
  _cfg -> add<int>("lumExt1",0);
  _cfg -> add<int>("lumExt2",0);
  _cfg -> add<int>("lumBlrPl",0);
  _cfg -> add<int>("lumDtPl",0);
  _cfg -> add<int>("lumAccd",0);
  _cfg -> add<int>("lumBlrSp",0);
  _cfg -> add<int>("lumDtSp",0);
  _cfg -> add<int>("lumBlrGu",0);
  _cfg -> add<int>("lumDtGu",0);
  _cfg -> add<int>("lumAccdGu",0);
  _cfg -> add<int>("lumSynPoint",0);
  _cfg -> add<int>("lumSscPoint",0);
  _cfg -> add<int>("lumExt1Point",0);
  _cfg -> add<int>("lumExt2Point",0);
  _cfg -> add<int>("lumBlrPlPoint",0);
  _cfg -> add<int>("lumDtPlPoint",0);
  _cfg -> add<int>("lumAccdPoint",0);
  _cfg -> add<int>("lumBlrSpPoint",0);
  _cfg -> add<int>("lumDtSpPoint",0);
  _cfg -> add<int>("lumBlrGuPoint",0);
  _cfg -> add<int>("lumDtGuPoint",0);
  _cfg -> add<int>("lumAccdGuPoint",0);
  _cfg -> add<int>("Adiabatic",0);
  _cfg -> add<int>("Syn",0);
  _cfg -> add<int>("Ssc",0);
  _cfg -> add<int>("Ext1",0);
  _cfg -> add<int>("Ext2",0);
  _cfg -> add<int>("BlrPl",0);
  _cfg -> add<int>("DtPl",0);
  _cfg -> add<int>("BlrSp",0);
  _cfg -> add<int>("DtSp",0);
  _cfg -> add<int>("Accd",0);
  _cfg -> add<int>("BlrGu",0);
  _cfg -> add<int>("DtGu",0);
  _cfg -> add<int>("AccdGu",0);
  _cfg -> add<int>("QAccd",0);
  _cfg -> add<int>("SABS",0);
  _cfg -> add<double>("AdiabaticABG",0.666667);
  _cfg -> add<std::string>("output","blazar_out");
  _cfg -> add<int>("saveElectrons",1);
  _cfg -> add<int>("saveElectronsAvg",1);
  _cfg -> add<int>("saveLum",1);
  _cfg -> add<int>("saveLumPoint",1);
  _cfg -> add<int>("root",1);
  _cfg -> add<int>("printRoot",1);
  
  _cfg -> updateAddedParameters( ); }

/**
   After reading parameters apply some necessary changes to them
   @param scfgp
   @param electrons
*/

void safetyConfigSwitches( scfgp* _cfg, electrons* _ele ) {
  /* some safety features */
  if( _ele -> ifEvol( ) &&  _cfg->get<int>("Syn") ) { _cfg->modify<int>("lumSyn",1); }
  if( _ele -> ifEvol( ) &&  _cfg->get<int>("Ssc") ) { _cfg->modify<int>("lumSsc",1); }

  if( _cfg->get<int>("lumSynPoint") ) { _cfg->modify<int>("lumSyn",1); }
  if( _cfg->get<int>("lumSscPoint") ) { _cfg->modify<int>("lumSsc",1); }
  if( _cfg->get<int>("lumExt1Point") ) { _cfg->modify<int>("lumExt1",1); }
  if( _cfg->get<int>("lumExt2Point") ) { _cfg->modify<int>("lumExt2",1); }
  if( _cfg->get<int>("lumBlrPlPoint") ) { _cfg->modify<int>("lumBlrPl",1); }
  if( _cfg->get<int>("lumDtPlPoint") ) { _cfg->modify<int>("lumDtPl",1); }
  if( _cfg->get<int>("lumAccdPoint") ) { _cfg->modify<int>("lumAccd",1); }
  if( _cfg->get<int>("lumBlrSpPoint") ) { _cfg->modify<int>("lumBlrSp",1); }
  if( _cfg->get<int>("lumDtSpPoint") ) { _cfg->modify<int>("lumDtSp",1); }
  if( _cfg->get<int>("lumBlrGuPoint") ) { _cfg->modify<int>("lumBlrGu",1); }
  if( _cfg->get<int>("lumDtGuPoint") ) { _cfg->modify<int>("lumDtGu",1); }
  if( _cfg->get<int>("lumAccdGuPoint") ) { _cfg->modify<int>("lumAccdGu",1); }

  if( _cfg->get<int>("lumSsc") ) { _cfg->modify<int>("Syn",1); }
  if( _cfg->get<int>("lumSsc") ) { _cfg->modify<int>("lumSyn",1); }

  if( _cfg->get<int>("lumSyn") ) { _cfg->modify<int>("Syn",1); }
  if( _cfg->get<int>("lumSsc") ) { _cfg->modify<int>("Ssc",1); }
  if( _cfg->get<int>("lumExt1") ) { _cfg->modify<int>("Ext1",1); }
  if( _cfg->get<int>("lumExt2") ) { _cfg->modify<int>("Ext2",1); }
  if( _cfg->get<int>("lumBlrPl") ) { _cfg->modify<int>("BlrPl",1); }
  if( _cfg->get<int>("lumDtPl") ) { _cfg->modify<int>("DtPl",1); }
  if( _cfg->get<int>("lumAccd") ) { _cfg->modify<int>("Accd",1); }
  if( _cfg->get<int>("lumBlrSp") ) { _cfg->modify<int>("BlrSp",1); }
  if( _cfg->get<int>("lumDtSp") ) { _cfg->modify<int>("DtSp",1); }
  if( _cfg->get<int>("lumBlrGu") ) { _cfg->modify<int>("BlrGu",1); }
  if( _cfg->get<int>("lumDtGu") ) { _cfg->modify<int>("DtGu",1); }
  if( _cfg->get<int>("lumAccdGu") ) { _cfg->modify<int>("AccdGu",1); }

  if( _cfg->get<int>("lumSscPoint") ) { _cfg->modify<int>("lumSynPoint",1); }
}

#endif /* _blazar_hpp */
