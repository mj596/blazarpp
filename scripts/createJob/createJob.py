import libConfigFile

conf = libConfigFile.configFile('test_AccdSphericalGu')
conf.setExe('/home/mjaniak/Soft/blazarpp/blazarpp')
conf.addOption('#PBS -q medium')
conf.addOption('#PBS -l walltime=24:0:0')
conf.addOption('#PBS -l mem=3000mb')
conf.setDir('/home/mjaniak/Work/pracaJetScan/test_accdSpehrcalGu',True)
conf.addIntParameter('Adiabatic',1)
conf.addIntParameter('lumBlrGuPoint',1)
conf.addIntParameter('lumDtGuPoint',1)
conf.addIntParameter('lumAccdGuPoint',1)
conf.addIntParameter('Ninj',10)
conf.addIntParameter('thetaJ',0.1)
conf.addIntParameter('thetaObs',0.1)
conf.addDoubleParameter('nSave',-1)
conf.addIntParameter('evol',1)
conf.addStringParameter('eleModel','steady')
conf.addStringParameter('lumModel','steady')
conf.addStringParameter('magModel','steady')
conf.addDoubleParameter('p1',0.5)
conf.addDoubleParameter('p2',2.5)
conf.addDoubleParameter('gammaMin',1.0)
conf.addDoubleParameter('gammaBreak',0.0)
conf.addDoubleParameter('gammaMax',1.0e4)
conf.addDoubleParameter('eleK',0.0)
conf.addDoubleParameter('mBH',1.0)
conf.addDoubleParameter('eEle',0.3)
conf.addDoubleParameter('eDiss',0.3)
conf.addDoubleParameter('eDisk',0.3)
conf.addDoubleParameter('eJet',0.3)
conf.addDoubleParameter('mDot',1.0)
conf.addDoubleParameter('Gamma',10.0)
conf.addDoubleParameter('sigmaB',0.1)
conf.addDoubleParameter('sigma',0.1)
conf.addDoubleParameter('NeNp',1.0)
conf.addDoubleParameter('blrGurext',0.0)
conf.addDoubleParameter('blrGuk',3.0)
conf.addDoubleParameter('blrGue',10.0)
conf.addDoubleParameter('blrGucf',0.1)
conf.addDoubleParameter('dtGue',0.6203)
conf.addDoubleParameter('dtGurext',0.0)
conf.addDoubleParameter('dtGuk',2.0)
conf.addDoubleParameter('dtGucf',0.1)
conf.addDoubleParameter('accdGurext',0.0)
conf.addIntParameter('SABS',1)
conf.addDoubleParameter('AdiabaticABG',0.66667)
conf.addIntParameter('saveElectrons',1)
conf.addIntParameter('saveElectronsAvg',1)
conf.addIntParameter('saveLum',1)
conf.addIntParameter('saveLumPoint',1)
conf.addIntParameter('saveLumPointAvg',1)
conf.addIntParameter('saveUpeVsR',1)
conf.addIntParameter('saveRmVsR',1)
conf.addIntParameter('saveExtPlVsRm',1)
conf.addIntParameter('root',0)
conf.addIntParameter('printRoot',0)
conf.addIntParameter('KN',0)
conf.addIntParameter('Nele',1000)
conf.addIntParameter('Nsyn',1000)
conf.addIntParameter('Nssc',1000)
conf.addIntParameter('NblrGu',1000)
conf.addIntParameter('NdtGu',1000)
conf.addIntParameter('NaccdGu',1000)
conf.addDoubleParameter('R0',1.0e15,1.0e20,100.0,'log')
conf.addDoubleParameter('RInjMax',2.0e15,2.0e20,100.0,'log')
conf.addDoubleParameter('RMax',1.0e16,1.0e21,100.0,'log')
conf.create()
conf.printSummary()
