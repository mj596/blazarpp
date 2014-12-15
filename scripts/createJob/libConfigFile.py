# linked parameters need to have the same size!

class configFile:
    import numpy as np
    import math
    import sys
    import os

    def __init__(self, _projectName):
        self.projectName = _projectName
        self.info = []
        self.addInfo('-------\nSummary\n-------')
        self.addInfo('Project name '+self.projectName)
        self.options = []
        self.outdir = './'
        self.overwrite = False
        self.params_text = []
        self.params = []
        self.params_single = []
        self.params_linked = []
        self.params_iter = []
        self.log = True

    def setExe(self, _pathToExe):
        self.exe = _pathToExe
        self.addInfo('executable file '+self.exe)

    def setAnalyseExe( self, _path_to_analyse_exe ):
        self.analyse_exe = _path_to_analyse_exe
        self.addInfo( 'analyse executable file '+self.analyse_exe)
    
    def printSummary(self):
        for entry in self.info:
            print entry

    def addInfo(self,_info):
        self.info.append(_info)

    def addOption(self,_option):
        self.options.append( _option )
        self.addInfo('option '+_option)

    def setNoLog(self):
        self.log = False

    def setDir(self,_outdir,_overwrite):
        self.outdir = _outdir+'/'+self.projectName
        self.overwrite = _overwrite
        self.addInfo('project directory '+self.outdir)
        self.addInfo(' out '+self.outdir+'/out')
        self.addInfo(' qsub '+self.outdir+'/qsub')
        self.addInfo(' log '+self.outdir+'/log')
        self.addInfo(' conf '+self.outdir+'/conf')
        if self.overwrite:
            self.addInfo('Existing files will be overwritten.')

    def addStringParameter(self,_name,_val):
        self.params_text.append([_name,_val])
        self.addInfo('parameter '+_name+' ... '+_val)

    def addDoubleParameter(self,_name,_val,_val2=0,_n=0,_scale='log',_mode='link'):
        if _scale != 'lin' and _scale != 'log':
            # scale is wrong; skip that parameter and move on
            print '*** Parameter',_name,'Please use \'lin\' or \'log\' scale.'
            print 'Skip parameter.'
        else:
            # scale is ok; continue
            if _val2 == 0 and _n == 0:
                # is single
                self.addInfo('parameter '+_name+' ... '+str(_val))
                self.params_single.append( _name )
                # add to main params list
                array_range = [0]
                self.params.append([_name,float(_val),float(_val2),float(_n),_scale,_mode,array_range])
            else:
                # is range
                self.addInfo('parameter '+_name+' ... '+_scale+' scale range ('+str(_val)+'-'+str(_val2)+';'+str(_n)+')')
                # add to main params list
                array_range = self.np.logspace(self.math.log10(_val),self.math.log10(_val2),_n)
                self.params.append([_name,float(_val),float(_val2),float(_n),_scale,_mode,array_range])
                # is iter
                if _mode == 'iter':
                    self.params_iter.append(_name)
                # is liked
                elif _mode == 'link':
                    self.params_linked.append(_name)
                else:
                    print '*** Parameter',_name,'Please use \'iter\' or \'link\' mode.'
                    print 'Skip parameter.'

# ---
    def addIntParameter(self,_name,_val,_val2=0,_n=0,_scale='log',_mode='link'):
        if _scale != 'lin' and _scale != 'log':
            # scale is wrong; skip that parameter and move on
            print '*** Parameter',_name,'Please use \'lin\' or \'log\' scale.'
            print 'Skip parameter.'
        else:
            # scale is ok; continue
            if _val2 == 0 and _n == 0:
                # is single
                self.addInfo('parameter '+_name+' ... '+str(_val))
                self.params_single.append( _name )
                # add to main params list
                array_range = [0]
                self.params.append([_name,int(_val),int(_val2),int(_n),_scale,_mode,array_range])
            else:
                # is range
                self.addInfo('parameter '+_name+' ... '+_scale+' scale range ('+str(_val)+'-'+str(_val2)+';'+str(_n)+')')
                # add to main params list
                array_range = self.np.logspace(self.math.log10(_val),self.math.log10(_val2),_n)
                self.params.append([_name,int(_val),int(_val2),int(_n),_scale,_mode,array_range])
                # is iter
                if _mode == 'iter':
                    self.params_iter.append(_name)
                # is liked
                elif _mode == 'link':
                    self.params_linked.append(_name)
                else:
                    print '*** Parameter',_name,'Please use \'iter\' or \'link\' mode.'
                    print 'Skip parameter.'


# ---


    def getParam(self,_name):
        for param in self.params:
            if param[0] == _name:
                return param

    def create(self):
        print '----------------\nCreate config...\n----------------'
        # first check if linked params have the same size
        if len(self.params_linked):
            params_linked_temp_size = self.getParam(self.params_linked[0])[3]
            for param in self.params_linked:
                if self.getParam(param)[3] != params_linked_temp_size:
                    print 'Linked parameters have different number of steps! Quit.'
                    self.sys.exit(1)

        def createDirs(self):
            os.makedirs(self.outdir)
            os.makedirs(self.outdir+'/out')
            os.makedirs(self.outdir+'/qsub')
            os.makedirs(self.outdir+'/conf')
            os.makedirs(self.outdir+'/log')
            print 'Done creating dirs.'

        print 'Creating directories ...',
        import os
        import sys
        if os.path.exists(self.outdir):
            print self.outdir, 'already exists.',
            if self.overwrite:
                print 'Overwritting.'
                import shutil
                shutil.rmtree(self.outdir)
                createDirs(self)
            else:
                print ' Quit.'
                sys.exit()
        else:
            createDirs(self)

        print 'Creating config files...'
        ranges=[]
        iter_names=[]
        print 'Used parameters'
        print 'Text',
        for param in self.params_text:
            print param[0],
        print ''        
        print 'Single',
        for param in self.params_single:
            print self.getParam(param)[0],
        print ''
        print 'Iter',
        for param in self.params_iter:
            print self.getParam(param)[0],
            if self.getParam(param)[3]:
                print '(range)',
            ranges.append(tuple(self.getParam(param)[6]))
        print ''
        print 'Linked',
        for param in self.params_linked:
            print self.getParam(param)[0],
            if self.getParam(param)[3]:
                print '(range)',
#            ranges.append(tuple(self.getParam(param)[6]))
        print ''

        ranges=tuple(ranges)
        import loop_rec
        iterations=[]
        loop_rec.loop_rec(ranges,iterations)

        # array for qsub_all.qsub
        runs_array = []

        # Make main loop over iter and linked parameter values
#        print "Parameter sets:"        
        for i in range(len(self.getParam(self.params_linked[0])[6])):
            data_linked=[]
            filename_linked = 'parameters'
            for linked in self.params_linked:
                filename_linked += '_'+linked+'_'+str((self.getParam(linked)[6])[i])
                data_linked.append([linked,str((self.getParam(linked)[6])[i])])
            for iter in iterations:
                filename_iter=''
                data_iter=[]
                for j in range(len(self.params_iter)):
                    filename_iter += '_'+self.params_iter[len(self.params_iter)-j-1]+'_'+str(iter[j])
                    data_iter.append([self.params_iter[len(self.params_iter)-j-1],str(iter[j])])
                filename = filename_linked+filename_iter
#                print filename
#                for i_linked in data_linked:
#                    print i_linked[0], i_linked[1]
#                for i_iter in data_iter:
#                    print i_iter[0], i_iter[1]
                # Here we can start writing files
                # qsub file
                qsub_filename = self.outdir+'/qsub'+'/'+filename+'.qsub' 
                file=open(qsub_filename,'w')
                runs_array.append(filename)

                for item in self.options:
                    file.write(item+'\n')

                file.write('source ~/.bashrc\n')
                file.write('cd '+self.outdir+'\n')
                if self.log == True:
                    file.write(self.exe+' '+self.outdir+'/conf'+'/'+filename+'.conf'+' > '+self.outdir+'/log'+'/'+filename+'.log'+' 2>&1')
                else:
                    file.write(self.exe+' '+self.outdir+'/conf'+'/'+filename+'.conf'+' > /dev/null')
                file.close()
                st = self.os.stat( qsub_filename )
                self.os.chmod( qsub_filename, 0o744 )
                
                # conf file
                file=open(self.outdir+'/conf'+'/'+filename+'.conf','w')
                # set output
                output_name = 'output'+' '+self.outdir+'/out/'+filename 
                file.write(output_name+'\n')
                # add all parameters
                for item in self.params_text:
                    file.write(item[0]+' '+str(item[1])+'\n')
                for item in self.params_single:
                    file.write(self.getParam(item)[0]+' '+str(self.getParam(item)[1])+'\n')
                for i_linked in data_linked:
                    file.write(i_linked[0]+' '+i_linked[1]+'\n')
                for i_iter in data_iter:
                    file.write(i_iter[0]+' '+i_iter[1]+'\n')
                
                file.write('END')
                file.close()
        print "Done creating files. All should be ready now."

        # create qsub file
        file=open(self.outdir+'/'+'qsub_all.sh','w')
        self.os.chmod( self.outdir+'/'+'qsub_all.sh', 0o744 )
        file.write('#! /bin/bash\n')
        file.write('project_dir='+self.outdir+'\n')
        file.write('runlist=( ')
        for line in runs_array:
            file.write('\"'+line+'\" ')
        file.write(')\n')
        file.write('for run in ${runlist[@]}\n')
        file.write('do\n')
        file.write('\t qsub \"${project_dir}'+'/qsub/'+'${run}.qsub\"\n')
        file.write('done\n')
        file.close()
        print 'Qsub all with\n\t'+self.outdir+'/'+'qsub_all.qsub'

        # create list all file
        file=open(self.outdir+'/'+'list_all.sh','w')
        self.os.chmod( self.outdir+'/'+'list_all.sh', 0o744 )
        file.write('#! /bin/bash\n')
        file.write('project_dir='+self.outdir+'\n')
        file.write('runlist=( ')
        for line in runs_array:
            file.write('\"'+line+'\" ')
        file.write(')\n')
        file.write('for run in ${runlist[@]}\n')
        file.write('do\n')
        file.write('\t echo ${project_dir}'+'/out/'+'${run}\n')
        file.write('done\n')
        file.close()
        print 'List all with\n\t'+self.outdir+'/'+'list_all.sh'

        # create run file
        file=open(self.outdir+'/'+'run_all.sh','w')
        self.os.chmod( self.outdir+'/'+'run_all.sh', 0o744 )
        file.write('#! /bin/bash\n')
        file.write('project_dir='+self.outdir+'\n')
        file.write('runlist=( ')
        for line in runs_array:
            file.write('\"'+line+'\" ')
        file.write(')\n')
        file.write('for run in ${runlist[@]}\n')
        file.write('do\n')
        file.write('\t \"${project_dir}'+'/qsub/'+'${run}.qsub\"\n')
        file.write('done\n')
        file.close()
        print 'Run all with\n\t'+self.outdir+'/'+'run_all.sh'

        # create analyse file
        file=open(self.outdir+'/'+'analyse_all.sh','w')
        self.os.chmod( self.outdir+'/'+'analyse_all.sh', 0o744 )
        file.write('#! /bin/bash\n')
        file.write('rm -v qw.dat\n')
        file.write('project_dir='+self.outdir+'\n')
        file.write('runlist=( ')
        for line in runs_array:
            file.write('\"'+line+'\" ')
        file.write(')\n')
        file.write('for run in ${runlist[@]}\n')
        file.write('do\n')
        file.write('         echo ${project_dir}/out/${run}\n')
        file.write('         r=$(echo ${project_dir}/out/${run} | awk -F\'R0_\' \'{ print $2 }\' | awk -F\'_RInjMax\' \'{ print $1 }\')\n')
        file.write('         printf \"%s \" ${r} >> qw.dat\n')
        file.write('         '+self.analyse_exe+' ${project_dir}/out/${run} LvPointAvg False >> qw.dat\n')
        file.write('done\n')
        file.close()
        print 'Analyse all with\n\t'+self.outdir+'/'+'analyse_all.sh'
