#!/usr/bin/env python 
#   Author: Christopher Bull. 
#   Affiliation:  Department of Geography and Environmental Sciences, 
#                 Northumbria University, Newcastle upon Tyne, UK
#   Contact: christopher.bull@northumbria.ac.uk
#   Date created: Thu, 18 Mar 2021 16:58:46
#   Machine created on: SB2Vbox
#
"""
Update of NEMO pre and post proc for ARCHER2, this was based on /lus/cls01095/work/n02/n02/chbull/repos/nemo_wed_analysis/ajtoy/configs/rnemoARCHER/run_nemo_GENERIC.py

This file now ONLY does post
"""
from cb2logger import *
import f90nml
import os
import glob
import shutil
import re
import sys
import numpy as np
import contextlib as ctx
from netCDF4 import Dataset

import rPARAMS as rP

rP_OCEANCORES=rP.rP_OCEANCORES
#rP_XIOCORES=rP.rP_XIOCORES
#rP_NODES=rP.rP_NODES
#rP_RHOURS=rP.rP_RHOURS
#rP_DEBUGJOB=rP.rP_DEBUGJOB
rP_STOCKDIR=rP.rP_STOCKDIR
rP_WORKDIR=rP.rP_WORKDIR+'/'
rP_RBUILD_NEMO=rP.rP_RBUILD_NEMO
rP_MKPSI=rP.rP_MKPSI
rP_PROJ=rP.rP_PROJ
rP_CONFIG=rP.rP_CONFIG
rP_CASE=rP.rP_CASE
rP_DESC=rP.rP_DESC
rP_YEAR0=int(rP.rP_YEAR0)
rP_YEAR_MAX=int(rP.rP_YEAR_MAX)
#rP_ARRAYJOB=rP.rP_ARRAYJOB
#rP_HISFTBL=rP.rP_HISFTBL

def mkdir(p):
    """make directory of path that is passed"""
    try:
       os.makedirs(p)
       lg.info("output folder: "+p+ " does not exist, we will make one.")
    except OSError as exc: # Python >2.5
       import errno
       if exc.errno == errno.EEXIST and os.path.isdir(p):
          pass
       else: raise

def rebuild_mesh_mask(handle):
    """function to rebuild the mesh_mask using the NEMO tool
    :handle: text file handle
    :returns: 
    """
    handle.write(''+'\n')
    handle.write('echo "Re-combining: '+'mesh_mask'+'"'+'\n')
    handle.write('cp '+rP_RBUILD_NEMO+' '+'rebuild_nemo'+'\n')
    handle.write('cp '+rP_RBUILD_NEMO+'.exe'+' '+'rebuild_nemo.exe'+'\n')
    handle.write('srun -n 1 '+'rebuild_nemo'+' '+'mesh_mask'+' '+str(rP_OCEANCORES)+'\n')
    handle.write(''+'\n')

    handle.write('if [ -f ' +'mesh_mask.nc'+' ]; then'+'\n')
    handle.write('   echo "File: '+'mesh_mask.nc'+' reassembled ok"'+'\n')
    handle.write('   rm '+'mesh_mask_*.nc'+'\n')
    handle.write('else'+'\n')
    handle.write('   echo "!@#$% PROBLEM WITH RE-ASSEMBLY OF FILE'+' mesh_mask"'+'\n')
    handle.write('   echo ">>>>>>>>>> STOP !"'+'\n')
    handle.write('   exit 1'+'\n')
    handle.write('fi'+'\n')
    return

if __name__ == "__main__": 

    ##########
    #  init  #
    ##########

    nmlpath="namelist_ref"
    # nmlpath_ice="namelist_ice_ref"

    os.chdir(rP_WORKDIR)

    ##########
    #  init  #
    ##########

    LogStart('',fout=False)

    ofol='/work/n02/n02/chbull/nemo/nemo_output/'

    lg.info("******************************************")
    lg.info("*   project "+rP_PROJ  +"                    ")
    lg.info("*   config  "+rP_CONFIG+"                    ")
    lg.info("*   case    "+rP_CASE  +"                    ")
    lg.info("*   Desc    "+rP_DESC  +"                    ")
    lg.info("*   Ofol    "+ofol  +"                    ")
    lg.info("*   OceanCores    "+str(rP_OCEANCORES)+"     ")
    #lg.info("*   XiosCores     "+str(rP_XIOCORES)  +"     ")
    lg.info("******************************************")

    YEAR=rP_YEAR0

    NDAYS=365
    nml = f90nml.read(nmlpath)
    RN_DT=nml['namdom']['rn_rdt']

    lineList=[]
    fileHandle = open(rP_WORKDIR+ 'prod_nemo.db',"r" )
    lineList.append(fileHandle.readlines()[0][:-1])
    fileHandle = open(rP_WORKDIR+ 'prod_nemo.db',"r" )
    lineList.append(fileHandle.readlines()[-1][:-1])
    fileHandle.close()

    NRUN,YEAR,NITENDM1=lineList[-1].strip().split(' ')
    #lg.info("./prod_nemo.db says, run to do is: " + str(lineList[-1]))
    NITEND=str(int(float(str(int(NITENDM1) + int(NDAYS) * 86400 / int(RN_DT)))))

    rnemo=rP_WORKDIR+'GoGGoNEMO_'+str(NRUN).zfill(4)+'.sh'
    with ctx.closing(open(rnemo,'w')) as handle:
        if int(str(NRUN))==1:
            handle.write('#!/bin/bash '+'\n')
            handle.write('cd '+rP_WORKDIR+'\n')
            handle.write('echo "Current directory is:"'+"\n")
            handle.write('pwd'+"\n")
            handle.write(' '+"\n")

            rebuild_mesh_mask(handle)

            handle.write('if [ -f ' +'output.init_0000.nc'+' ]; then'+'\n')
            handle.write('   '+'\n')
            handle.write('   echo "Re-combining: '+'output.init'+'"'+'\n')
            handle.write('   srun -n 1 '+'rebuild_nemo'+' '+'output.init'+' '+str(rP_OCEANCORES)+'\n')
            handle.write('   '+'\n')

            handle.write('   if [ -f ' +'output.init.nc'+' ]; then'+'\n')
            handle.write('      echo "File: '+'output.init.nc'+' reassembled ok"'+'\n')
            handle.write('      rm '+'output.init_*.nc'+'\n')
            handle.write('   else'+'\n')
            handle.write('      echo "!@#$% PROBLEM WITH RE-ASSEMBLY OF FILE'+' output.init"'+'\n')
            handle.write('      echo ">>>>>>>>>> STOP !"'+'\n')
            handle.write('      exit 1'+'\n')
            handle.write('   fi'+'\n')

            handle.write('else'+'\n')
            handle.write('   echo "no output.init'+'"'+'\n')
            handle.write('fi'+'\n')

            handle.write(''+'\n')
            #handle.write('#create psi'+'\n')
            #handle.write('mkdir psi'+'\n')
            #handle.write('cp -v mesh_mask.nc psi'+'\n')

        rone=rP_CONFIG+'_'+rP_CASE+'_'+NITEND.zfill(8)+'_'+'restart.nc'
        # rone_ice=rP_CONFIG+'_'+rP_CASE+'_'+NITEND.zfill(8)+'_'+'restart_ice.nc'

        rone_star=rP_CONFIG+'_'+rP_CASE+'_'+NITEND.zfill(8)+'_'+'restart_????.nc'
        # rone_ice_star=rP_CONFIG+'_'+rP_CASE+'_'+NITEND.zfill(8)+'_'+'restart_ice_????.nc'

        for f,f_star in zip([rone],[rone_star]):
            handle.write(''+'\n')
            handle.write('echo "Re-combining: '+f+'"'+'\n')
            handle.write('mkdir tempo '+'\n')
            handle.write('mv '+f_star +' tempo/'+'\n')
            handle.write('cp '+'rebuild_nemo*' +' tempo/'+'\n')
            handle.write('cd tempo '+'\n')

            handle.write('srun -n 1 '+'rebuild_nemo'+' '+f[:-3]+' '+str(rP_OCEANCORES)+'\n')

            handle.write('if [ -f ' +f+' ]; then'+'\n')
            handle.write('   echo "File: '+f+' reassembled ok"'+'\n')
            handle.write('   mv '+f+' ..'+'\n')
            handle.write('   rm '+f_star+'\n')
            handle.write('else'+'\n')
            handle.write('   echo "!@#$% PROBLEM WITH RE-ASSEMBLY OF FILE '+f+'"\n')
            handle.write('   echo ">>>>>>>>>> STOP !"'+'\n')
            handle.write('   exit 1'+'\n')
            handle.write('fi'+'\n')

            handle.write('cd ..'+'\n')
            handle.write('rm -r tempo/'+'\n')

    subprocess.call('chmod u+x '+rnemo,shell=True)

    lg.info("Launching NEMO with the following script: ")
    lg.info(rnemo)
    lg.info("")
    subprocess.call(rnemo,shell=True) #cbnow


    # ##-- compress output files:
    
    odir=rP_WORKDIR+'OUTNEMO_'+str(NRUN).zfill(4)+'/'
    rdir=odir + 'restarts/'
    initdir=odir + 'inits/'
    if not os.path.exists(odir):
        mkdir(odir)
    if not os.path.exists(rdir):
        mkdir(rdir)

    ofiles=sorted(glob.glob(rP_WORKDIR+rP_CONFIG+'_'+rP_CASE+'_[1-5][dhmy]_*nc'))
    rfiles=sorted(glob.glob(rP_WORKDIR+rP_CONFIG+'_'+rP_CASE+'_*_restart*.nc'))

    #something bad maybe happened let's see if there's a NEMO error..
    if ofiles==[] or rfiles==[]:
        lg.error("")
        lg.error("No output or restart files, looking for a NEMO error...")
        p = subprocess.Popen("grep -A 4 '\''E R R'\'' ocean.output", stdout=subprocess.PIPE, shell=True)
        (output, err) = p.communicate()
        p_status = p.wait()
        lg.info(output)

    assert(ofiles!=[]),"E R R O R: Didn't find any NEMO output files, STOP! "+'Glob was: '+rP_WORKDIR+rP_CONFIG+'_'+rP_CASE+'_[1-5][dhmy]_*nc'
    assert(rfiles!=[]),"E R R O R: Didn't find any NEMO restart files,STOP! "+'Glob was: '+rP_WORKDIR+rP_CONFIG+'_'+rP_CASE+'_*_restart*.nc'

    #move the output/restart and initfiles (if they existed)
    for f in ofiles:
        shutil.move(f, odir+os.path.basename(f))

    for f in rfiles:
        shutil.move(f, rdir+os.path.basename(f))

    # initfiles=sorted(glob.glob(rP_WORKDIR+'output.init*.nc'))
    # if initfiles==[]: lg.warning("Didn't find any NEMO init files.")
    # if initfiles!=[]:
        # if not os.path.exists(initdir):
            # mkdir(initdir)

        # for f in initfiles:
            # shutil.move(f, initdir+os.path.basename(f))

    ##########################################################
    ##-- prepare next run if every little thing went all right (if not, we should have crashed by now via the globs)
    m = re.search(rP_WORKDIR+rP_CONFIG+'_'+rP_CASE+"_(.*)_restart",rfiles[0])

    if m:
        LAST_RESTART_NIT=m.groups()[0]
    else:
        lg.error("Couldn't work out (regexfail) last restart tstep..")

    # YEARf=`expr $YEAR + 1`
    YEARf=str(int(YEAR)+1)
    MONTHf=1
    DAYf=1


    lg.info("Last restart created at ocean time step "+LAST_RESTART_NIT)
    lg.info("  ---> writting this date in prod_nemo.db")
    lg.info(" ")
    subprocess.call("echo "+LAST_RESTART_NIT+" > restart_nit.txt", shell=True)

    #fileHandle = open (rP_WORKDIR+ 'prod_nemo.db',"r" )
    #lineList = fileHandle.readlines()
    #lineList=[r.rstrip() for r in lineList]
    #lineList = filter(None, lineList)

    with open(rP_WORKDIR+ 'prod_nemo.db') as f:
        lineList = [line.rstrip() for line in f]


    lineList[-1]=' '.join(lineList[-1].split(' ')[:-1])+' '+LAST_RESTART_NIT+'\n'

    #specify the starting point for the next run (because the last one finished okay)
    NexRUN=str(int(NRUN)+1).zfill(4)
    next_run=str(NexRUN)+' '+str(YEARf)+' '+LAST_RESTART_NIT
    lineList.append(next_run)

    fileHandle_two = open (rP_WORKDIR+ 'prod_nemo_two.db',"w" )
    for line in lineList:
        if '\n' not in line:
            line=line+'\n'
        fileHandle_two.write(line) 

    fileHandle.close()
    fileHandle_two.close()
    shutil.move(rP_WORKDIR+'prod_nemo_two.db',rP_WORKDIR+'prod_nemo.db')

    # ##########################################################
    # ##-- send outputs and restarts to storage disk
    
    # shutil.move('namelist','namelist.'+NRUN)
    shutil.copyfile(nmlpath, odir+'namelist.'+str(NRUN).zfill(4) ) #we need it for the next year..
    #shutil.copyfile(nmlpath[:-3]+'ref', odir+'namelist_ref.'+str(NRUN).zfill(4) ) #we need it for the next year..
    # shutil.copyfile(nmlpath_ice, odir+'namelist_ice_ref.'+str(NRUN).zfill(4) ) #we need it for the next year..
    shutil.move('ocean.output',odir+'ocean.output.'+str(NRUN).zfill(4))
    shutil.move(rnemo,odir+'GoGGoNEMO_'+str(NRUN).zfill(4)+'.sh')

    # update location of ofiles / rfiles/ initfiles
    ofiles=[odir+os.path.basename(f) for f in ofiles]
    rfiles=[rdir+os.path.basename(f) for f in rfiles]
    # if initfiles!=[]:
        # initfiles=[odir+os.path.basename(f) for f in initfiles]

    # #qsub compress_nemo_${NRUN}.sh
    cleannemo=rP_WORKDIR+'cnemo_'+str(NRUN).zfill(4)+'.sh'
    with ctx.closing(open(cleannemo,'w')) as handle:
        handle.write('#!/bin/bash '+'\n')
        #handle.write('#PBS -l select=serial=true:ncpus=1'+'\n')
        #handle.write('#PBS -l walltime=01:00:00'+'\n')
        #handle.write('#PBS -A n02-FISSA'+'\n')
        #handle.write('# Make sure any symbolic links are resolved to absolute path'+'\n')
        #handle.write('export PBS_O_rP_WORKDIR=$(readlink -f $PBS_O_rP_WORKDIR)'+'\n')
        #handle.write('#load modules'+'\n')
        #handle.write('module load cray-netcdf-hdf5parallel/4.4.1.1 '+'\n')
        #handle.write('module load cray-hdf5-parallel/1.10.0.1'+'\n')
        #handle.write('module swap PrgEnv-cray PrgEnv-intel'+'\n')

        handle.write(''+'\n')
        handle.write('cd '+rP_WORKDIR+'\n')

        #compression step ..
        handle.write('#compress ofiles'+ '\n')
        for f in ofiles:
            handle.write(''+ '\n')
            handle.write('#Doing file: '+os.path.basename(f)+ '\n')
            tmpf=' '+os.path.dirname(f)+'/tmp.nc'
            handle.write('nccopy -d 5 '+f+ tmpf +' \n')

            handle.write('if [ -f ' +tmpf+' ]; then'+'\n')
            handle.write('   mv -f '+tmpf+' '+f+'\n')
            handle.write('else'+'\n')
            handle.write('   echo "!@#$% PROBLEM WITH COMPRESSION OF FILE '+f+'"\n')
            handle.write('   echo ">>>>>>>>>> STOP !"'+'\n')
            handle.write('   exit'+'\n')
            handle.write('fi'+'\n')

        # handle.write('# ------ '+ '\n')
        # handle.write('#compress rfiles'+ '\n')
        # for f in rfiles:
            # handle.write(''+ '\n')
            # handle.write('#Doing file: '+os.path.basename(f)+ '\n')
            # tmpf=' '+os.path.dirname(f)+'/tmp.nc'
            # handle.write('nccopy -d 5 '+f+ tmpf +' \n')

            # handle.write('if [ -f ' +tmpf+' ]; then'+'\n')
            # handle.write('   mv -f '+tmpf+' '+f+'\n')
            # handle.write('else'+'\n')
            # handle.write('   echo "!@#$% PROBLEM WITH COMPRESSION OF FILE"'+f+'\n')
            # handle.write('   echo ">>>>>>>>>> STOP !"'+'\n')
            # handle.write('   exit'+'\n')
            # handle.write('fi'+'\n')

        handle.write(''+ '\n')
        handle.write('# ------ '+ '\n')
        #write to RDF
        handle.write('echo "Finished compressing, WRITE to now-not-RDF"'+ '\n')

        handle.write('mkdir -p '+ofol+rP_CONFIG+'_'+rP_CASE+'/'+str(YEAR).zfill(4)+'/'+'\n')
        handle.write('mv -v '+ odir + '*.nc '+ofol+rP_CONFIG+'_'+rP_CASE+'/'+str(YEAR).zfill(4)+'/'+'\n')
        handle.write('mv -v '+ odir + 'namelist* '+ofol+rP_CONFIG+'_'+rP_CASE+'/'+str(YEAR).zfill(4)+'/'+'\n')
        handle.write('mv -v '+ odir + 'ocean* '+ofol+rP_CONFIG+'_'+rP_CASE+'/'+str(YEAR).zfill(4)+'/'+'\n')
        handle.write('mv -v '+ odir+'GoGGoNEMO_'+str(NRUN).zfill(4)+'.sh '+ofol+rP_CONFIG+'_'+rP_CASE+'/'+str(YEAR).zfill(4)+'/'+'\n')

        if int(str(NRUN))==1:
            handle.write('mv -v '+ rP_WORKDIR + 'mesh_mask.nc '+ofol+rP_CONFIG+'_'+rP_CASE+'/'+'\n')
            handle.write('cp -v '+ rP_WORKDIR + 'README '+ofol+rP_CONFIG+'_'+rP_CASE+'/'+'\n')
            handle.write('cp -v '+ rP_WORKDIR + 'bathy_meter.nc '+ofol+rP_CONFIG+'_'+rP_CASE+'/'+'\n')
            handle.write('cp -v '+ rP_WORKDIR + 'TS_init.nc '+ofol+rP_CONFIG+'_'+rP_CASE+'/'+'\n')
            handle.write('cp -v '+ rP_WORKDIR + 'rPARAMS.py '+ofol+rP_CONFIG+'_'+rP_CASE+'/'+'\n')
            handle.write('cp -v '+ rP_WORKDIR + 'setupNEMO_ARC2.py '+ofol+rP_CONFIG+'_'+rP_CASE+'/'+'\n')
            handle.write('cp -v '+ rP_WORKDIR + 'domaincfg/namelist_ref '+ofol+rP_CONFIG+'_'+rP_CASE+'/domaincfg_namelist_ref'+'\n')
            handle.write('cp -v '+ rP_WORKDIR + 'domaincfg/namelist_cfg '+ofol+rP_CONFIG+'_'+rP_CASE+'/domaincfg_namelist_cfg'+'\n')
            handle.write('cp -v '+ rP_WORKDIR + 'domain_cfg.nc '+ofol+rP_CONFIG+'_'+rP_CASE+'/'+'\n')
            
            #handle.write('cp -v '+ rP_WORKDIR + 'production_nemo_ARCHER.py '+ofol+rP_CONFIG+'_'+rP_CASE+'/'+'\n')

            handle.write('   if [ -f ' +'output.init.nc'+' ]; then'+'\n')
            handle.write('      mv -v '+ rP_WORKDIR + 'output.init.nc '+ofol+rP_CONFIG+'_'+rP_CASE+'/'+str(YEAR).zfill(4)+'/'+'\n')
            handle.write('   else'+'\n')
            handle.write('      echo "no output.init.nc so not going to put it on the RDF"'+'\n')
            handle.write('   fi'+'\n')

        
    subprocess.call('chmod u+x '+cleannemo,shell=True)
    #subprocess.call('qsub '+cleannemo,shell=True)

    lg.info("Cleaning up NEMO with the following script: ")
    lg.info(cleannemo)
    lg.info("")
    subprocess.call(cleannemo,shell=True) #cbnow

    # we will keep the restarts in the run dir until we're happy with all the outputs..
    if not os.path.exists(rP_WORKDIR+'mv_restarts_rdf.sh'):
        append_write = 'w' # make a new file if not
        frestart= open(rP_WORKDIR+'mv_restarts_rdf.sh',append_write)
        frestart.write('#!/bin/bash'+'\n')
        frestart.write('#this is a manual step to be done later'+'\n')
        frestart.write(''+'\n')
        
    else:
        append_write = 'a' # append if already exists
        frestart= open(rP_WORKDIR+'mv_restarts_rdf.sh',append_write)

    frestart.write(''+'\n')
    frestart.write('mkdir -p '+ofol+rP_CONFIG+'_'+rP_CASE+'/'+str(YEAR).zfill(4)+'/restarts/'+'\n')
    frestart.write('mv -v '+rdir+'*.nc '+ofol+rP_CONFIG+'_'+rP_CASE+'/'+str(YEAR).zfill(4)+'/restarts/'+'\n')
    frestart.write('#mv -v '+ofol+rP_CONFIG+'_'+rP_CASE+'/'+str(YEAR).zfill(4)+'/restarts/ '+rdir+' # to undo \n')
    frestart.write(''+'\n')
    subprocess.call('chmod u+x '+rP_WORKDIR+'mv_restarts_rdf.sh',shell=True) 
    frestart.close()

    #kick on for another year.
    #YEAR=str(int(YEAR)+1)

    lg.info('')
    lg.info('All done! Consider moving your restarts to the RDF, see: '+rP_WORKDIR+'mv_restarts_rdf.sh')
    lg.info('')
    localtime = time.asctime( time.localtime(time.time()) )
    lg.info("Local current time : "+ str(localtime))
    lg.info('SCRIPT ended')
