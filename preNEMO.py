#   Author: Christopher Bull. 
#   Affiliation:  Department of Geography and Environmental Sciences, 
#                 Northumbria University, Newcastle upon Tyne, UK
#   Contact: christopher.bull@northumbria.ac.uk
#   Date created: Tue, 16 Mar 2021 16:37:13
#   Machine created on: SB2Vbox
#
"""
Update of NEMO pre and post proc for ARCHER2, this was based on /lus/cls01095/work/n02/n02/chbull/repos/nemo_wed_analysis/ajtoy/configs/rnemoARCHER/run_nemo_GENERIC.py

This file now ONLY does pre
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

#no longer inherited due to NEMO4 - only changes to the run aspects are modified here (e.g. start/stop/restarts etc)
#rP_nml_patch=rP.rP_nml_patch

rP_nml_patch={}

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

    ofol='/fs2/n02/n02/chbull/nemo/nemo_output/'

    lg.info("******************************************")
    lg.info("*   project "+rP_PROJ  +"                    ")
    lg.info("*   config  "+rP_CONFIG+"                    ")
    lg.info("*   case    "+rP_CASE  +"                    ")
    lg.info("*   Desc    "+rP_DESC  +"                    ")
    lg.info("*   Ofol    "+ofol  +"                    ")
    #lg.info("*   OceanCores    "+str(rP_OCEANCORES)+"     ")
    #lg.info("*   XiosCores     "+str(rP_XIOCORES)  +"     ")
    lg.info("******************************************")

    YEAR=rP_YEAR0

    #currently using python2
    #https://medium.com/techmacademy/leading-zeros-in-python-a-hidden-mystery-revealed-ee3e3065634d
    #https://stackoverflow.com/questions/13142347/how-to-remove-leading-and-trailing-zeros-in-a-string-python
    #rP_HISFTBL=float(rP_HISFTBL.lstrip("0"))

    # while [  $YEAR -lt $rP_YEAR_MAX ]; do
    while int(YEAR) <= int(rP_YEAR_MAX):
        # lg.info("Currently working on year: " + str(YEAR))

        NDAYS=5
        # NDAYS=32
        NDAYS=365
        lg.info("We are running with NDAYS: "+ str(NDAYS))

        # ##-- calculate corresponding number of time steps for NEMO:
        # NIT000=`echo "$NITENDM1 + 1" | bc`
        # NITEND=`echo "$NITENDM1 + ${NDAYS} * 86400 / ${RN_DT}" | bc`

        # RN_DT=`grep "rn_rdt " namelist_nemo_GENERIC_${rP_CONFIG} |cut -d '=' -f2 | cut -d '!' -f1 | sed -e "s/ //g"`
        nml = f90nml.read(nmlpath)
        RN_DT=nml['namdom']['rn_rdt']

        nml_patch={}
        if not os.path.exists(rP_WORKDIR+'prod_nemo.db'):
            lg.info("First NEMO run, creating ./prod_nemo.db. Starting: "+"1 "+str(rP_YEAR0)+" 0")

            NITENDM1=0
            NIT000=str(int(NITENDM1) + 1)
            NITEND=str(int(NITENDM1) + int(NDAYS) * 86400 / int(RN_DT))

            fileHandle = open (rP_WORKDIR+ 'prod_nemo.db',"w" )
            fileHandle.write("0001 "+str(rP_YEAR0)+" 0 \n")
            fileHandle.close()

            YEAR=rP_YEAR0
            MONTH=1
            DAY=1

            NRUN=1

            # warning: can't do integer with leading zeros...
            # nml_patch['namrun']={'ln_rstart':False,'cn_exp':rP_CONFIG+'_'+rP_CASE,'nn_date0':00010101,'nn_rstctl':0,'nn_it000':int(NIT000),'nn_itend':int(NITEND),'nn_leapy':30}
            nml_patch['namrun']={'ln_rstart':False,'cn_exp':rP_CONFIG+'_'+rP_CASE,'nn_rstctl':0,'nn_it000':int(NIT000),'nn_itend':int(float(NITEND)),'nn_leapy':30}
            nml_patch['namdom']={'ln_meshmask':True}

            # thickness of the top boundary layer           (Losh et al. 2008)
            #nml_patch['namsbc_isf']={'rn_hisf_tbl':rP_HISFTBL}

        else:
            # ncj's..
            # read NRUN YEAR MONTH DAY NITENDM1 << EOF
            # `tail -1 prod_nemo.db`

            #old python2 way
            #fileHandle = open(rP_WORKDIR+ 'prod_nemo.db',"r" )
            #lineList = fileHandle.readlines()

            #lineList=[r.rstrip() for r in lineList]

            ##remove empty elements
            #lineList = filter(None, lineList)
            #fileHandle.close()

            #python3 (I hope!)
            lineList=[]
            fileHandle = open(rP_WORKDIR+ 'prod_nemo.db',"r" )
            lineList.append(fileHandle.readlines()[0][:-1])
            fileHandle = open(rP_WORKDIR+ 'prod_nemo.db',"r" )
            lineList.append(fileHandle.readlines()[-1][:-1])
            fileHandle.close()

            NRUN,YEAR,NITENDM1=lineList[-1].strip().split(' ')
            lg.info("./prod_nemo.db says, run to do is: " + str(lineList[-1]))

            NIT000=str(int(NITENDM1) + 1)
            NITEND=float(str(int(NITENDM1) + int(NDAYS) * 86400 / int(RN_DT)))
            dummy,rP_YEAR0,dummy2=(lineList[0].split(' '))[0:3]
            dummy,dummy2=None,None
            nml_patch['namrun']={'ln_rstart':True,'cn_exp':rP_CONFIG+'_'+rP_CASE,'nn_rstctl':2,'nn_it000':int(NIT000),'nn_itend':int(NITEND),'cn_ocerst_indir':rP_WORKDIR+'OUTNEMO_'+str(int(NRUN)-1).zfill(4)+'/restarts/','nn_leapy':30,'cn_ocerst_in':rP_CONFIG+'_'+rP_CASE+'_'+str(NITENDM1)+'_restart'}

            # thickness of the top boundary layer           (Losh et al. 2008)
            #nml_patch['namsbc_isf']={'rn_hisf_tbl':rP_HISFTBL}

            #check here to see if the restart files exist
            rfiles=sorted(glob.glob(rP_WORKDIR+'OUTNEMO_'+str(int(NRUN)-1).zfill(4)+'/restarts/'+rP_CONFIG+'_'+rP_CASE+'_*_restart*.nc'))

            assert(rfiles!=[]),"E R R O R: Didn't find any NEMO restart files,STOP!"
	    #lg.info("Found "+str(len(rfiles)) + " restart files, e.g. "+os.path.basename(rfiles[0]))
            rfiles=None

            if int(YEAR)>int(rP_YEAR_MAX):
                lg.warning("We are stopping (as there's nothing to do, we've already done the required years)!")
                sys.exit()

            lg.info("Running from restart files: "+rP_WORKDIR+'OUTNEMO_'+str(int(NRUN)-1).zfill(4)+'/restarts/'+rP_CONFIG+'_'+rP_CASE+'_'+str(NITENDM1)+'_restart*.nc')

        # print NIT000,NITEND 

        if rP_nml_patch!={}:
            for run_management in rP_nml_patch.keys():
                nml_patch[run_management]=rP_nml_patch[run_management]
                lg.warning('We are overwriting the default namelist variables with: '+ run_management + ' '+str(rP_nml_patch[run_management]))
        else:
            lg.info('We are not overwriting the detault namelist variables.')

        lg.info("")
        lg.info("Resulting patch for namelist_ref: "+str(nml_patch))
        f90nml.patch(nmlpath, nml_patch, out_path=nmlpath+'_new')
        shutil.move(nmlpath+'_new', nmlpath)
        shutil.copy2(nmlpath, 'namelist_cfg')


        ###########################################################
        ###-- run
        
        lg.info("")
        lg.info("LAUNCHING the simulation for YEAR "+str(YEAR))
        break


    lg.info('')
    #lg.info('All done! Consider moving your restarts to the RDF, see: '+rP_WORKDIR+'mv_restarts_rdf.sh')
    lg.info('')
    localtime = time.asctime( time.localtime(time.time()) )
    lg.info("Local current time : "+ str(localtime))
    lg.info('SCRIPT ended')
