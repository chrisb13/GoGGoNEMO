#!/usr/bin/env python 
#   Author: Christopher Bull. 
#   Affiliation:  Department of Geography and Environmental Sciences, 
#                 Northumbria University, Newcastle upon Tyne, UK
#   Contact: christopher.bull@northumbria.ac.uk
#   Date created: Mon, 15 Feb 2021 09:43:25
#   Machine created on: SB2Vbox
#
#   Update (15/02/2021): modified to work on ARCHER2 with nemo3.6 STABLE 
"""
 This is the script to set-up NEMO runs

 Update (16/05/2022): removed a lot of fluff to put it on github

 Update (13/11/2023): hacked with Louis S for psmn e5
"""
import sys,os
import datetime
import contextlib as ctx
import shutil
import glob
import subprocess
import f90nml
import numpy as np

class bcolors:
    """
    stolen: https://stackoverflow.com/questions/287871/how-to-print-colored-text-in-python
    """
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def mkdir(p):
    """make directory of path that is passed"""
    try:
       os.makedirs(p)
       print("output folder: "+p+ " does not exist, we will make one.")
    except OSError as exc: # Python >2.5
       import errno
       if exc.errno == errno.EEXIST and os.path.isdir(p):
          pass
       else: raise

def mkdomcfg(rP_WORKDIR,workfol,DOMAINcfg,BFILE,rP_nml_patch,killisf):
    """function to do leg work for making domain_cfg.nc
    """
    os.mkdir(rP_WORKDIR+'/domaincfg/')
    os.chdir(rP_WORKDIR+'/domaincfg/')
    wdir=workfol+'domaincfg/'

    ifiles=sorted(glob.glob(wdir+'*'))
    assert(ifiles!=[]),"glob didn't find anything!"
    for f in ifiles:
        if os.path.isfile(f):
            shutil.copy2(f,rP_WORKDIR+'/domaincfg/'+os.path.basename(f))

    os.symlink(DOMAINcfg, 'make_domain_cfg.exe')
    if not killisf:
        os.symlink(BFILE, 'isf_draft_meter.nc')
        os.symlink(BFILE, 'bathy_meter.nc')
    else:
        shutil.copy2(BFILE,rP_WORKDIR+'/domaincfg/bathy_meter_template.nc')
        subprocess.call('module load nco;ncks -x -v isf_draft bathy_meter_template.nc bathy_meter_middle.nc;ncrename -v Bathymetry_isf,Bathymetry bathy_meter_middle.nc bathy_meter.nc;',shell=True)
        if os.path.exists(rP_WORKDIR+'/domaincfg/bathy_meter.nc'):
            print("We have sucessfully killed the ice-shelf!") 
        else:
            print("We have not sucessfully killed the ice-shelf!") 
            import pdb;pdb.set_trace()
        os.remove(rP_WORKDIR+'/domaincfg/bathy_meter_middle.nc')
        os.remove(rP_WORKDIR+'/domaincfg/bathy_meter_template.nc')

    #need bathymetry and isf file

    #if rP_nml_patch is not None:
    #    nmlpath=rP_WORKDIR+'/domaincfg/namelist_cfg'
    #    assert(os.path.exists(nmlpath)),"can't find domaincfg/namelist_cfg"
    #    print(bcolors.HEADER+ "WARNING: "+bcolors.ENDC   +'We are patching the DOMAINCFG namelist '+nmlpath+' with: '+  bcolors.HEADER +str(rP_nml_patch) +bcolors.ENDC)
    #    f90nml.patch(nmlpath, rP_nml_patch, out_path=nmlpath+'_new')
    #    shutil.move(nmlpath+'_new', nmlpath)


def mkmesh(workfol,rP_WORKDIR,DOMAINcfg,BFILE,rP_nml_patch=None,killisf=False,extdomaincfg=['','']):
    """NEMO 4 in some cases, needs you to make a domain_cfg.nc file

    Note to make the namelist_cfg that worked for NEMO 4.0.4 we had to, hacking from:
        /mnt/lustre/a2fs-work2/work/n02/n02/chbull/nemo/nemo_output/aj_ts_melt_off_dynlin_014/0001/namelist.0001
        -change jp_cfg into an integer
        -remove various elements from namrun
        -remove a few things from namdom (everything between rn_rdt and ln_crs)
        -added ln_teos10    = .true.
    
    :DOMAINcfg: executable to use
    :extdomaincfg: 3 options
        ['',''] - default, we will use parameters provided or modified via namelist (alone)
        ['hacked','/PATH/TO/pythonscripttohack_domain_cfg.nc'] - uses default but will additionally hack the domain_cfg.nc output using passed file
        ['rogue','/PATH/TO/domain_cfg.nc'] - total rogue, i.e., if you'd rather use an externally provided domaincfg (made elsewhere).
    :returns: @todo
    """
    if extdomaincfg[0]=='':
        mkdomcfg(rP_WORKDIR,workfol,DOMAINcfg,BFILE,rP_nml_patch,killisf)
    elif extdomaincfg[0]=='hacked':
        mkdomcfg(rP_WORKDIR,workfol,DOMAINcfg,BFILE,rP_nml_patch,killisf)
        print(bcolors.HEADER+ "WARNING: "+bcolors.ENDC   +'We are going to hack the domaincfg, using '+os.path.dirname(extdomaincfg[1])+'/'+bcolors.HEADER+os.path.basename(extdomaincfg[1])+ bcolors.ENDC)
    elif extdomaincfg[0]=='rogue':
        assert(os.path.exists(extdomaincfg[1])),"domain_cfg passed does not exist!" + extdomaincfg[1]
        #shutil.copy2(extdomaincfg[1],rP_WORKDIR+'/domain_cfg.nc')
        os.symlink(extdomaincfg[1], rP_WORKDIR+'/domain_cfg.nc')
        print(bcolors.HEADER+ "WARNING: "+bcolors.ENDC   +'We are using an externally provided domaincfg '+os.path.dirname(extdomaincfg[1])+'/'+bcolors.HEADER+os.path.basename(extdomaincfg[1])+ bcolors.ENDC)
    else:
        pass
        print(bcolors.WARNING + "ERROR: "  +'Do not know what to do with your passed argument for extdomaincfg'+  bcolors.ENDC )
        import pdb;pdb.set_trace()

    return



def main(workfol,rP_CONFIG,rP_CONFIG_TYPE,rP_CASE,NEMOexe,WCONFIG,BFILE,TSFILE,rP_nml_patch=None,FLXFCE='',extdomaincfg=['','']):
    """@todo: Docstring for main
    
    :DOMAINcfg: executable to use
    :extdomaincfg: 3 options
        ['',''] - default, we will use parameters provided or modified via namelist (alone)
        ['hacked','/PATH/TO/pythonscripttohack_domain_cfg.nc'] - uses default but will additionally hack the domain_cfg.nc output using passed file
        ['rogue','/PATH/TO/domain_cfg.nc'] - total rogue, i.e., if you'd rather use an externally provided domaincfg (made elsewhere).
    :returns: @todo
    """
    # check here that the symlinks exist
    print("Check forcing files (bathymetry/ts/executable) exist..")
    assert(os.path.exists(BFILE)),"can't find file: "+BFILE
    assert(os.path.exists(TSFILE)),"can't find file: "+TSFILE
    assert(os.path.exists(NEMOexe)),"can't find file: "+NEMOexe
    print(".. Forcing files exist!")


    RDIR="/work/n02/n02/chbull/nemo/run"
    RDIR="/work/n02/n02/asmou/NEMORC/cfgs/ROSS025_4.2/TESTING"
    RDIR="/scratch/E5N/lsaddier/nemo/run"
    rP_WORKDIR=RDIR+'/'+rP_CONFIG+'_'+rP_CASE

    rP_PROJ='n02-PROPHET'
    #change me alethea to tippaccs (check spelling)
    rP_PROJ='n02-PROPHET'

    rP_OCEANCORES=20
    if rP_CONFIG_TYPE=="ASF":
        rP_OCEANCORES=1024
    elif rP_CONFIG_TYPE=="ROSS025":
        #change me alethea
        rP_OCEANCORES=300
    elif rP_CONFIG_TYPE=="LSisomipplus":
        rP_OCEANCORES=20

    #avoid weird char'
    ## dodgy hack for desc
    #import xarray as xr
    #ifile=xr.open_dataset(FLXFCE[:-3]+'T.nc')
    #rP_DESC=ifile.attrs['experiment'] + '. Now with ice shelf turned off'
    #rP_DESC=ifile.attrs['experiment'] 

    rP_DESC='RUND ESCRIPTION'
    #change me alethea
    rP_DESC='RUND ESCRIPTION'
    rP_DESC='RUN DESCRIPTION. Second run using GoGGoNEMO now with netcdf4 installed in LouisS conda env, hopefully this time pre and postnemo work'
    rP_DESC='RUN DESCRIPTION. Third run using GoGGoNEMO now with netcdf4 installed in LouisS conda env, hopefully this time pre and postnemo work. Now trying out goNEMOlong..'
    rP_DESC='RUN DESCRIPTION. Third run using GoGGoNEMO now with netcdf4 installed in LouisS conda env, hopefully this time pre and postnemo work. Now trying out goNEMOlong.. Testing the queue limit length, ask then beg for forgiveness'
    rP_DESC='RUN DESCRIPTION. Third run using GoGGoNEMO now with netcdf4 installed in LouisS conda env, hopefully this time pre and postnemo work. Now trying out goNEMOlong..  Now having fixed rebuidnemo executable link'
    rP_DESC='RUN DESCRIPTION. Third run using GoGGoNEMO now with netcdf4 installed in LouisS conda env, hopefully this time pre and postnemo work. Now trying out goNEMOlong..  Now having fixed rebuidnemo executable link. Fixing, i think, cwd before the creation of time.year.step in goNEMOlong. NOw moved to rebuild_nemo.exe with a namelist file to get rid of the dependency on ksh shell. Trying again now on the restart file as well'


    rP_YEAR0=1
    rP_YEAR_MAX=20

    #rP_STOCKDIR="/nerc/n02/n02/chbull/RawData/NEMO"  #- restart and output directory on rdf

    #change me alethea - run dir for output of runs
    #rP_STOCKDIR="/work/n02/n02/chbull/nemo/nemo_output" #- restart and output directory; now that rdf is offline
    rP_STOCKDIR="/scratch/E5N/lsaddier/nemo/nemo_output" #- restart and output directory; now that rdf is offline

    #WCONFIG=/work/n02/n02/chbull/nemo/bld_configs/input_MISOMIP/NEMO_TYP
    #WCONFIG='/work/n02/n02/chbull/nemo/bld_configs/input_ajtoy'

    FORCING='/work/n01/shared/core2'

    #if rP_nml_patch!={}:
        #print("")
        #print(bcolors.HEADER+ "WARNING: "  +'We are patching the namelist with: '+  bcolors.ENDC +str(rP_nml_patch) )
        #print("")

    #make sure you've compiled this!
    #rP_RBUILD_NEMO=NEMOdir+'/TOOLS/REBUILD_NEMO/rebuild_nemo'
    #rP_RBUILD_NEMO='/mnt/lustre/a2fs-work2/work/n02/n02/chbull/nemo/models/NEMO4/tools/REBUILD_NEMO/rebuild_nemo'
    #rP_RBUILD_NEMO='/scratch/E5N/lsaddier/nemo/models/nemo_4.2.1/tools/REBUILD_NEMO/rebuild_nemo'
    rP_RBUILD_NEMO='/scratch/E5N/lsaddier/nemo/models/nemo_4.2.1/tools/REBUILD_NEMO/BLD/bin/rebuild_nemo.exe'
    rP_MKPSI='/work/n02/n02/chbull/nemo/bld_configs/input_ajtoy/ncj_psi/post_grid_UV'

    DATE=str(datetime.datetime.now())

    ##-- User's choices END

    ##############################################
    ##-- initializations

    print("****************************************************")
    print("*          NEMO SIMULATION                         " )
    print("*   project      "+rP_PROJ                              )
    print("*   config       "+rP_CONFIG                          )
    print("*   config_type  "+bcolors.WARNING+rP_CONFIG_TYPE                     +bcolors.ENDC  )
    print("*   NEMOexe      "+NEMOexe                           )
    print("*   wconfig      "+WCONFIG                           )
    print("*   case         "+bcolors.WARNING+rP_CASE                            +bcolors.ENDC  )
    print("*   desc         "+rP_DESC                              )
    print("****************************************************")

    #############################################################
    ###-- prepare run dir

    print( "  ")
    print( "NOTE: we are using wconfig:"+WCONFIG)
    print( "  ")

    ##create run folder
    if os.path.exists(rP_WORKDIR):
        print(bcolors.FAIL+"WARNING: "+rP_WORKDIR+" already exists, so we will STOP"+bcolors.ENDC  )
        sys.exit("")
    else:
        print("WARNING: "+rP_WORKDIR+" does not exist , so we will try "+bcolors.WARNING+"CREATE"+bcolors.ENDC +' it.' )
        os.mkdir(rP_WORKDIR)
        os.chdir(rP_WORKDIR)

        with ctx.closing(open(rP_WORKDIR+'/README','w')) as handle:
            handle.write('*****************************************************'+"\n")
            handle.write('*          NEMO SIMULATION                          *'+"\n")
            handle.write('*   project   '+rP_PROJ            +'                     *'+"\n")
            handle.write('*   config    '+rP_CONFIG          +'                     *'+"\n")
            handle.write('*   NEMOexe   '+NEMOexe         +'                     *'+"\n")
            handle.write('*   case      '+rP_CASE            +'                     *'+"\n")
            handle.write('*   Desc      '+rP_DESC            +'                     *'+"\n")
            handle.write('*   Date      '+DATE            +'                     *'+"\n")
            handle.write('*                                                   *'+"\n")
            handle.write('*   Users choices:                                 *'+"\n")
            handle.write('*   rP_PROJ        '+rP_PROJ        + '                      *'+"\n")
            handle.write('*   rP_CONFIG      '+rP_CONFIG      + '                      *'+"\n")
            handle.write('*   rP_CONFIG_TYPE '+rP_CONFIG_TYPE + '                      *'+"\n")
            handle.write('*   BFILE       '+BFILE       + '                      *'+'\n')
            handle.write('*   TSFILE      '+TSFILE      + '                      *'+'\n')
            if FLXFCE!='':
                handle.write('*   FLXFILE      '+FLXFCE      + '                      *'+'\n')
                #handle.write('*   ENAME      '+ifile.attrs['ename']      + '                      *'+'\n')
            handle.write('*   rP_CASE        '+rP_CASE        + '                      *'+'\n')
            handle.write('*   rP_YEAR0       '+str(rP_YEAR0)      + '                      *'+'\n')
            handle.write('*   rP_YEAR_MAX    '+str(rP_YEAR_MAX)   + '                      *'+'\n')
            handle.write('*   RDIR        '+RDIR        + '                      *'+'\n')
            handle.write('*   rP_WORKDIR     '+rP_WORKDIR     + '                      *'+'\n')
            handle.write('*   rP_STOCKDIR    '+rP_STOCKDIR    + '                      *'+'\n')
            handle.write('*   rP_CONFIG     '+rP_CONFIG     + '                      *'+'\n')
            handle.write('*   FORCING     '+FORCING     + '                      *'+'\n')
            handle.write('*   NEMOdir     '+NEMOdir     + '                      *'+'\n')
            handle.write('*****************************************************'+'\n')

        #cat << EOF > ./env_rec
        with ctx.closing(open(rP_WORKDIR+'/env_rec','w')) as handle:
            handle.write('*   Date      '+DATE            +'*'+"\n")
            handle.write('*   record of current env:      '+"\n")
            ENV= os.getenv('ENV')
            handle.write(str(ENV)+"\n")
            handle.write(''+"\n")
            handle.write('*   record of current path:      '+"\n")
            ENV= os.getenv('PATH')
            handle.write(ENV+"\n")


        #nemo and xios
        #os.symlink(src, dst)
        os.symlink(NEMOexe, 'nemo.exe')
        #os.symlink('/work/n02/n02/asmou/NEMORC/cfgs/ROSS025_4.2/EXP02/domain_cfg.nc', 'domain_cfg.nc')

        os.symlink(WCONFIG+'/domain_cfg.nc', 'domain_cfg.nc')

        #ln -s /work/n02/n02/chbull/nemo/models/XIOSv1/bin/xios_server.exe 
        #os.symlink('/work/n02/n02/chbull/nemo/models/XIOSv1_arc2/bin/xios_server.exe', 'xios_server.exe')

        ##cp template: namelists, *.xml etc
        # Copy src to dst. (cp src dst)
        # see pro/cons at 
        # https://stackoverflow.com/questions/123198/how-do-i-copy-a-file-in-python
        #ifiles=sorted(glob.glob('/home/chris/VBoxSHARED/repos/nemo_wed_analysis/ajtoy/configs/rnemoARCHER/*'))


        #change me alethea -- add all namelist ref and cfg files, AND xml files for xios into GoGGoNEMO
        ifiles=sorted(glob.glob(workfol+'*'))
        assert(ifiles!=[]),"glob didn't find anything!"
        for f in ifiles:
            if os.path.isfile(f):
                shutil.copy2(f,rP_WORKDIR+'/')


        #change me alethea -- add symlinks for forcing (see example below)

        #if rP_CONFIG_TYPE=="AJTOY":
        #    shutil.move(rP_WORKDIR+'/'+'namelist_ref_ajtoy',rP_WORKDIR+'/'+'namelist_ref')
        #    print(bcolors.WARNING + "WARNING: Using namelist_ref_ajtoy" + bcolors.ENDC)
        #elif rP_CONFIG_TYPE=="ASF":
        #    shutil.move(rP_WORKDIR+'/'+'namelist_ref_asf',rP_WORKDIR+'/'+'namelist_ref')
        #    print(bcolors.WARNING + "WARNING: Using namelist_ref_asf" + bcolors.ENDC)

        #    print(bcolors.WARNING + "WARNING: Using xml file and field defs that will output namdyntrd" + bcolors.ENDC)
        #    shutil.move(rP_WORKDIR+'/'+'file_def_nemo-oce_asfmo.xml',rP_WORKDIR+'/'+'file_def_nemo-oce.xml')
        #    shutil.move(rP_WORKDIR+'/'+'field_def_nemo-oce_asfmo.xml',rP_WORKDIR+'/'+'field_def_nemo-oce.xml')

        #    #print(bcolors.WARNING + "WARNING: Using xml file and field defs that will output spin up diogs" + bcolors.ENDC)
        #    #shutil.move(rP_WORKDIR+'/'+'file_def_nemo-oce_spin.xml',rP_WORKDIR+'/'+'file_def_nemo-oce.xml')

        #    os.mkdir(rP_WORKDIR+'/flxfce/')
        #    if FLXFCE!='':
        #        print()
        #        print(bcolors.WARNING + "WARNING: Using flxfce " +  FLXFCE + bcolors.ENDC)
        #        print()
        #        assert(os.path.exists(FLXFCE[:-3]+'T.nc')),"can't find flux forcing grid_T file: "+FLXFCE[:-3]+'T.nc'
        #        assert(os.path.exists(FLXFCE[:-3]+'U.nc')),"can't find flux forcing grid_T file: "+FLXFCE[:-3]+'U.nc'
        #        assert(os.path.exists(FLXFCE[:-3]+'V.nc')),"can't find flux forcing grid_T file: "+FLXFCE[:-3]+'V.nc'

        #        for yy in np.arange(rP_YEAR0,rP_YEAR_MAX+52):
        #            yy=str(yy).zfill(4)+'.nc'
        #            os.symlink(FLXFCE[:-3]+'T.nc', rP_WORKDIR+'/flxfce/flxforce_grid_T_y'+yy)
        #            os.symlink(FLXFCE[:-3]+'U.nc', rP_WORKDIR+'/flxfce/flxforce_grid_U_y'+yy)
        #            os.symlink(FLXFCE[:-3]+'V.nc', rP_WORKDIR+'/flxfce/flxforce_grid_V_y'+yy)
        #else:
        #    os.remove(rP_WORKDIR+'/'+'namelist_ref_ajtoy')

        if rP_nml_patch is not None:
            nmlpath=rP_WORKDIR+'/namelist_ref'
            assert(os.path.exists(nmlpath)),"can't find namelist_ref"
            print(bcolors.HEADER+ "WARNING: "+bcolors.ENDC   +'We are patching the main NEMO namelist '+nmlpath+' with: '+  bcolors.HEADER +str(rP_nml_patch) +bcolors.ENDC + " (this shortly gets copied to be namelist_cfg, which then gets patched by preNEMO.py ...)")
            f90nml.patch(nmlpath, rP_nml_patch, out_path=nmlpath+'_new')
            shutil.move(nmlpath+'_new', nmlpath)

        #os.symlink('namelist_ref', 'namelist_cfg')
        #NEMO4 doesn't seem to like this! (symlink)
         #STOP
           #===>>>> : bad opening file: namelist_cfg

        shutil.copy2(rP_WORKDIR+'/'+'namelist_ref',rP_WORKDIR+'/'+'namelist_cfg')

        print(bcolors.WARNING + "WARNING: Using bathymetry and ice-shelf draft" + bcolors.ENDC+" from: "+BFILE)
        print("WARNING: ")

        os.symlink(BFILE, 'isf_draft_meter.nc')
        os.symlink(BFILE, 'bathy_meter.nc')

        print("WARNING: ")
        print(bcolors.WARNING + "WARNING: Using restoring and non-standard temp init"+ bcolors.ENDC+": "+TSFILE)
        os.symlink(TSFILE, 'TS_init.nc')
        os.symlink(TSFILE, 'resto.nc')
        print("WARNING: ")


        #ln -s /work/n02/n02/chbull/nemo/bld_configs/input_MISOMIP/NEMO_TYP/nemo_base_WARM-NEWFIX.nc TS_init.nc

    #file that will be imported by GoGoNEMO.py for all run parameters
    with ctx.closing(open(rP_WORKDIR+'/rPARAMS.py','w')) as handle:
            handle.write(''+'\n')

            handle.write('"""'+'\n')
            handle.write("This is the script that captures all the key run parameters and is imported by GoGoNEMO.py"+"\n")
            handle.write('"""'+'\n')
            handle.write(''+'\n')

            handle.write("rP_OCEANCORES='"  +str(rP_OCEANCORES)+"'"+'\n')
            handle.write("rP_STOCKDIR='"  +str(rP_STOCKDIR)+"'"+'\n')
            handle.write("rP_WORKDIR='"  +str(rP_WORKDIR)+"'"+'\n')
            handle.write("rP_RBUILD_NEMO='"  +str(rP_RBUILD_NEMO)+"'"+'\n')
            handle.write("rP_MKPSI='"  +str(rP_MKPSI)+"'"+'\n')
            handle.write("rP_PROJ='"  +str(rP_PROJ)+"'"+'\n')
            handle.write("rP_CONFIG='"  +str(rP_CONFIG)+"'"+'\n')
            handle.write("rP_CASE='"  +str(rP_CASE)+"'"+'\n')
            handle.write("rP_DESC='"  +str(rP_DESC)+"'"+'\n')
            handle.write("rP_YEAR0='"  +str(rP_YEAR0)+"'"+'\n')
            handle.write("rP_YEAR_MAX='"  +str(rP_YEAR_MAX)+"'"+'\n')

            #handle.write("rP_nml_patch="  +str(rP_nml_patch)+" "+'\n')

            handle.write(''+'\n')
            handle.write("if __name__ == '__main__':" +"\n")
            handle.write("    print('This script is designed to be imported...')"+"\n")




    with ctx.closing(open(rP_WORKDIR+'/goNEMOquick.sh','w')) as handle:
        handle.write('#!/bin/bash'+'\n')
        handle.write('#SBATCH --job-name=TestISOMIP_quick'+'\n')
        handle.write('#SBATCH -o ./%x.%j.%N.out           # output file'+'\n')
        handle.write('#SBATCH -e ./%x.%j.%N.err           # errors file'+'\n')
        handle.write('#'+'\n')
        handle.write('#SBATCH -p E5'+'\n')
        handle.write('#SBATCH --nodes=2'+'\n')
        handle.write('#SBATCH --ntasks-per-node=10'+'\n')
        handle.write('#SBATCH --cpus-per-task=1'+'\n')
        handle.write('#SBATCH --time=0-00:30:00           # day-hours:minutes:seconds'+'\n')
        handle.write('#SBATCH --mem-per-cpu=2G'+'\n')
        handle.write('###SBATCH --ntasks=20'+'\n')
        handle.write('#SBATCH --exclusive'+'\n')
        handle.write('#'+'\n')
        handle.write('echo "The job ${SLURM_JOB_ID} is running on these nodes:"'+'\n')
        handle.write('echo ${SLURM_NODELIST}'+'\n')
        handle.write('echo'+'\n')
        handle.write('#'+'\n')
        handle.write('cd $SLURM_SUBMIT_DIR    # go to the work / submission directory'+'\n')
        handle.write('#'+'\n')
        handle.write('#'+'\n')
        handle.write('module purge'+'\n')
        handle.write('module use /applis/PSMN/debian11/E5/modules/all/'+'\n')
        handle.write('module load netCDF-Fortran/4.6.1-gompi-2023a'+'\n')

        handle.write('export PYTHONPATH=/home/lsaddier/miniconda3/pkgs;export PATH=/home/lsaddier/miniconda3/bin:$PATH;source activate root'+'\n')

        handle.write(''+'\n')

        handle.write('#'+'\n')
        handle.write("#This is the script allows you to run NEMO real quick..."+"\n")
        handle.write(''+'\n')

        handle.write('python preNEMO.py'+'\n')
        handle.write('mpirun -np $SLURM_NTASKS ./nemo.exe'+'\n')
        handle.write('python postNEMO.py'+'\n')


        #CB old version, only works on ARCHER (UK)
        #handle.write('#!/bin/bash'+'\n')
        #handle.write('#SBATCH --qos=short'+'\n')
        #handle.write('#SBATCH --job-name=nemo_test'+'\n')
        #handle.write('#SBATCH --time=00:20:00'+'\n')
        #if rP_CONFIG_TYPE=="ASF":
        #    handle.write('#SBATCH --nodes=8'+'\n')
        #elif rP_CONFIG_TYPE=="ROSS025":
        #    handle.write('#SBATCH --nodes=3'+'\n')
        #else:
        #    handle.write('#SBATCH --nodes=1'+'\n')

        #handle.write('#SBATCH --ntasks='+str(rP_OCEANCORES)+'\n')
        #handle.write('#SBATCH --account=n02-PROPHET'+'\n')
        #handle.write('#SBATCH --partition=standard'+'\n')
        #handle.write('module restore'+"\n")
        #handle.write('module load cray-hdf5-parallel'+"\n")
        #handle.write('module load cray-netcdf-hdf5parallel'+"\n")
        #handle.write('module load xpmem'+"\n")
        #handle.write('module load perftools-base'+"\n")
        #handle.write('export OMP_NUM_THREADS=1'+'\n')
        #handle.write('export PYTHONPATH=/work/n02/n02/chbull/anaconda3/pkgs;export PATH=/work/n02/n02/chbull/anaconda3/bin:$PATH;source activate root'+'\n')
        #handle.write('#'+'\n')
        #handle.write("#This is the script allows you to run NEMO real quick..."+"\n")
        #handle.write(''+'\n')


        #handle.write('python preNEMO.py'+'\n')
        #if rP_CONFIG_TYPE=="ASF":
        #    handle.write('srun --ntasks='+str(rP_OCEANCORES)+' ./nemo.exe'+'\n')
        #else:
        #    handle.write('srun --ntasks='+str(rP_OCEANCORES)+' --mem-bind=local --cpu-bind=v,map_cpu:00,0x1,0x2,0x3,0x4,0x5,0x6,0x7,0x8,0x9,0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,0x18,0x19,0x20,0x21,0x22,0x23, ./nemo.exe'+'\n')
        #handle.write('python postNEMO.py'+'\n')
        #handle.write('#or more simply'+'\n')
        #handle.write('#srun --distribution=block:block --hint=nomultithread -n 4 ./nemo.exe'+'\n')

        #handle.write('#rebuild outputs '+'\n')
        #handle.write('#/mnt/lustre/a2fs-work2/work/n02/n02/chbull/nemo/models/NEMO4/tools/REBUILD_NEMO/rebuild_nemo mesh_mask '+str(rP_OCEANCORES)+'\n')


    subprocess.call('chmod u+x '+rP_WORKDIR+'/goNEMOquick.sh',shell=True)

    with ctx.closing(open(rP_WORKDIR+'/goNEMOlong.sh','w')) as handle:
        handle.write('#!/bin/bash'+'\n')
        handle.write('#SBATCH --job-name=TestISOMIP_long'+'\n')
        handle.write('#SBATCH -o ./%x.%j.%N.out           # output file'+'\n')
        handle.write('#SBATCH -e ./%x.%j.%N.err           # errors file'+'\n')
        handle.write('#'+'\n')
        handle.write('#SBATCH -p E5'+'\n')
        handle.write('#SBATCH --nodes=2'+'\n')
        handle.write('#SBATCH --ntasks-per-node=10'+'\n')
        handle.write('#SBATCH --cpus-per-task=1'+'\n')
        handle.write('#SBATCH --time=0-10:00:00           # day-hours:minutes:seconds'+'\n')
        handle.write('#SBATCH --mem-per-cpu=2G'+'\n')
        handle.write('###SBATCH --ntasks=20'+'\n')
        handle.write('#SBATCH --exclusive'+'\n')
        handle.write('#'+'\n')
        handle.write('echo "The job ${SLURM_JOB_ID} is running on these nodes:"'+'\n')
        handle.write('echo ${SLURM_NODELIST}'+'\n')
        handle.write('echo'+'\n')
        handle.write('#'+'\n')
        handle.write('cd $SLURM_SUBMIT_DIR    # go to the work / submission directory'+'\n')
        handle.write('#'+'\n')
        handle.write('#'+'\n')
        handle.write('module purge'+'\n')
        handle.write('module use /applis/PSMN/debian11/E5/modules/all/'+'\n')
        handle.write('module load netCDF-Fortran/4.6.1-gompi-2023a'+'\n')

        handle.write('export PYTHONPATH=/home/lsaddier/miniconda3/pkgs;export PATH=/home/lsaddier/miniconda3/bin:$PATH;source activate root'+'\n')

        handle.write(''+'\n')

        handle.write('#'+'\n')
        handle.write("#This is the script allows you to run NEMO one job at a time..."+"\n")
        handle.write('echo "Current Date and Time is (start): " '+'`date`'+"\n")
        handle.write(''+'\n')

        #if extdomaincfg[0]=='' or extdomaincfg[0]=='hacked':
        #    handle.write('#create domain_cfg.nc'+'\n')
        #    handle.write('cd '+'domaincfg'+"\n")
        #    handle.write('srun -n 1 ./make_domain_cfg.exe '+"\n")
        #    handle.write("if [ -f domain_cfg.nc ]; then"+"\n")
        #    handle.write('   mv domain_cfg.nc ../'+"\n")
        #    handle.write('else'+"\n")
        #    handle.write('   echo "ERROR: domain_cfg NOT created, stopping."'+"\n")
        #    handle.write('   exit 1'+"\n")
        #    handle.write("fi"+"\n")
        #    handle.write('cd '+'..'+"\n")
        #    handle.write(''+'\n')
        #    if extdomaincfg[0]=='hacked':
        #        handle.write('echo "Current Date and Time is (for domaincfg hacking timer): " '+'`date`'+"\n")
        #        handle.write('python ' +extdomaincfg[1] +' '+rP_WORKDIR+'/domain_cfg.nc '+rP_WORKDIR+'/domain_cfg_hckd.nc '+'\n')
        #        handle.write("if [ -f domain_cfg_hckd.nc ]; then"+"\n")
        #        handle.write('   echo "YAY: hacked domain_cfg created."'+"\n")
        #        handle.write('   rm domain_cfg.nc'+"\n")
        #        handle.write('   mv domain_cfg_hckd.nc domain_cfg.nc'+"\n")
        #        handle.write('else'+"\n")
        #        handle.write('   echo "ERROR: hacked domain_cfg NOT created, stopping."'+"\n")
        #        handle.write('   exit 1'+"\n")
        #        handle.write("fi"+"\n")
        #        handle.write('echo "Current Date and Time is (for domaincfg hacking timer): " '+'`date`'+"\n")
        #        handle.write(''+"\n")


        handle.write(' '+"\n")
        #bug? cwd was empty on e5 (and archer?!)
        #handle.write('cd $cwd'+"\n")

        handle.write('    pwd'+"\n")
        handle.write('    cd '+rP_WORKDIR+"\n")
        handle.write('    echo "Current directory is:"'+"\n")
        handle.write('    pwd'+"\n")

        handle.write("if [ -f time.year.step ]; then"+"\n")
        handle.write("   source time.year.step"+"\n")
        handle.write('else'+"\n")
        handle.write("   echo 'year=1' > time.year.step"+"\n")
        handle.write("   source time.year.step"+"\n")
        handle.write(' '+"\n")
        handle.write("fi"+"\n")
        handle.write('echo "Running year: "'+"\n")
        handle.write('echo "Running year: "$year" out of '+str(rP_YEAR_MAX)+'"'+"\n")
        handle.write('echo " "'+"\n")
        handle.write("while [ $year -lt "+str(int(rP_YEAR_MAX)+1)+" ]"+"\n")
        #handle.write('for i in {1..'+rP_YEAR_MAX+'}'+"\n")
        #handle.write('for i in {1..'+rP_YEAR_MAX+'}'+"\n")
        handle.write('do'+"\n")

        handle.write('    echo "Current Date and Time is ("$year" start): " '+'`date`'+"\n")
        handle.write('    echo "Run NEMO"'+"\n")

        handle.write('    cd '+rP_WORKDIR+"\n")
        handle.write('    echo "Current directory is:"'+"\n")
        handle.write('    pwd'+"\n")
        handle.write('    python preNEMO.py'+"\n")
        #handle.write('    srun --ntasks='+str(rP_OCEANCORES)+' ./nemo.exe'+'\n')
        handle.write('    mpirun -np $SLURM_NTASKS ./nemo.exe'+'\n')
        handle.write(' '+"\n")

        handle.write('    echo "Clean up NEMO"'+"\n")
        #handle.write('    srun --ntasks='+'1'+' --tasks-per-node='+'1'+' --cpus-per-task=1 python postNEMO.py '+"\n")
        handle.write('    python postNEMO.py'+'\n')

        handle.write('    wait'+"\n")
        handle.write('    echo "Current Date and Time is ("$year" end): " '+'`date`'+"\n")

        handle.write('    '+"\n")

        handle.write('    echo "And repeat."'+"\n")

        handle.write("    year=$[$year+1]"+"\n")
        handle.write('    cd $cwd'+"\n")
        handle.write("    echo 'year='$year > time.year.step"+"\n")
        handle.write('done'+"\n")

        handle.write('echo "Current Date and Time is (end): " '+'`date`'+"\n")
        handle.write('echo "Its over..."'+"\n")

    subprocess.call('chmod u+x '+rP_WORKDIR+'/goNEMOlong.sh',shell=True)


    #so we can run lots at once
    mkdir(workfol+'rfiles/')
    if not os.path.exists(workfol+'rfiles/runme'):
        with ctx.closing(open(workfol+'rfiles/runme','w')) as handle:
            handle.write(str(rP_YEAR_MAX)+''+'\n')

    with ctx.closing(open(workfol+'rfiles/runme','a')) as handle:
        handle.write(rP_WORKDIR+'/ '+'\n')

    return rP_WORKDIR


if __name__ == "__main__":
    #for idx,offset in zip([6,7,8,9],[0.02,0.04,0.06,0.08]):
    workfol='/work/n02/n02/chbull/repos/nemo_wed_analysis/ajtoy/configs/rnemoARCHER2/'
    workfol='/work/n02/n02/asmou/NEMORC/cfgs/ROSS025_4.2/GoGGoNEMO/'
    workfol='/home/lsaddier/runnemo/'

    IFX='-30'
    CASE='04'
    CASE='01'
    #CASE='s01'
    rP_CASE=CASE.zfill(3)

    rP_CONFIG_TYPE='SVM'
    #rP_CONFIG_TYPE='AJTOY'
    rP_CONFIG_TYPE='ASF'
    rP_CONFIG_TYPE="SR_ML"
    rP_CONFIG_TYPE="ROSS025"

    rP_CONFIG_TYPE="LSisomipplus"

    mkillisf=False
    rP_nml_patch={}


    if rP_CONFIG_TYPE=="SVM":
        pass
    elif rP_CONFIG_TYPE=="AJTOY":
        pass
    elif rP_CONFIG_TYPE=="ASF":
        pass
    elif rP_CONFIG_TYPE=="SR_ML":
        pass
    elif rP_CONFIG_TYPE=="ROSS025":
        NUM=1
        rP_CONFIG='expstuff'+str(NUM).zfill(3)
        WCONFIG='/mnt/lustre/a2fs-work2/work/n02/shared/shrr/output_tests_v2'
        NEMOdir='/mnt/lustre/a2fs-work2/work/n02/n02/chbull/nemo/models/NEMO4/'
        NEMOexe='/mnt/lustre/a2fs-work2/work/n02/n02/chbull/nemo/models/NEMO4/tests/slopeVmelt/BLD/bin/nemo.exe'
        #sebastian's wacky geometries
        #BFILE='/mnt/lustre/a2fs-work2/work/n02/shared/chbull/SR_64x64_test/bathy_meter_Exp00003.nc'
        #TSFILE='/mnt/lustre/a2fs-work2/work/n02/shared/chbull/SR_64x64_test/TS_init_Exp00003.nc'
        BFILE =WCONFIG+'/bathy_meter_ROSS025.nc'
        TSFILE=WCONFIG+'/istate_TS_ROSS025.nc'

        WCONFIG='/work/n02/n02/asmou/NEMORC/cfgs/ROSS025_4.2/TESTING'
        NEMOdir='/work/n02/n02/asmou/NEMORC/cfgs/ROSS025_4.2/EXP02/'
        NEMOexe='/work/n02/n02/asmou/NEMORC/cfgs/ROSS025_4.2/BLD/bin/nemo.exe'
        BFILE =WCONFIG+'/bathy_meter_ROSS025.nc'
        TSFILE=WCONFIG+'/istate_TS_ROSS025.nc'
        rP_nml_patch=None
        rP_nml_domcfgpatch={}
        rP_nml_domcfgpatch['namcfg']={'jpiglo':int(64),'jpidta':int(64),'jpjglo':int(64),'jpjdta':int(64)}
        rP_nml_domcfgpatch['namdom']={'pphmax':float(2000.0)}
        flxfce=''
        customdomcfg=['','']
    elif rP_CONFIG_TYPE=="LSisomipplus":
        NUM=4
        #NUM=idx
        NUM=14
        rP_CONFIG='iso'+str(NUM).zfill(3)

        WCONFIG='/scratch/E5N/lsaddier/nemo/bld_configs/isomipplus'
        NEMOdir='/scratch/E5N/lsaddier/nemo/models/nemo_4.2.1/'
        NEMOexe='/scratch/E5N/lsaddier/nemo/models/nemo_4.2.1/tests/MY_ISOMIP+_2/BLD/bin/nemo.exe'

        BFILE =WCONFIG+'/isomip+_NEMO_242_geom_ocean4.nc'
        TSFILE=WCONFIG+'/nemo_base_WARM.nc'

        #LS: modify NEMO namelist 
        rP_nml_patch=None
        #rP_nml_patch={}
        #rP_nml_patch['namisf']={'rn_gammat0':0.0215,'cn_gammablk':'spe'}
        #rP_nml_patch['namisf']={'rn_gammat0':0.0215 + offset}

        rP_nml_domcfgpatch={}
        rP_nml_domcfgpatch['namcfg']={'jpiglo':int(64),'jpidta':int(64),'jpjglo':int(64),'jpjdta':int(64)}
        rP_nml_domcfgpatch['namdom']={'pphmax':float(2000.0)}
        flxfce=''
        customdomcfg=['','']
    else:
        print(bcolors.FAIL + "We are using config_type: "+rP_CONFIG_TYPE + bcolors.ENDC)
        print(bcolors.FAIL + "E R R O R: I don't know what to do with this config type."+ bcolors.ENDC)
        sys.exit()

    print(bcolors.WARNING + "WARNING: "  +'Hard overwriting'+  bcolors.ENDC + ' of NEMO version')

    rP_WORKDIR=main(workfol,rP_CONFIG,rP_CONFIG_TYPE,rP_CASE,NEMOexe,WCONFIG,BFILE,TSFILE,FLXFCE=flxfce,rP_nml_patch=rP_nml_patch,extdomaincfg=customdomcfg)

    #domaincfg='/scratch/E5N/lsaddier/nemo/models/nemo_4.2.1/tools/DOMAINcfg/BLD/bin/make_domain_cfg.exe'
    #mkmesh(workfol,rP_WORKDIR,domaincfg,BFILE,rP_nml_patch=rP_nml_domcfgpatch,killisf=mkillisf,extdomaincfg=customdomcfg)
    #subprocess.call('cd '+rP_WORKDIR+' ; sbatch '+rP_WORKDIR+'/goNEMOquick.sh',shell=True)
    subprocess.call('cd '+rP_WORKDIR+'; sbatch '+rP_WORKDIR+'/goNEMOlong.sh',shell=True)
