# -*- coding: utf-8 -*-
"""
Functions used to read PEST control file parameter values from an Excel spreadsheet
and write the PEST control file.  Note that these functions do not check for the
internal consistency of the control file (e.g., number of COMLINE = NUMCOM),
instead leaving notification of those errors to the execution of PEST.
"""

import pandas as pd
import os

# -------------------------------------

def wl(param_list):
    '''Writes the user-specified line of parameters to the control file.  Note that
    some parameters are defined as lists of parameter values.  In this case the command
    iterates through the list.'''
    
    for iparam in param_list:
        if (pd.isnull(control_param_dict[iparam]) == True):
            fout.write('%s ' %('  '))
        else:
            fout.write('%s ' %(control_param_dict[iparam]))
    fout.write('\n')
    return

def get_count(idict):
    '''Returns the number of groups and number of observations or parameters in
    the obs_group_dict and the param_group_dict.'''
    
    ngrp = len(idict)
    icount = 0
    for igrp in idict:
        if 'PARAMS' in idict[igrp]:
            icount += len(idict[igrp]['PARAMS'])
        else:
            icount += len(idict[igrp]['OBS'])

    return ngrp,icount

def write_control_section(NPARGP,NPAR,NOBSGP,NOBS):
    '''Writes the control section. NOTE: Assumes that each parameter group
    has its own template file and each observation group has its own instruction file.'''

    fout.write('pcf\n* control data\n')
    wl(['RSTFLE','PESTMODE'])
    
    fout.write('%-10i%-10i%-10i%-10i%-10i%-10i\n'%(NPAR,NOBS,NPARGP,control_param_dict['NPRIOR'],NOBSGP,control_param_dict['MAXCOMPDIM']))
    fout.write('%-10i%-10i' %(NPARGP,NOBSGP))

    wl(['PRECIS','DPOINT','NUMCOM','JACFILE','MESSFILE'])
    wl(['RLAMBDA1','RLAMFAC','PHIRATSUF','PHIREDLAM','NUMLAM','JACUPDATE'])
    wl(['RELPARMAX','FACPARMAX','FACORIG','IBOUNDSTICK','UPVECBEND'])
    wl(['PHIREDSWH','NOPTSWITCH','SPLITSWH','DOAUI','DOSENREUSE'])
    wl(['NOPTMAX','PHIREDSTP','NPHISTP','NPHINORED','RELPARSTP','NRELPAR','PHISTOPTHRESH','LASTRUN','PHIABANDON'])
    wl(['ICOV','ICOR','IEIG','IRES','JCOSAVEITN','REISAVEITN'])
    return()

def write_aui(NPAR):
    
    fout.write('* automated user intervention\n')  

    MAXAUI = control_param_dict['mult_MAXAUI'] * NPAR
    NAUINOACCEPT = control_param_dict['mult_NAUINOACCEPT'] * MAXAUI
    
    fout.write('%-10i' %(MAXAUI))
    
    wl(['AUISTARTOPT','NOAUIPHIRAT','AUIRESTITN'])
    wl(['AUISENSRAT','AUIHOLDMAXCHG','AUINUMFREE'])
    
    fout.write('%-10.2f%-10.2f%-10i\n' %(control_param_dict['AUIPHIRATSUF'],control_param_dict['AUIPHIRATACCEPT'],NAUINOACCEPT))    
        
    return()
    
def write_svd(NPAR):
    
    fout.write('* singular value decomposition\n')    
    wl(['SVDMODE'])
    
    fout.write('%-10i' %(NPAR))
    
    wl(['EIGTHRESH'])
    wl(['EIGWRITE'])    
    return()

def write_lsqr():
    
    fout.write('* lsqr\n')
    wl(['LSQRMODE'])
    wl(['LSQR_ATOL','LSQR_BTOL','LSQR_CONLIM','LSQR_ITNLIM'])
    wl(['LSQR_WRITE'])        
    return

def write_svda():
    
    fout.write('* svd assist\n')    
    wl(['BASEPESTFILE'])
    wl(['BASEJACFILE'])
    wl(['SVDA_MULBPA','SVDA_SCALADJ','SVDA_EXTSUPER','SVDA_SUPDERCALC'])    
    return
    
def write_sensreuse():
    
    fout.write('* sensitivity reuse\n')    
    wl(['SENRELTHRESH','SENMAXREUSE'])
    wl(['SENALLCALCINT','SENPREDWEIGHT','SENPIEXCLUDE'])    
    return

def write_param_groups(param_group_dict):
    '''Writes the parameter group information for each parameter group.  Note
    that the current version allows for multiple parameter groups but assumes
    that each parameter group uses the same parameter group values.'''
    
    fout.write('* parameter groups\n')
    for igrp in param_group_dict:
        fout.write(igrp + ' ')
        wl(['INCTYP','DERINC','DERINCLB','FORCEN','DERINCMUL','DERMTHD','SPLITTHRESH','SPLITRELDIFF','SPLITACTION'])
    return

def write_param_data(param_group_dict):
    '''Reads the parameter names for each group from a list for that group.'''
    
    fout.write('* parameter data\n')    
    for igrp in param_group_dict:
            
        init_dict = param_group_dict[igrp]['INIT']    # All parameter initialization parameters except for the initial value
        init_vals = init_dict['PARVAL1']        # A dictionary of initial parameter values         
    
        for iparam in param_group_dict[igrp]['PARAMS']:
            fout.write('%-15s' %(iparam))
            for jparam in ['PARTRANS','PARCHGLIM']:
                fout.write('%10s' %(init_dict[jparam]))
                
            fout.write('%10s' %(init_vals[iparam]))
            
            for jparam in ['PARLBND','PARUBND']:
                fout.write('%10s' %(init_dict[jparam]))            
            
            fout.write('%10s' %(igrp))
            for jparam in ['SCALE','OFFSET','DERCOM','PARTIED']:
                fout.write('%10s' %(init_dict[jparam])) 
            fout.write('\n')
    return
    
def write_obs_groups(obs_group_dict):
    '''Writes the observation group names.'''
    
    fout.write('* observation groups\n')
    for igrp in obs_group_dict:
        fout.write(igrp + '\n')
        
    return
    
def write_obs_data(obs_group_dict):
    '''Reads the observation name, data, and weight from a dataframe. Assumes
    particular column names for dataframe.'''    
    
    fout.write('* observation data\n')
            
    for igrp in obs_group_dict:
        for idx,irow in obs_group_dict[igrp]['OBS'].iterrows():
            
            iname,ival,iweight = irow['ObsName'],irow['Obs'],irow['Weight']
            fout.write('%-20s%20.8e%20.8e%10s\n' %(iname,ival,iweight,igrp))
    
    return
    
def write_derivs():
    
    fout.write('* derivatives command line\n')
    wl(['DERCOMLINE'])
    wl(['EXTDERFLE'])
    return

def write_model_command(): # IMPROVE THIS - PUT a check to require coordination with numcom
        
    fout.write('* model command line\n')
    
    for icommand in control_param_dict['COMLINE'].split(','):
        fout.write(icommand.strip() + '\n')
    
    return

def write_model_in_out(param_group_dict,obs_group_dict):
    
    fout.write('* model input/output\n')
    
    for igrp in param_group_dict:
        itpl,ipvl = param_group_dict[igrp]['TPL'],param_group_dict[igrp]['PVAL']
        fout.write('%s %s \n' %(itpl,ipvl))

    for igrp in obs_group_dict:
        iinst,isimout = obs_group_dict[igrp]['INST'],obs_group_dict[igrp]['SIMOUT']
        fout.write('%s %s \n' %(iinst,isimout))        

    return
    
def write_prior():
    
    fout.write('* prior information\n')
    return
    
def write_predict():

    fout.write('* predictive analysis\n')
    wl(['NPREDMAXMIN','PREDNOISE'])
    wl(['ABSPREDLAM','RELPREDLAM','INITSCHFAC','MULSCHFAC','NSEARCH'])
    wl(['ABSPREDSWH','RELPREDSWH'])
    wl(['NPREDNORED','ABSPREDSTP','RELPREDSTP','NPREDSTP'])
    return    

def write_reg():
    
    fout.write('* regularisation\n')
    wl(['PHIMLIM','PHIMACCEPT','FRACPHIM','MEMSAVE'])
    wl(['WFINIT','WFMIN','WFMAX','LINREG','REGCONTINUE'])
    wl(['WFFAC','WFTOL','IREGADJ','NOPTREGADJ','REGWEIGHTRAT','REGSINGTHRESH'])
    return

def get_control_params(pest_param_xls,control_sheet,param_name_col,param_val_col):
    '''Reads PEST control file parameter values from an Excel spreadsheet 
    and returns a dictionary with key=parameter name and value=parameter value.'''
    
    param_xls = pd.ExcelFile(pest_param_xls)
    param_df = param_xls.parse(control_sheet)[[param_name_col,param_val_col]]
    control_param_dict = pd.Series(param_df[param_val_col].values,index=param_df[param_name_col]).to_dict()
    
    return control_param_dict
        
# ---------------------------------------
# THIS FUNCTION CALLS THE OTHER FUNCTIONS
# ---------------------------------------

def write_control(param_group_dict,obs_group_dict,pest_control_file,pest_param_xls,control_sheet,param_name_col,param_val_col):
    '''The managing function that is called externally.'''

    NPARGP,NPAR = get_count(param_group_dict)
    NOBSGP,NOBS  = get_count(obs_group_dict)
    
    global control_param_dict 
    control_param_dict = get_control_params(pest_param_xls,control_sheet,param_name_col,param_val_col)       
    
    global fout
    with open(pest_control_file,'w') as fout:        
        
        write_control_section(NPARGP,NPAR,NOBSGP,NOBS)
        if (control_param_dict['DOAUI'] == 'aui'): write_aui(NPAR)
        if (control_param_dict['SVDMODE'] == 1): write_svd(NPAR)
        if (control_param_dict['LSQRMODE'] == 1): write_lsqr()
        if (control_param_dict['use_SVDA'] == True): write_svda()
        if (control_param_dict['DOSENREUSE'] == 'senreuse'): write_sensreuse()
        
        # These read template, instruction, and/or UCODE files   
        write_param_groups(param_group_dict)
        write_param_data(param_group_dict)
        write_obs_groups(obs_group_dict)
        write_obs_data(obs_group_dict)
        
        if (control_param_dict['use_DER'] == True): write_derivs()
            
        write_model_command()
        write_model_in_out(param_group_dict,obs_group_dict)
        
        # Make the following conditional
        if (control_param_dict['use_PRIOR'] == True): write_prior()
        if (control_param_dict['use_PREDICT'] == True): write_predict()
        if (control_param_dict['use_REG'] == True): write_reg()
            
    return
    
if __name__ == '__main__':
    write_control()