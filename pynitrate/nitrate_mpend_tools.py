# -*- coding: utf-8 -*-
"""
Functions for operating on nitrate solute concentrations
simulated with MODPATH.

Created on Fri Sep 19 16:10:39 2014
@author: Wesley Zell
"""

# --- STOP PARAM SET ----

def partition_nitrate(row):
    """Partitions available nitrate to either (1) leachate to water table for time
    step k or (2) residual nitrate in soil column for time step k+1.
    
    Rules, as a function of mean fall (Oct-Dec) precipitation P_bar:
    if P_bar > flood_threshold, all nitrate to leachate
    if P_bar < drought_threshold, all nitrate remains in soil column
    otherwise, Fx(P_bar) to leachate, 1-Fx(P_bar) remains in soil column.
    
    NOTE: this function is called by df.apply()
    
    returns: """
    
    if   (row.Fall_Precip_Percentile > flood_threshold):
        to_leach,to_soil = 1.,0.
    elif (row.Fall_Precip_Percentile < drought_threshold):
        to_leach,to_soil = 0.,1.
    else:
        to_leach = row.Fall_Precip_Percentile
        to_soil  = 1 - to_leach
        
    return to_soil,to_leach

def get_uptake_fxn(base_uptake,improve_uptake,year_improve):
    '''Writes the uptake function.  Constant until improve year, then
    linear improvement.'''
    
    uptake_dict = dict()
    
    for iyear in loading_ts_dict.keys():
        
        if iyear <= year_improve:
            iuptake = base_uptake
        else:
            iuptake = base_uptake + improve_uptake * (iyear - year_improve)
            
        uptake_dict[iyear] = iuptake
        
    return uptake_dict