# -*- coding: utf-8 -*-
import numpy as np

from curp.twobody import lib_fmm


class FMMCalculatorBase:
    
    def __init__(self, setting, natom, charges, interact_table, gname_iatoms_pairs, gpair_table):
        self.__n_crit = setting.curp.coulomb_fmm_cell_contains
        self.__theta = setting.curp.coulomb_fmm_direct_parm
        self.__gname_iatoms_pairs = gname_iatoms_pairs
        self.__gpair_table = gpair_table
        self.__gnames = []
        
        # for debug
        # from curp.twobody import test
        # obj = test.test_class()
        # x = obj.add_struct(9, 4)
        # if x[0] == 13:
        #     raise ImportError("x=5. Please check your installation.")
        
        self.__mod_fmm = lib_fmm.cal_fmm()
        self.__mod_fmm.setup(int(natom), int(self.__n_crit), float(self.__theta),
                            charges, interact_table, self.__gname_iatoms_pairs, self.__gpair_table)
         
    def initialize(self, crd):
        self.__crd = crd
        self.__mod_fmm.initialize(crd)
        
    def get_n_crit(self):
        return self.__n_crit
    
    def get_theta(self):
        return self.__theta
    
    def get_gname_iatoms_pairs(self):
        return self.__gname_iatoms_pairs
    
    def get_gpair_table(self):
        return self.__gpair_table
    
    def get_mod_fmm(self):
        return self.__mod_fmm
    
    def get_crd(self):
        return self.__crd
    
####################################################################################################################   

class FMMCellMaker(FMMCalculatorBase):
    
    def __init__(self, setting, natom, charges, interact_table, gname_iatoms_pairs, gpair_table):
        FMMCalculatorBase.__init__(self, setting, natom, charges, interact_table, gname_iatoms_pairs, gpair_table)
    
    def make_cells(self):
        
        # build tree
        all_cells = self.build_all_tree(self.get_gname_iatoms_pairs(), self.get_crd(), self.get_n_crit())
        
        return all_cells
    
    def build_all_tree(self, group_atoms, crd, n_crit):
        return self.get_mod_fmm().setup_all_cells()

###################################### Calculator ##########################################################################

class FMMCellCalculator(FMMCellMaker):

    def __init__(self, setting, natom, charges, interact_table, gname_iatoms_pairs, gpair_table):
        FMMCellMaker.__init__(self, setting, natom, charges, interact_table, gname_iatoms_pairs, gpair_table)
        self.__mod_fmm = self.get_mod_fmm()
        
    def cal_fmm(self, all_cells):
        
        self.__mod_fmm.get_all_cells(all_cells)
        
        self.__mod_fmm.cal_force_fmm()
        
        # # get multipole arrays
        # self.__mod_fmm.cal_p()
        
        # # upward sweep
        # self.__mod_fmm.cal_M2M()

        # # evaluate potential
        # self.__mod_fmm.cal_force()
        
        atomwise = dict(atomwise_i = self.__mod_fmm.atomwise_i,
                        atomwise_j = self.__mod_fmm.atomwise_j, 
                        atomwise_f = self.__mod_fmm.atomwise_f,
                        atomwise_r = self.__mod_fmm.atomwise_r)
        cellwise = dict(cellwise_i = self.__mod_fmm.cellwise_i,
                        cellwise_J = self.__mod_fmm.cellwise_J, 
                        cellwise_f = self.__mod_fmm.cellwise_f,
                        cellwise_r = self.__mod_fmm.cellwise_r)
        # raise ValueError('type atomwise_i:{}, atomwise_j:{}, atomwise_f:{}, atomwise_r:{}\n'
        #                  'type cellwise_i:{}, cellwise_J:{}, cellwise_f:{}, cellwise_r:{}'.format(
        #     type(atomwise['atomwise_i']), type(atomwise['atomwise_j']),
        #     type(atomwise['atomwise_f']), type(atomwise['atomwise_r']),
        #     type(cellwise['cellwise_i']), type(cellwise['cellwise_J']),
        #     type(cellwise['cellwise_f']), type(cellwise['cellwise_r'])))
        return dict(atomwise=atomwise, cellwise=cellwise)
        

####################################################################################################################

def check_setting(setting):
    """Check setting parameters for FMM calculation."""

    if setting.curp.method != 'energy-flux' or 'heat-flux':
        raise ValueError('The method should be "energy-flux" or "heat-flux"')
    else:
        if setting.curp.flux_grain != 'group':
            raise ValueError('The flux grain should be "group"')
        else:
            return True

