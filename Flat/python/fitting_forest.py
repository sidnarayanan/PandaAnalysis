'''PandaAnalysis.Flat.fitting_forest

Converts large flat trees into smaller trees used for fitting
''' 

import PandaCore.Tools.root_interface as root_interface 
from PandaCore.Tools.Misc import tAND 
from PandaCore.Utils.logging import logger
import numpy as np 
import ROOT as root 
from itertools import chain 

_trees = {} # used to prevent gc and multiple openings of a file
_files = {}

treename = 'events' # input tree name

class Process:
    def __init__(self, name, tree, cut, variables, weights):
        '''        
        Arguments:
            name {str} -- name of this process
            tree {ROOT.TTree} -- input TTree
            cut {str} -- selection to apply
            variables {dict} -- map : output branch name -> input formula
            weights {dict} -- map : systematic shift name -> input formula
        '''
        self.name = name 
        self.tree = tree
        self.cut = cut 
        self.variables = variables
        self.__outputs = {}
        if not weights:
            self.weights = {}
            self.nominal_weight = '1'
        elif type(weights)==str:
            self.nominal_weight = weights 
            self.weights = {}
        else:
            self.nominal_weight = weights['nominal']
            self.weights = weights.copy()
            del self.weights['nominal']
        self.all_branches = self.variables.copy()
        self.all_branches.update(self.weights)
        self.all_branches['nominal'] = self.nominal_weight
    def __write_out(self, f_out, xarr, fields, postfix):
        repl = {self.all_branches[f] : (f if f in self.variables else 'weight') for f in fields}
        varr = xarr[list(set(repl.keys()))]
        root_interface.rename_dtypes(varr, repl)
        root_interface.array_as_tree(xarr = varr, 
                                     treename = self.name+postfix, 
                                     fcontext = f_out)
    def run(self, f_out):
        logger.info('fitting_forest.Process.run', 'Running '+self.name)
        branches = sorted(self.all_branches.values())
        try:
            xarr = root_interface.read_tree(tree = self.tree, 
                                            branches = branches, 
                                            cut = self.cut)
            fields = self.variables.keys()+['nominal']
            self.__write_out(f_out, xarr, fields, '')
            for shift, weight in self.weights.iteritems():
                fields = self.variables.keys()+[shift]
                self.__write_out(f_out, xarr, fields, '_'+shift)
        except ValueError as e:
            logger.error('fitting_forest.Process.run', str(e))
            return


class RegionFactory:
    def __init__(self, name, cut, variables, mc_variables, mc_weights):
        '''        
        Arguments:
            name {str} -- name of this region
            cut {str} -- selection to apply
            variables {dict} -- map : output branch name -> input formula
            mc_variables {dict} -- map : output branch name -> input formula (gen info)
            mc_weights {dict} -- map : systematic shift name -> input formula
        '''
        self.name = name 
        self.cut = cut 
        self.variables = variables
        self.mc_variables = mc_variables
        self.mc_weights = mc_weights
        self.__mc_procs = [] # procs that get mc_weights
        self.__data_procs = [] # procs that don't get mc_weights
    def add_process(self, input_info, pname, is_data=False, extra_weights=None, extra_cut=None):
        global _trees,  _files,  treename 
        if type(input_info)==str: # assume it's a filepath
            if input_info in _trees:
                tree = _trees[input_info]
            else:
                _files[input_info] = root.TFile.Open(input_info)
                tree = _files[input_info].FindObjectAny(treename)
                _trees[input_info] = tree 
        else: # assume it's a TTree
            tree = input_info
        weights = None 
        if is_data:
            if extra_weights:
                weights = extra_weights
        else:
            weights = self.mc_weights.copy()
            if extra_weights:
                weights.update(extra_weights)
        cut = tAND(self.cut, extra_cut) if extra_cut else self.cut 
        variables_ = self.variables.copy()
        if not is_data:
            variables_.update(self.mc_variables)
        pname += '_' + self.name
        proc = Process(pname, tree, cut, variables_, weights)
        if is_data:
            self.__data_procs.append(proc)
        else:
            self.__mc_procs.append(proc)
    def run(self, f_out_path):
        f_out = root.TFile.Open(f_out_path, 'RECREATE')
        for proc in chain(self.__data_procs, self.__mc_procs):
            proc.run(f_out)
        f_out.Close() 
        logger.info('fitting_forest.RegionFactory.run', 'Created output in %s'%f_out_path)

