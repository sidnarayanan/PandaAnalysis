#!/usr/bin/env python 

import argparse
from sys import argv
from os import getenv
parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str)
parser.add_argument('--infile',metavar='infile',type=str,default=None)
args = parser.parse_args()
if not args.infile:
    args.infile = getenv('PANDA_FITTING') + '/fittingForest.root'

argv = []
import ROOT as root
from math import sqrt
from collections import namedtuple 
from array import array
from PandaCore.Utils.load import Load
from PandaCore.Tools.root_interface import draw_hist, read_tree

Load('PandaCoreTools')
Load('PandaCoreDrawers')

# create some global variables
plot = root.HistogramDrawer()
plot.SetTDRStyle()
plot.InitLegend()
plot.SetLumi(35.8)

f_input = root.TFile(args.infile)
fztow = root.TFile(getenv('PANDA')+'/data/theory_uncs/wtoz_unc.root')
fztoa = root.TFile(getenv('PANDA')+'/data/theory_uncs/atoz_unc.root')
f_input.cd()

recoil_bins = array('f',[250,280,310,350,400,450,600,1000])
hbase = root.TH1D('dummy','',len(recoil_bins)-1,recoil_bins)
hbase.GetXaxis().SetTitle('U [Gev]')
hbase.GetYaxis().SetTitle('')

Region = namedtuple('Region',['name','label','main_proc','main_bkg','theory_uncs'])
regions = {
        'photon' : Region('photon','#gamma','Pho','QCD',False),
        'dimuon' : Region('dimuon','Z#rightarrow#mu#mu','Zll','ttbar',True),
        'dielectron' : Region('dielectron','Z#rightarrowee','Zll','ttbar',True),
        'singlemuonw' : Region('singlemuonw','W#rightarrow#mu#nu','Wlv','ttbar',False),
        'singleelectronw' : Region('singleelectronw','W#rightarrowe#nu','Wlv','ttbar',False),
        'singlemuontop' : Region('singlemuontop','t#rightarrowb#mu#nu','ttbar',None,False),
        'singleelectrontop' : Region('singleelectrontop','t#rightarrowbe#nu','ttbar',None,False),
}

ztow_uncs = {
        'wrenscale_up' : fztow.Get("znlo1_over_wnlo1_renScaleUp"),
        'wrenscale_down' : fztow.Get("znlo1_over_wnlo1_renScaleDown"),
        'wfacscale_up' : fztow.Get("znlo1_over_wnlo1_facScaleUp"),
        'wfacscale_down' : fztow.Get("znlo1_over_wnlo1_facScaleDown"),
        'wpdf_up' : fztow.Get("znlo1_over_wnlo1_pdfUp"),
        'wpdf_down' : fztow.Get("znlo1_over_wnlo1_pdfDown"),
        'wewk_up' : fztow.Get('w_ewkcorr_overz_Upcommon'),
        'wewk_down' : fztow.Get('w_ewkcorr_overz_Downcommon'),
        }
ztoa_uncs = {
        'renscale_up' : fztoa.Get("znlo1_over_anlo1_renScaleUp"),
        'renscale_down' : fztoa.Get("znlo1_over_anlo1_renScaleDown"),
        'facscale_up' : fztoa.Get("znlo1_over_anlo1_facScaleUp"),
        'facscale_down' : fztoa.Get("znlo1_over_anlo1_facScaleDown"),
        'pdf_up' : fztoa.Get("znlo1_over_anlo1_pdfUp"),
        'pdf_down' : fztoa.Get("znlo1_over_anlo1_pdfDown"),
        'ewk_up' : fztoa.Get('a_ewkcorr_overz_Upcommon'),
        'ewk_down' : fztoa.Get('a_ewkcorr_overz_Downcommon'),
        }

cuts = {
        'loose' : '0.1<top_ecf_bdt && top_ecf_bdt<0.45',
        'tight' : 'top_ecf_bdt>0.45',
        }


# helper functions
trees = {}
def get_tree(key):
    global trees, f_input
    if key not in trees:
        trees[key] = f_input.Get(key)
    return trees[key]

def build_ratio(hnum,hden):
    hratio = hnum.Clone()
    hratio.Divide(hden)
    return hratio

def subtract(hdata,hbkg):
    hsub = hdata.Clone()
    hsub.Add(hbkg,-1)
    return hsub

def draw(tree,weight,cut):
    h = hbase.Clone()
    xarr = read_tree(tree,['met',weight],cut)
    draw_hist(h,xarr,['met'],weight)
    return h

def build_unc(tree,branch_name,hist):
    if hasattr(tree,branch_name):
        return 
    ba = root.BranchAdder()
    ba.verbose = False
    ba.formula = 'genBosonPt'
    ba.newBranchName = branch_name
    ba.AddBranchFromHistogram(tree,hist)
    return

# main plotting function
def plot_ratio(num_region,den_region,cat,flat_uncs=[],shape_uncs={}):
    num = regions[num_region]
    den = regions[den_region]
    cut = cuts[cat]

    # first get the data
    hnum_data = draw(get_tree('Data_%s'%num_region),'1',cut)
    hden_data = draw(get_tree('Data_%s'%den_region),'1',cut)
    
    # now get the leading bkg
    tnum_mainproc = get_tree('%s_%s'%(num.main_proc,num_region))
    hnum_mainproc = draw(tnum_mainproc,'weight',cut)
    tden_mainproc = get_tree('%s_%s'%(den.main_proc,den_region))
    hden_mainproc = draw(tden_mainproc,'weight',cut)

    # subleading background to subtract (if necessary)
    if num.main_bkg:
        hnum_mainbkg = draw(get_tree('%s_%s'%(num.main_bkg,num_region)),'weight',cut)
        hnum_data = subtract(hnum_data,hnum_mainbkg)
    if den.main_bkg:
        hden_mainbkg = draw(get_tree('%s_%s'%(den.main_bkg,den_region)),'weight',cut)
        hden_data = subtract(hden_data,hden_mainbkg)
    
    hdata = build_ratio(hnum_data,hden_data)
    hmc = build_ratio(hnum_mainproc,hden_mainproc)

    all_uncs = flat_uncs[:]
    t_theory = None
    if num.theory_uncs:
        t_theory = tnum_mainproc
    elif den.theory_uncs:
        t_theory = tden_mainproc
    if t_theory:
        for k,v in shape_uncs.iteritems():
            build_unc(t_theory,k,v)
            hist_shift = draw(t_theory,'weight*%s'%(k),cut)
            if num.theory_uncs:
                hratio_shift = build_ratio(hist_shift,hden_mainproc)
            else:
                hratio_shift = build_ratio(hnum_mainproc,hist_shift)
            all_uncs.append(hratio_shift)
    hmc_err = hmc.Clone();
    for ib in xrange(1,hmc_err.GetNbinsX()+1):
        stat_err = hmc_err.GetBinError(ib)
        err_up = pow(stat_err,2)
        err_down = pow(stat_err,2)
        for u in all_uncs:
            if type(u)==float:
                err_up += pow(u,2)
                err_down += pow(u,2)
            else:
                diff = u.GetBinContent(ib) - hmc.GetBinContent(ib)
                if diff>0:
                    err_up += pow(diff,2)
                else:
                    err_down += pow(diff,2)
        err_up = sqrt(err_up)
        err_down = sqrt(err_down)
        hmc_err.SetBinError(ib,(err_up+err_down)/2)


                        
    root.gStyle.SetOptStat(0)

    plot.cd()
    plot.Reset()

    hmc_err.SetMinimum(0) 
    hmc_err.SetMaximum(2*max(hdata.GetMaximum(),hmc.GetMaximum()))
    hmc.SetLineWidth(2)
    hmc.SetLineColor(2)
    hmc_err.SetLineWidth(0)
    hmc_err.SetLineColor(root.kGray)
    hmc_err.SetFillColor(root.kGray)
    plot.AddCMSLabel()
    plot.AddLumiLabel(True)
    plot.AddPlotLabel('#frac{%s}{%s}'%(num.label,den.label),0.18,0.77,False,42,0.05,11)
    plot.AddHistogram(hmc_err,'Stat+sys unc',root.kExtra1,root.kGray,'e2')
    plot.AddHistogram(hmc,'Prediction',root.kExtra2,root.kRed,'el')
    plot.AddHistogram(hdata,'Data',root.kData,root.kBlack,'elp')

    plotname = 'ratio_%s_%s_%s'%(cat,num_region,den_region)
    plot.Draw(args.outdir+'/',plotname)

plot_ratio('photon','dimuon','tight',flat_uncs=[0.02,0.01],shape_uncs=ztoa_uncs)
plot_ratio('photon','dimuon','loose',flat_uncs=[0.02,0.01],shape_uncs=ztoa_uncs)

plot_ratio('photon','dielectron','tight',flat_uncs=[0.02,0.01],shape_uncs=ztoa_uncs)
plot_ratio('photon','dielectron','loose',flat_uncs=[0.02,0.01],shape_uncs=ztoa_uncs)

plot_ratio('singlemuonw','dimuon','tight',flat_uncs=[0.02,0.01],shape_uncs=ztow_uncs)
plot_ratio('singlemuonw','dimuon','loose',flat_uncs=[0.02,0.01],shape_uncs=ztow_uncs)

plot_ratio('singleelectronw','dielectron','tight',flat_uncs=[0.02,0.01],shape_uncs=ztow_uncs)
plot_ratio('singleelectronw','dielectron','loose',flat_uncs=[0.02,0.01],shape_uncs=ztow_uncs)

plot_ratio('dimuon','dielectron','tight',flat_uncs=[0.02,0.04])
plot_ratio('dimuon','dielectron','loose',flat_uncs=[0.02,0.04])

plot_ratio('singlemuonw','singleelectronw','tight',flat_uncs=[0.01,0.02])
plot_ratio('singlemuonw','singleelectronw','loose',flat_uncs=[0.01,0.02])


fztoa.Close()
fztow.Close()

