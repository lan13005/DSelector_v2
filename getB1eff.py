import root_pandas as rp
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import os
import ROOT
from decimal import Decimal
import pandas as pd

import mplhep as hep
mpl.rcParams.update(mpl.rcParamsDefault)
plt.style.use([hep.styles.ATLAS])
SMALL_SIZE = 20
MEDIUM_SIZE = 22
BIGGER_SIZE = 24
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=17)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


target=1.22*10**(-9) #nanobarns

fig,ax=plt.subplots(2,3,figsize=(18,9))
ax=ax.flatten()
binedges=np.linspace(6.6,11.4,21)

xerrs=(binedges[1]-binedges[0])/2
bincenters=(binedges[1:]+binedges[:-1])/2
print("bin edges: {}".format(binedges))
print("bin centers: {}".format(bincenters))
binMin=6
binMax=10
if binMax==-1:
    includingBinMax=len(bincenters)
else:
    includingBinMax=binMax#+1

emin=binedges[binMin]
emax=binedges[binMax]

getCS=True # should we get the array of cross sections
getFlux=True # should we get the array of fluxes 
getEff=True # should we get the array of efficiencies
getYields=True # should we get the array of yields

rerunFlux=False # should we get the flux again

outputFolder="b1_efficiency/"
weightBranch="weightASBS" # "weightASBS", "AccWeight"
os.system("mkdir -p "+outputFolder)

####################################
# Obtain Cross Section of b1
####################################
if getCS:
    print("\n\nGetting cross sections\n------------------------")
    crossSections=[1.37,1.37,1.35,1.36,1.35,1.3,1.26,1.18,1.1,1.1,1.03,1.04,1,1,0.95,0.95,0.9,0.9,0.82,0.82]
    crossSections=np.array(crossSections)*1000 # Foda has it in microbarns
    crossSectionErrs=[0.03,0.03,0.03,0.03,0.02,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01]
    crossSectionErrs=np.array(crossSectionErrs)
    assert((len(binedges)-1==len(crossSections)) and (len(crossSections)==len(crossSectionErrs)))
    ax[4].errorbar(bincenters,crossSections,yerr=crossSectionErrs,xerr=xerrs,fmt="ro",ecolor="black")
    ax[4].set_ylabel("b1 Cross Section (nb)",size=16)
    ax[4].set_xlabel("Beam Energy (GeV)",size=16)
    ax[4].set_ylim(bottom=0)
    ax[4].axvline(emin,c="black",linestyle="--")
    ax[4].axvline(emax,c="black",linestyle="--")
    print("Got it!")

####################################
# Obtain Tagged Flux
####################################
if getFlux:
    print("\n\nGetting tagged flux\n------------------------")
    fluxCounts=[]
   
    runs=["2017","2018_1","2018_8"]
    runStarts=[30274,40856,50677]
    runEnds=[31057,42577,51768]
    #rcdbQueries=[""," --rcdb-query='@is_production'"," --rcdb-query='@is_2018production and @status_approved and beam_on_current > 49'"]
    rcdbQueries=["","",""]
    rcdbQueries=[""," --rcdb-query='@is_2018production and @status_approved'"," --rcdb-query='@is_2018production and @status_approved and beam_on_current > 49'"]
    for i in range(3):
        fluxCounts.append([])
        cmd_base="/d/grid13/gluex/gluex_top/hd_utilities/hd_utilities-1.17/psflux/plot_flux_ccdb.py --begin-run="+str(runStarts[i])+" --end-run="+str(runEnds[i])
        cmd_bins="--num-bins="+str(len(binedges)-1)
        cmd_lowE="--energy-min="+str(binedges[0])
        cmd_uppE="--energy-max="+str(binedges[-1])
        cmds=[cmd_base,cmd_bins,cmd_lowE,cmd_uppE]
        cmd=" ".join(cmds)
        cmd+=rcdbQueries[i]
        print("Running following command:")
        print(cmd)
        if rerunFlux:
            os.system(cmd)
            os.system("mv flux_"+str(runStarts[i])+"_"+str(runEnds[i])+".root "+outputFolder)
    
        fluxFile=ROOT.TFile.Open(outputFolder+"flux_"+str(runStarts[i])+"_"+str(runEnds[i])+".root")
        fluxHist=fluxFile.Get("tagged_flux")
        for j in range(fluxHist.GetNbinsX()):
            count=fluxHist.GetBinContent(j+1)
            print("Bin{0} counts: {1}".format(j,count))
            fluxCounts[i].append(count)
    fluxCounts=np.array(fluxCounts)
    fluxCounts=fluxCounts.sum(axis=0)

    fluxErrors=np.sqrt(fluxCounts)
    ax[3].errorbar(bincenters,fluxCounts,yerr=fluxErrors,xerr=xerrs,fmt="ro",ecolor='black')
    ax[3].set_ylabel("Tagged Flux",size=16)
    ax[3].set_xlabel("Beam Energy (GeV)",size=16)
    ax[3].set_ylim(bottom=0)
    ax[3].axvline(emin,c="black",linestyle="--")
    ax[3].axvline(emax,c="black",linestyle="--")

####################################
# Get reconstruction efficiency of b1 to etapi
####################################
if getEff:
    print("\n\nGetting efficiencies\n------------------------")
    thrownEnergy="ThrownBeam__GeneratedEnergy"
    dataEnergy="Ebeam_thrown" # "Ebeam" if we want the reconstructed
    baseFolder="/d/grid17/ln16/myDSelector/"
    columns=[dataEnergy,"AccWeight","weightASBS","insideEllipse","Mpi0eta"]
    recon=rp.read_root(baseFolder+"degALL_b1vps_as_4g_mEllipse_8288_tLT1_chi13_treeFlat_DSelector.root",columns=columns)
    thrown=rp.read_root(baseFolder+"degALL_b1vps_as_4g_gen_trees_DSelector.root",columns=[thrownEnergy])
    #recon=recon[recon.insideEllipse]
    
    dat_counts, edges, _ = ax[0].hist(recon[dataEnergy],bins=binedges,weights=recon[weightBranch],label="8.2<\n$E_{beam}$\n<8.8 GeV")
    thrown_counts, edges, _ = ax[1].hist(thrown[thrownEnergy],bins=binedges,label="8.2<\n$E_{beam}$\n<8.8 GeV")
    ax[0].legend()
    ax[1].legend()
    ax[0].set_ylabel("Reconstructed Yield",size=16)
    ax[1].set_ylabel("Thrown Yield",size=16)
    ax[0].set_xlabel("Beam Energy (GeV)",size=16)
    ax[1].set_xlabel("Beam Energy (GeV)",size=16)
    #ax[0].set_ylim(bottom=0)
    ax[1].set_ylim(bottom=0)
    ax[0].axvline(emin,c="black",linestyle="--")
    ax[1].axvline(emin,c="black",linestyle="--")
    ax[0].axvline(emax,c="black",linestyle="--")
    ax[1].axvline(emax,c="black",linestyle="--")
    
    
    # We will skip the bins that have exactly 0 entries, cant really calculate an efficiency and error from them due to div-by-zero
    skipSinceZero=[False if dat_count==0 and thrown_count==0 else True for dat_count,thrown_count in zip(dat_counts,thrown_counts)]
    dat_counts=dat_counts[skipSinceZero]
    thrown_counts=thrown_counts[skipSinceZero]
    print(dat_counts)
    print(thrown_counts)
    _bincenters=bincenters[skipSinceZero]
    efficiencies=[dat_count/thrown_count for dat_count,thrown_count in zip(dat_counts,thrown_counts)]
    
    # Calculate efficiencies
    dat_counts_err=np.sqrt(abs(dat_counts)) # since we oversubtract sometimes we will have negative yields
    thrown_counts_err=np.sqrt(thrown_counts)
    efficiencies_err=efficiencies*np.sqrt( (dat_counts_err/dat_counts)*(dat_counts_err/dat_counts) + (thrown_counts_err/thrown_counts)*(thrown_counts_err/thrown_counts) )
    
    print("efficiency:")
    print(efficiencies)
    print("efficiencies error:")
    print(efficiencies_err)
    
    ax[2].errorbar(_bincenters,efficiencies,yerr=efficiencies_err,xerr=xerrs,fmt="ro",ecolor='black')
    ax[2].set_ylabel("Efficiency",size=16)
    ax[2].set_xlabel("Beam Energy (GeV)",size=16)
    ax[2].ticklabel_format(style='sci')
    #ax[2].set_ylim(bottom=0)
    ax[2].axvline(emin,c="black",linestyle="--")
    ax[2].axvline(emax,c="black",linestyle="--")


####################################
# Calculating Yield
####################################
if getYields:
    print("\n\nGetting Expected Yields\n------------------------")
    cs=crossSections[binMin:includingBinMax]
    cserr=crossSectionErrs[binMin:includingBinMax]
    fc=fluxCounts[binMin:includingBinMax]
    fcerr=fluxErrors[binMin:includingBinMax]
    eff=efficiencies
    efferr=efficiencies_err
    energies=bincenters[binMin:includingBinMax]
    
    expectedYields=cs*eff*fc*target
    expectedYieldErrs=expectedYields*np.sqrt((cserr/cs)*(cserr/cs)+(fcerr/fc)*(fcerr/fc)+(efferr/eff)*(efferr/eff))
    
    ax[5].errorbar(energies,expectedYields,yerr=expectedYieldErrs,xerr=xerrs,fmt="ro",ecolor='black')
    ax[5].set_ylabel("Expected Yield",size=16)
    ax[5].set_xlabel("Beam Energy (GeV)",size=16)
    #ax[5].set_ylim(bottom=0)

    nsigs=3
    lower3SigEstimate='%.2E' % Decimal(sum(expectedYields-nsigs*expectedYieldErrs))
    upper3SigEstimate='%.2E' % Decimal(sum(expectedYields+nsigs*expectedYieldErrs))
    # might swap upper and lower estimates later but maximumSeparation is always the upper. The most we will subtract or add is determined from upper
    maximumSeparation=sum(expectedYields+nsigs*expectedYieldErrs) 
    if sum(expectedYields+nsigs*expectedYieldErrs)<sum(expectedYields-nsigs*expectedYieldErrs): # if upper limit < lower limit
        lower3SigEstimate,upper3SigEstimate=upper3SigEstimate,lower3SigEstimate
    ax[5].set_title(r"3$\sigma$ Integral: ["+lower3SigEstimate+", "+upper3SigEstimate+"]",size=16,fontweight='bold') 

    #integral='%.2E' % Decimal(sum(expectedYields))
    #ax[5].set_title("Integral: "+integral,size=16,fontweight='bold')


###################################
# Plotting
###################################
print("\n\nPlotting final plot!\n------------------------")
plt.tight_layout()
plt.savefig(outputFolder+"efficiencyB1_"+weightBranch+".png")

fig,ax=plt.subplots(1,1,figsize=(10,8))
dataBaseFolder="/d/grid17/ln16/myDSelector/amptools/zPhase1_t0103061_e79828890/baseFiles_v2/"
datas=[]
for run in ["2017","2018_1","2018_8"]:
    datas.append(rp.read_root(dataBaseFolder+"degALL_data_"+run+"_mEllipse_8288_tLT1_treeFlat_DSelector.root",columns=["Mpi0eta","AccWeight","weightASBS"]))
datas=pd.concat(datas)
edges = np.histogram(datas["Mpi0eta"],weights=datas[weightBranch],bins=100)[1]
ax.hist(datas["Mpi0eta"],weights=datas[weightBranch],bins=edges,histtype='step',color='black',linewidth=2,label="Phase1 Data")
ax.hist(recon["Mpi0eta"],weights=recon[weightBranch]*maximumSeparation/recon[weightBranch].sum(),bins=edges,histtype='step',
        color='red',linewidth=2,label=r"Upper 3$\sigma$ limit b1 leakage")
ax.set_xlabel(r"$M(4\gamma)$ GeV")
ax.set_ylabel("Entries / {0:0.3f} GeV".format(edges[1]-edges[0]))
ax.legend()
plt.savefig(outputFolder+"expectedYield_"+weightBranch+".png")














