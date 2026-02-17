# DielectronTRDReconstruction
For dielectron final state reconstruction in GlueX-II Spring 2025 using the large-scale GEM-TRD prototype

################################  
### Original Workflow from Lubomir  
On JLab ifarm, set up gluex env:
```
cd /work/halld2/home/$USERNAME/dsoft/halld_recon/  
source setup_gluex.csh  
source setenv.csh   
   
```
In working directory `run/`, use the config file that points to the Reaction Filter plugin (`jana.config`) to process the data files. The skimmed data files with e+e- FS in all acceptance live in the `data/` directory:
``` 
cd run/  
hd_root --loadconfigs jana.config /work/halld2/home/lkasper/dsoft/halld_recon/data/skim_phi/jpsi_skims/*.evio -POUTPUT_FILENAME=hd_root_JPsi_default.root  
   
```
This will output 2 files: `hd_root_JPsi_default.root` and `tree_epem__B3.root`.  
Open `tree_epem__B3.root`, which contains 4 vector information for the events stored after using the Reaction Filter. Run the first TSelector (`jpsi_lp_selector.C`) over this, which takes the 4 vectors from the Reaction Filter output and applies the KinFit, boosts in different coord. system, i.e. transforms the variables into something more useful/accessible.
``` 
root -l tree_epem__B3.root  
epem__B3_Tree->Process(“jpsi_lp_selector.C++”)  

```
This will output another root file, `jpsi_lp_out2.root`, which contains a new TTree `JP` with these transformed variables, etc.  
To see, for example, the inv. mass plot for the J/Psi, do:
``` 
root -l jpsi_lp_out2.root  
JP->Draw(“Minv”)  

```
At this point, no cuts have been applied. A second TSelector, `JP_selector.C` is used for E/p, chi^2 cut on kinematic fit, etc. Inside this file is also marked with a comment for `important cuts for TRD`.  
``` 
JP->Process(“JP_selector.C++”)

```
This will output a file named `s3de6.root`.  
To see the inv. mass plot for the J/Psi with these selections now applied, do:  
``` 
root -l s3de6.root  
dminva->Draw()  

```
To use a macro with RooFit to fit this distribution, first copy the root file and then execute the macro `simfit_de_AN_ok.C`:  
``` 
cp s3de6.root s3de6_roo.root  
root -l simfit_de_AN_ok.C
  
``` 
This will output the pdf `TRD_JPsi_default.pdf`.  
Another macro, `GetJpsis_in_TRD.C`, can be used to visualize the e+e- pairs in all acceptance at the Z plane of the TRD:  
```
root -l GetJpsis_in_TRD.C  

``` 
