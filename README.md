This repository contains code developed to be used at the CMS Experiment, specifically for retrieving recorded hits from atomic particles on the CMS Detector within the Large Hadron Collider. 

When high-energy particles move away from the center of the detector post-collision, they leave energy signatures on silicon based micro-chips placed along the layers of the detector. Different layers of the detector are specialized for different putposes, such as the pixel tracker (one of the innermost layers) is specialized for smaller particles that are easily stopped whereas the outer layers of the detector are specialized for heavier, higher momentum particles that are harder to stop.

The detector can be thought of as a huge 'camera' that takes about 40 million 'images'/second. These images are energy maps that tell us with some error at what points we detected a particle. As you might imagine, it becomes extremely hard to tell this with millions of particles leaving their tracks every second. In order to reconstruct these tracks efficiently, we use graph neural networks. This repository handles the data preprocessing and visualization while [DeepJet](https://github.com/SwapneelM/DeepJetCore) handles the machine learning model training and bookkeeping in the CMS Software Environment (CMSSW).

The plots below show how the detector geometry is through extracting the energy signatures from the particles moving through the layers aka 'hits'. They show how the hundreds of thousands of hits (rechits, simHits, stereo rechits) tend to highlight the entire layout of the silicon plates along each layer (1 through 6) on the CMS detector. There are separate plots for inner and outer layers for clarity. 

[Outer Detector Geometry](plot-1.png)
[Inner Detector Geometry](plot-2.png)
[Not all types of Particles 'hit' all Layers](plot-3.png)

To skip some theory, we prototype how we can reconstruct a subset of tracks and a qualitative comparison of how well a random sample of tracks fits the data (rechits). Already with simple clustering in a lower dimensional space, we can see the tracks fitting the data relatively well. Note that these tracks are actually helical in structure, but are represented using a subset of their [track parameters](https://www-cdf.fnal.gov/physics/new/qcd/ue_escan/etaphi.html) (transverse momentum, pseudo-rapidity, azimuthal angle) for the clustering and plotted in 2 dimensions.

[Track Reconstruction using Clustering](scatter-1.png)

# TrackingNtuples
Custom Plugin for the CMSSW Repository

* This serves as a sample [EDAnalyzer](https://twiki.cern.ch/twiki/bin/view/Main/CMSSWatFNALANALYZER) 
to access event data from root files in the CMSSW Environment.

* It runs within the [CMS Patatrack Environment](https://github.com/cms-patatrack/cmssw/).

* For the Data Preprocessing, Plots, and Visualizations check the [Notebook](https://github.com/SwapneelM/TrackingNtuples/blob/master/TrackingNtuples/scripts/Track%20Reconstruction%20-%20Plotting%20Data.ipynb)
