# Gravity Machine data analysis 

## What it is
This repo contains python and jupyter notebook scripts used for [Gravity Machine](https://gravitymachine.org) data analysis. 

To first visualize Gravity Machine data use [this repository](https://github.com/deepakkrishnamurthy/gravitymachine-analysis-gui).

The repo implements:
- Simple plotting and analysis functionality for 3D trajectories (like those from [Gravity Machine data](https://gravitymachine.org/gallery)).
- Mean-Squared-Displacement calculations.
- [Weighted-Least-Squares (Including correlation in error)](https://github.com/impaktor/wlsice.git) fitting for simple model building from trajectories.
- Calculating and plotting of trajectory statistics. 
- Particle-Image-Velocimetry (using [OpenPIV](https://github.com/OpenPIV/openpiv-python.git)) of [Gravity Machine image data](https://gravitymachine.org/gallery). Including using measured fluid velocities to compute true displacements relatuve to fluid. 


## Getting started
To install the software dependencies follow the instructions in [installation.md](https://github.com/deepakkrishnamurthy/gravitymachine-analysis-scripts/blob/2ace1f3f49892c10a9a2dd20e4d2d5d86999c447/installation.md)

## Selected Publications
1. Krishnamurthy, Deepak, Hongquan Li, François Benoit du Rey, Pierre Cambournac, Adam G. Larson, Ethan Li, and Manu Prakash. "Scale-free vertical tracking microscopy." Nature Methods (2020): 1-12. [Weblink](https://www.nature.com/articles/s41592-020-0924-7)
2. Krishnamurthy, Deepak, Hongquan Li, François Benoit du Rey, Pierre Cambournac, Adam Larson, and Manu Prakash. "Scale-free Vertical Tracking Microscopy: Towards Bridging Scales in Biological Oceanography." bioRxiv (2019): 610246. [Weblink](https://doi.org/10.1101/610246)

## To cite this work
	@article{krishnamurthy2020scale,
	  title={Scale-free vertical tracking microscopy},
	  author={Krishnamurthy, Deepak and Li, Hongquan and du Rey, Fran{\c{c}}ois Benoit and Cambournac, Pierre and Larson, Adam G and Li, Ethan and Prakash, Manu},
	  journal={Nature Methods},
	  volume={17},
	  number={10},
	  pages={1040--1051},
	  year={2020},
	  publisher={Nature Publishing Group}
	}
