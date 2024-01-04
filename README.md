# Python tools for OpenFOAM simulations of filament shapes in embedded 3D printing

## Authors
- Leanne M. Friedrich
    - National Institute of Standards and Technology, MML
    - Leanne.Friedrich@nist.gov
    - https://github.com/leanfried
    - ORCID: 0000-0002-0382-3980
- Ross T. Gunther
    - Georgia Institute of Technology
    - https://github.com/RossGunther
    - ORCID: 0000-0002-6442-5396

## Contact
- Leanne Friedrich
    - Leanne.Friedrich@nist.gov
    - https://github.com/leanfried

## Description

In embedded 3D printing, a nozzle is embedded into a support bath and extrudes filaments or droplets into the bath. Using OpenFOAM, we simulated the extrusion of filaments and droplets into a moving bath, for single filaments, single filaments being disturbed by a nozzle, and printing pairs of filaments. OpenFOAM is an open source computational fluid dynamics solver. This repository contains the following Python tools:

- Tools for generating input files for OpenFOAM v1912 or OpenFOAM v8 tailored to a nozzle extruding a filament into a static support bath.
- Tools for monitoring the status of OpenFOAM simulations and aborting them if they are too slow.
- Tools for moving output files between storage locations. (For example, it can automatically move all files to a server, but only necessary files to your hard drive)
- Tools for generating images and tables from the 3D time series.
- Tools for compiling images into videos.
- Tools for analyzing, summarizing, and plotting data.

## Disclaimer

Certain equipment, instruments, software, or materials are identified in this code in order to specify the experimental procedure adequately.  Such identification is not intended to imply recommendation or endorsement of any product or service by NIST, nor is it intended to imply that the materials or equipment identified are necessarily the best available for the purpose.


## Change log

|version|Timeframe|Scope|
|-------|---------|-----|
|1.0.0  |April 2020 - April 2021|Single filaments, cylindrical nozzle|
|1.1.0  |April 2021 - April 2022|Single filaments, cylindrical or conical nozzle of varying size|
|1.2.0  |April 2022 - December 2023|Pairs of filaments and running the nozzle next to a filament|


## Data Use Notes


This code is publicly available according to the NIST statements of copyright,
fair use and licensing; see 
https://www.nist.gov/director/copyright-fair-use-and-licensing-statements-srd-data-and-software

You may cite the use of this code as follows:
> Friedrich, L.M. & Gunther, R.T. (2023), OpenFOAM simulations of filament shapes in embedded 3D printing, Version 1.2.0, National Institute of Standards and Technology, doi:[10.18434/mds2-3128] (Accessed XXXX-XX-XX)

The OpenFOAM input file tools are compatible with OpenFOAM v1912 and OpenFOAM v8. Between the two versions, OpenFOAM did change the syntax on a couple of input variables, so you may need to update the initialize files for future versions of OpenFOAM. Between the two versions, OpenFOAM also modified the way they output paraview files. In v1912, OpenFOAM generated a .vtm file and folder for each time step, and compiled all of the timesteps into a .vtm.series file that could easily be imported into Paraview. In v8, OpenFOAM generates a .vtk file for each time step but does not generate a .series file. file_handling.py has a function generateVTKSeries that generates these .series files, which make it easy to load the whole time series in Paraview.


## References


This code is described in the following papers:
> Friedrich, L.M., & Seppala, J.E. (2021) Simulated filament shapes in embedded 3D printing, Soft Matter, doi: 10.1039/D1SM00731A

> Friedrich, L.M., Gunther, R.T. & Seppala, J.E. (2022) Simulated stress mitigation strategies in embedded 3D bioprinting, Physics of Fluids, doi: 10.1063/5.0102573

> Friedrich, L.M., & Gunther, R.T. (2023) Simulated inter-filament fusion in embedded 3D printing, submitted for publication.

The datasets generated by this code are stored at:
> Friedrich, L.M., & Seppala, J.E. (2021) Simulated filament shapes in embedded 3D printing, Version 1.0.0, National Institute of Standards and Technology, doi: 10.18434/mds2-2391

> Friedrich, L.M., Gunther, R.T. & Seppala, J.E. (2022) OpenFOAM simulations of stress mitigation strategies in embedded 3D bioprinting, Version 1.0.0, National Institute of Standards and Technology, doi: 10.18434/mds2-2604

> Friedrich, L.M. & Gunther, R.T. (2023) OpenFOAM simulations of filament fusion in embedded 3D printing, Version 1.0.0, National Institute of Standards and Technology, doi:10.18434/mds2-3129


## Usage

For full functionality, you will need to make the following files:

- Copy `configs/config_template.yml` and call it `configs/config.yml`. Change the path names to the paths to your data.
- To use the paraview scripts, you will need to establish a virtual environment:

  ```
  python -m virtualenv ./env
  .\env\scripts\activate
  python3 -m pip install -r requirements.txt || conda install --file requirements.txt
  deactivate
  ```

This repository is linked with the paper and dataset linked at the top of this document. Typical workflow for a Windows desktop for a single parent folder (e.g. HBHBsweep) was as follows:

1.  
    1. In JupyterLab, generate input files using `noz3dscript.ipynb`. For each simulation, this generates `legend.csv`, `geometry.csv`, `case`, `0`, `constant`, `system`, bash files, slurm files. For the parent folder, this generates mesh and geometry files. 

    2. (optional file generation) Move files to cluster, other computers, etc.

2.  Run simulations. To add these to our slurm queue, we use sbatch run.slurm.

3.  Move complete simulations back to the local server.

4.  Generate tables of points and images using paraview. Modify `comboscript.py` to include relevant folders and desired image types. In Command Prompt, run (`[path]\pvpython.exe [path]\comboscript.py`).

5.  
    1. In JupyterLab, use `interfacemetrics_XXX.ipynb` to generate `sliceSummaries.csv`, `steadyPositions.csv`, and `steadyTimes.csv` for each folder.
    
    2. In JupyterLab, use `interfacemetrics_XXX.ipynb` to make plots.

## Data Overview


This data follows the basic storage structure established by OpenFOAM.

The files included in this publication use the following hierarchy:

- *README.md*

- *LICENSE*

- *requirements.txt*  
    List of required packages, for use with virtual environments. This file can be used with virtualenv for easy import of dependencies.

- **configs/**
    - *config.yml*, *config_template.yml*  
        establish paths to data and tolerable timing for simulations. Copy *config_template.yml* to *config.yml* and edit file locations in *config.yml* before proceeding.

    - *logging.yml*  
        set up log generation
        
- **jupyter_notebooks/**  
    - *convergence.ipynb* (no longer compatible)   
         Jupyter notebook for tracking convergence over the solve. 
    
    - *interfacemetrics_adjacent_vicsweep.ipynb* (no longer compatible)  
        Jupyter notebook for analyzing preliminary results for the adjacent filament dataset (Friedrich, L.M., & Gunther, R.T. (2023) Simulated inter-filament fusion in embedded 3D printing, submitted for publication.)   
    
    - *interfacemetrics_conicalNozzle.ipynb* (no longer compatible)  
        Jupyter notebook for analyzing OpenFOAM simulation data, for the stress mitigation dataset (Friedrich, L.M., Gunther, R.T. & Seppala, J.E. (2022) Simulated stress mitigation strategies in embedded 3D bioprinting, Physics of Fluids)
        
    - *interfacemetrics_fusion.ipynb*  
        Jupyter notebook for analyzing results for the adjacent filament dataset (Friedrich, L.M., & Gunther, R.T. (2023) Simulated inter-filament fusion in embedded 3D printing, submitted for publication.)
        
    - *interfacemetrics_LapRD.ipynb* (no longer compatible)  
        Jupyter notebook for analyzing OpenFOAM simulation data, for the dataset that investigates models fitted to Laponite gels (unpublished)

    - *interfacemetrics_viscsweep.ipynb*  (no longer compatible)   
        Jupyter notebook for analyzing OpenFOAM simulation data, for the viscosity sweep dataset (Friedrich, L.M., & Seppala, J.E. (2021) Simulated filament shapes in embedded 3D printing, Soft Matter, doi: 10.1039/D1SM00731A)
        
    - *interfacemetrics_yielding.ipynb*  (no longer compatible)  
        Jupyter notebook for analyzing OpenFOAM simulation data, for the dataset that investigates yielded behaviors (unpublished)
        
    - *noz3dscript.ipynb*    
        Jupyter notebook for generating OpenFOAM input files.
        
    - *point_sort.ipynb* 
        Jupyter notebook for sorting cross-section points into discrete objects
        
    - *videoCombine.ipynb*    
        Jupyter notebook for combining all of the time series for all simulations into one big video.

- **logs/**  
    for holding logs. This folder is created automatically when ```exportLog=True``` is added as an argument to ```fp.openLog(...)``` in any Jupyter notebook.

- **paraviewscripts/**  
    pvpython scripts for generating images and tables. These run on pvpython.exe. Some scripts use packages that are not native to pvpython.exe and require a virtual environment. We use virtualenv to implement this.
    
    - *\__init\__.py*  
        Metadata about this project
        
    - *colorBars.py*  
        Functions for creating color bars to be displayed on Paraview images.
    
    - *comboscript.py*    
        Script for collecting interface points into csvs from vtk files and generating images from vtk files. Scripting for many folders and many images and tables.

    - *paraview_csv.py*  
        Functions for collecting interface points and slices inside the nozzle into csvs from vtk files

    - *paraview_general.py*  
        Functions for importing vtk files of simulated filaments.

    - *paraview_line.py*  
        Functions for collecting line traces through the bath in vtk files

    - *paraview_screenshots.py*  
        Functions for generating images of filaments from vtk files.
        
    - *sim_metrics.py* (no longer compatible)  
        Functions for measuring cross-sections in paraview.

- **py/**  
    python tools for generating and analyzing OpenFOAM files. These are written for python3.
    
    - *\__init\__.py*  
        Metadata about this project
        
    - *add_units.py*  
        For adding units to tables generated in Paraview.

    - *folder_parser.py*  
        Tools for ending simulations early at a designated time.
        
    - *folder_scraper.py*    
        Tools for reading legends quickly.
        
    - *folder_stats.py*  
        A class that scrapes metadata from a folder and keeps track of units
        
    - *scrape_tools.py*  
        Functions for scraping values from legends.
        
    - *scrape.py*  
        A class for generating and reading legends for OpenFOAM simulations of embedded 3D printing of single filaments. Written for OpenFOAM v1912 and OpenFOAM 8. Scrapes input files for input variables.
        
    - **cluster_tools/**
        Files that handle simulations while they are running.
    
        - *folder_loop.py*  
            A tool for looping through all folders, running a function, and collecting errors into a table.
        
        - *donescript.py*    
            Functions for moving folders between computers, servers, for OpenFOAM simulations of embedded 3D printing of single filaments.
        
        - *foldermover.py*  
            Functions for moving folders between computers, servers, for OpenFOAM simulations of embedded 3D printing of single filaments.
        
    - **file/**
        For handling files and folders
    
        - *backwards_read.py*   
            Function for reading text files from the bottom
            
        - *file_export.py*  
            Functions for exporting and importing tables and dictionaries
            
        - *file_handling.py*  
            Functions that handle files that are in simulation folders
            
        - *file_names.py*  
            Functions for generating descriptive file names
            
        - *plainIm.py*  
            Functions for exporting and importing tables as pandas dataframes
            
    - **initialize/**
        For creating OpenFOAM initialization files
            
        - *block_points.py*  
            For generating coordinates of blocks used for initializing the mesh
            
        - *block.py*  
            A class that describes blocks to initialize the mesh
            
        - *boundary_input.py*  
            A class that describes boundaries with defined properties
            
        - *cd_vars.py*  
            For handling variables that go into the controlDict
            
        - *compile_XXX.py*  
            For generating specific files in the initialization folder
            
        - *creator_adjacent.py*  
            Functions for generating input files for disturbed filaments or fusion between filaments at a specific set of conditions.
            
        - *creator.py*  
            Functions for generating all files.
            
        - *dict_list.py*  
            A class useful for creating OpenFOAM files in the correct format
            
        - *export.py*  
            Tools for exporting files
            
        - *file_creator.py*  
            A class that holds all input parameters
            
        - *file_group.py*  
            A class that holds all strings to be exported as files
            
        - *file_plotter.py*  
            A class that creates plots representing the input boundaries
            
        - *fluid.py*  
            A class that holds rheological data about a fluid
            
        - *fv_sol_grp.py*  
            A class that holds fvsolution variables for a given solve variable (e.g. interFoam)
            
        - *fv_vars.py*  
            A class that holds all fvsolution and fvschemes variables
            
        - *geometry_file.py*  
            A class that holds a table summarizing the geometry for export
            
        - *initialize_tools.py*  
            A class that creates common headers and footers for files
            
        - *mesh_vars.py*  
            A class that describes the meshing variables
            
        - *ncreate3d.py*  
            Imports dependencies
            
        - *noz_vars.py*  
            Describes the geometry of the nozzle
            
        - *real_boundaries.py*  
            Creates coordinates and boundary conditions for all boundaries
            
        - *transport_group.py*  
            Classes that describe rheology, density, etc. of fluids
            
    - **plot/**
        Tools for plotting results
        
        - *cells_plot.py*  
            Plot the number of cells over time
            
        - *colors.py*  
            Handle colors for all plots
            
        - *combo_plot.py*  
            A class for plotting the same variables across multiple axes, for different sets of data. Reshapes input data so x and y data that is unevenly spaced is plotted categorically instead.
            
        - *convergence_plot.py*  
            Plot residuals over time
            
        - *folder_plots.py*  
            A class for plotting many folders on many plots
            
        - *legend.py*  
            For creating several types of legends
            
        - *markers.py*  
            For determining what markers should be plotted
              
        - *measurement_plot.py*  
            For plotting measurements of cross-sections across different simulations
            
        - *meta_plot.py*  
            A class for plotting metadata
            
        - *multi_plot.py*  
            A class for plotting measurements directly against an input variable.
            
        - *pic_plot.py*  
            For plotting pictures collected using ParaView
            
        - *plot_adjacent.py*  
            Tools for plotting common variables used in simulations of fused filaments and disturbed filaments
            
        - *rate_plot.py*  
            For plotting the simulation rate as circles as a function of input variables
            
        - *sizes.py*  
            For determining the size of the plot, font sizes, based on default settings
             
        - *super_summary_plot.py*  
            For plotting measurements of cross-sections at a specific slice and time, directly against input variables
            
        - *time_plot.py*  
            For plotting simulation time as text as a function of input variables
            
        - *trace_plots.py*  
            For plotting measurements of cross-sections against position in the filament
            
        - *txt_plot.py*  
            For plotting simulation names as a function of input variables
            
        - *value_plot.py*  
            A class for plotting data as squares or circles
            
        - *varPlots.py*  
            Generic class for creating plots and handling the variables to plot against
            
        - *xs_plot.py*  
            Plot cross-sections
            
    - **points/**
        For handling tables of interface points
        
        - *folder_points.py*  
            For importing and exporting points files
            
        - *points_tools.py*  
            For selecting and sorting points
            
        - *slice_points.py*  
            For sorting and labeling points belonging to different segments
            
    - **summarize/**
        For measuring cross-sections and summarizing data across folders
        
        - *ideals.py*  
            Ideal values of cross-section measurements
            
        - *legend_summary.py*  
            For scraping all legends together into one table.
            
        - *log_reader.py*  
            For scraping values from logs
            
        - *steady.py*  
            For determining when simulations reach steady state based on measurements of cross-sections
            
        - *sum_and_steady.py*  
            For summarizing data and finding steady state
            
        - *summarizer_adjacent.py*  
            For measuring cross-sections of disturbed and fused filaments
            
        - *summarizer_single.py*  
            For measuring cross-sections of single filaments
            
        - *summarizer.py*  
            For measuring cross sections, generically
            
        - *super_summary.py*  
            For collecting summary data at a slice and time for all folders and combining into a single table
            
    - **tools/**
        General tools
          
        - *config.py*  
            script for importing packages and setting environmental variables, e.g. folders
            
        - *logs.py*  
            For generating log files
            
        - *strings.py*  
            For converting variable names between their English, code, shortened, and symbolic forms
            
        - *val_tools.py*  
            For converting variables to floats

    - **video/**
        For generating videos
        
        - *video_funcs.py*  
            Tools generating videos
