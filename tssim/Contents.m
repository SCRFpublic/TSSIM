% TSSIM
%
% Files
%   border_nans                     - border_nans Find NaNs connected to the DEM border
%   center_line_proxy               - Considers channel center_line (A) and trend map from sediment source,
%   centerline_generator            - This function generates centerline connecting 'head' and 'tail'
%   channel_shape_generator         - This function generate channel thickness map using a parabola function
%   cont                            - gets contour of matrix a, 2D geobody.
%   convert_to_binary               - This function converts an ascii TI file to binary to be used in SGEMS
%   curvature                       - 8-connected neighborhood curvature of a digital elevation model 
%   dem_flow                        - dem_flow Downslope flow direction for a DEM 
%   dependence_map                  - dependence_map Dependence map for pixel flow in a DEM
%   dssim_lvm                       - This program does direct sequential simulation (via SGEMS) using simple kriging
%   erosion                         - Function that computes the erosion considering the gradient magnitud, 
%   facet_flow                      - facet_flow Facet flow direction
%   fill_sinks                      - fill_sinks Fill interior sinks in DEM
%   find_closest                    - area is a image with zeros (non-area) and ones (area of interest)
%   find_nearest_point              - Find nearest points on centerline defined in d to every point on grid 
%   flow_matrix                     - flow_matrix System of linear equations representing pixel flow
%   geoeas2matlab                   - Transform GeoEas array into Matlab image
%   get_angle                       - gets angle used in the object-based simulation in the lobe orientation. The 
%   get_core_data                   - gets core data located in the simulation area
%   get_core_data_in_SA             - retrieves information regarding core in SA.
%   get_dirname                     - Returns fullpath
%   get_fgu                         - Gets fine-grained unit considering an altitude component. Thickner at low 
%   get_model                       - computes grid model for geoeas software (SGeMS) usage 
%   get_model_slice                 - gets grid model for slice number nslice
%   get_new_start_point_indep       - This function uses a probability map to get the next lobe start point
%   get_simulation_area             - function that gives simulations area.
%   get_simulation_area_conditional - function that gives simulations area that fits the hard data (data)
%   get_source_location             - Finds sediment source location in the range.
%   gradient8                       - 8-connected neighborhood gradient and aspect of a digital elevation model
%   influence_map                   - influence_map Influence map for pixel flow in a DEM
%   lmv_surface                     - smoothes out surface E using local moving windows of windows_size.
%   lobe_cond_thickness_assignment  - Assigns conditional thickness to binary 2D geobody (closeBW).
%   lobe_generator                  - Functions that generates lobate part of the geobody. 
%   lobe_thickness_assignment       - Assigns thickness to binary 2D geobody (closeBW). Returns thickness map
%   make_grid                       - Calculate grid dimensions
%   make_pfield_for_startpoint_tau2 - computes p-field for anchor point.
%   matlab2geoeas                   - Transform Matlab array into GeoEAS column vector
%   model_movie                     - gives animation of the model slice
%   modify_cdf                      - Used by tssim to modify cdf by extrapolating endpoints.
%   noisy_data                      - Writes noisy data called for .par file
%   pixel_flow                      - pixel_flow Downslope flow direction for DEM pixels
%   plot_model_slice                - plots slice number n
%   plotsim                         - Plots surface cube output from tssim
%   plotting_commands               - function plotting_commands
%   postprocess_plateaus            - postprocess_plateaus Replace upslope areas for plateaus with mean value 
%   read_cdf_data                   - Load pattern parameter cdf data file
%   read_earthvision_grid           - This function reads a grid output by Petrel in EarthVision format
%   read_harddata                   - reads hard data from text files in tssim\Well_data\
%   read_SGEMS_output               - Read in SGEMS simulation
%   river_direction                 - Returns river flow direction according to the following pattern:
%   river_length                    - Used by tssim to compute river length
%   run_SGEMS_sim                   - This function runs SGEMS in command mode from MATLAB
%   statistics                      - script to compute some statistics of surface cube output from tssim
%   trend_map_maker                 - Program generates trend map from sediment source until b.
%   tssim                           - Main script for surface-based turbidite lobe system simulation
%   upslope_area                    - upslope_area Upslope area measurements for a DEM
%   vis_dem_flow                    - vis_dem_flow Visual pixel flow directions and magnitude on a DEM
%   vis_map                         - vis_map Visualize influence or dependence map on a DEM
%   warpImage                       - Warps Image so that orinal marked points are at the desired marked points and all other points 
%   write_dimension_parfile         - Writes .par file for sgsim to add noise to the simulation. Only for non-
%   write_dimension_parfile_depo    - Writes .par file for sgsim to add noise to the simulation. Only for non-
%   write_dimension_parfile_erosion - Writes .par file for sgsim to add noise to the simulation. Only for non-
%   write_dimension_parfile_fgu     - Writes .par file for sgsim to add noise to the simulation. Only for non-
%   write_gslib_grid                - write a grid in GEOEAS/GSLIB format to a named file
