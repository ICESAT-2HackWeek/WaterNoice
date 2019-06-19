# WaterNoice

### Data Science Objective
Do the residual histograms of ATL06 data vary across different hydrologic features (meltwater ponds, 
firn aquifer seeps, blue ice megadunes, etc.)?
And if so, characterize this variability through space and time.

### Datasets
- ICESat-2 [ATL06](https://nsidc.org/data/atl06?qt-data_set_tabs=3#qt-data_set_tabs) (40 m resolution)
- REMA: 2 m digital elevation data for Antarctica
- Sentinel-1 and LANDSAT imagery data for identifying perennial surface features

### Validation Datasets
ArcticDEM, REMA, or other local DEM can serve as a comparison benchmark.  How do the ATL06 data differ from these baselines, in both a mean sense and (especially) the variability within a single flyover (histogram width)?

### Tools
adapt jupyter notebooks provided by instructors into a python workflow

### Tasks
- Learn how to download the ICESat-2 data by lat lon bounding box
- explore parameters associated with level 6 data - gain familarity 
- develop pipeline for subsetting data and downloading ICESat-2 ATL06
- develop first order filtering for each region (keep in mind if specific filtering parameters vary for different complex       surfaces)
- exploratory statistics that describe differences in key parameters between different water features
- develop python code for plotting results
- what are the big takeaways? Keep a running document of key findings to share with hackweek group at end of week

### Study Sites:
- Store Glacier, west Greenland
- Firn aquifer region in southeast Greenland (near Helheim)
- Amery ice shelf melt water ponds
- RBIS (east Antarctica)
-

