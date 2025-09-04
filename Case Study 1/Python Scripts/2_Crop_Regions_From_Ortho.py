import os
from typing import List
import geopandas as gpd
import rasterio
from rasterio.mask import mask

def crop_regions_from_ortho(
    flights: List[str],
    shp_folder: str,
    shp_base: str,
    tif_folder: str,
    tif_base: str,
    output_path: str,
    plot_id_colname: str
) -> None:
    """
    Crops regions from a set of orthomosaic TIFFs based on corresponding shapefiles.
    
    Parameters
    ----------
    flights : List[str] [required]
        A list of flight dates (e.g., ['20220505', '20220512', '20220520']).
    shp_folder : str [required]
        Path to the folder containing the shapefiles.
    shp_base : str [required]
        Common part of the shapefile name after the flight date (e.g., '_Texas_Plots.shp'). This allows for iteration
        through multiple shapefiles that are named the same way (e.g., 20220505_Texas_Plots.shp, 20220512_Texas_Plots.shp)
    tif_folder : str [required]
        Path to the folder containing the TIFFs.
    tif_base : str [required]
        Common part of the orthomosaic TIFF name after the flight date (e.g., '_Texas_Ortho.tif').
    output_path : str [required]
        Folder path where the output cropped files will be saved.
    plot_id_colname : str [required]
        ESRI shapefiles are in essence data frames; this argument should specify the column within the data frame
        that indicates the plot ID or plant ID. It's a unique identifier for each experimental unit in the shapefile.
        If you're unsure what this column name is called, you can drag and drop a shapefile into QGIS, right click
        on its name in the "Layers" tab, then open the "Attributes Table". This will allow you to view all the columns
        present in the shapefile.
    
    Returns
    -------
    None
    """
    
    # Ensure the output directory exists
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    for flight in flights:
        # Construct shapefile path
        shapefile_path = os.path.join(shp_folder, f"{flight}{shp_base}")
        
        # Read the shapefile
        shapefile = gpd.read_file(shapefile_path)
        
        # Construct raster path
        raster_path = os.path.join(tif_folder, f"{flight}{tif_base}")
        
        # Open the raster
        with rasterio.open(raster_path) as src:
            print(f"Processing flight: {flight}")
            print("Raster CRS: ", src.crs)
            print("Shapefile CRS: ", shapefile.crs)

            # If needed, reproject shapefile to match raster CRS
            if shapefile.crs != src.crs:
                 shapefile = shapefile.to_crs(src.crs)
                 print("Shapefile reprojected to match raster CRS.")
            
            # Loop through each geometry in the shapefile
            for _, row in shapefile.iterrows():
                
                # Perform masking/cropping
                out_image, out_transform = mask(
                    dataset=src,
                    shapes=[row['geometry']],
                    crop=True
                )
                
                # Copy the metadata and update height, width, and transform
                out_meta = src.meta.copy()
                out_meta.update({
                    "driver": "GTiff",
                    "height": out_image.shape[1],
                    "width": out_image.shape[2],
                    "transform": out_transform
                })

                # Build the output file path
                plot_id = row[plot_id_colname]
                output_file = os.path.join(output_path, f"{flight}_{plot_id}.tif")
                
                # Write the cropped raster to disk
                with rasterio.open(output_file, "w", **out_meta) as dest:
                    dest.write(out_image)
                
                print(f"Saved cropped image for plot {plot_id} at {output_file}")

# Folder where TIFs are stored
tif_folder = 'path/to/TIFs'
# Obtain all flight dates (if named with the flight dates at the start of the file name)
flight_list = [f[:8] for f in os.listdir(tif_folder) if f.lower().endswith(".tif")]
# Alternatively, manually type out the flight dates in a list
#flight_list  = ['20210606', '20210615']
# Folder where shapefiles are stored
shp_folder = 'path/to/shapefiles'
# Portion of the shapefile file name that remains constant
shp_base = '_Rust.shp'
# Portion of the TIF file name that remains constant
tif_base = '_Rust.tif'
# Path where you'd like to save the cropped plot images
output_path = 'path/for/saving/cropped/plots'
# Column name (within the shapefile) that provides the unique plot identifier
plot_id_colname = 'Range_Row'

# Use the function
crop_regions_from_ortho(
    flights=flight_list,
    shp_folder=shp_folder,
    shp_base=shp_base,
    tif_folder=tif_folder,
    tif_base=tif_base,
    output_path=output_path,
    plot_id_colname=plot_id_colname)
