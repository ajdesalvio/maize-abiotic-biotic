import geopandas as gpd
from shapely.geometry import Polygon
import shapely.affinity
import pandas as pd
import math
import os

def create_shapefile_from_csv(
    csv_file: str,
    shapefile_output: str,
    A: tuple,
    B: tuple,
    rangelength: float,
    rangebuffer: float,
    rowwidth: float,
    rowbuffer: float,
    utm_epsg: int,
    plot_id_colname: str,
    range_colname: str,
    row_colname: str,
    ft_toggle: float = 0.3048,
    shift_rows_to_one: bool = True,   
    remap_rows: bool = True
) -> gpd.GeoDataFrame:
    """
    Create a shapefile from a CSV file that defines plot dimensions.
    
    Parameters
    ----------
    csv_file : str [required]
        Path to the CSV file containing Range and Row definitions. Range is the Y axis of the field and Row is the X axis. There should also be
        a column name (e.g., Plot_ID, Plant_ID) that containes UNIQUE identifiers for each experimental unit.
    shapefile_output : str [required]
        Output path for the shapefile (e.g., 'C:/Users/username/output.shp').
    A : tuple of float [required]
        Coordinates (x, y) of the starting point (A) for the AB line. Obtain these from QGIS by right clicking on the A point.
        This most likely needs to be coordinates in the same EPSG code that you will enter below for the utm_epsg argument.
        For example, College Station, TX, is in UTM zone 14N, which has an EPSG code of 32614. Right clicking in QGIS
        will bring up a dialog box that says "Copy Coordinate" and multiple options will be presented. In this case, select
        the option that corresponds to the correct EPSG code. As a specific example, the option for UTM zone 14N would be to copy
        the coordinates for "WGS 84 / UTM zone 14N", and then paste them within parentheses such as: (752856.1989,3390459.0102).
        If your data were exported from your photogrammetry software and were not saved in your specific UTM zone, go back and save them 
        according to the UTM zone you collected the data in. An example of how to find your UTM zone is given in the utm_epsg argument below.
    B : tuple of float [required]
        Coordinates (x, y) of the ending point (B) for the AB line. Obtain these from QGIS by right clicking on the B point as with the A point.
    rangelength : float [required]
        The TOTAL length of each plot. This is end-to-end assuming there's no space between plots.
    rangebuffer : float [required]
        Reduction on each side of each plot along the "range" dimension. E.g., if your total range length is 25 feet end-to-end but you have
        4 feet of border between the end of one plot and the beginning of the next plot, this means you have a 2-food reduction on each side of
        your plots. You would enter 2 for the rangebuffer parameter.
    rowwidth : float [required]
        The TOTAL width of each plot along the "row" dimension. This is side-to-side assuming there's no space between plots.
    rowbuffer : float [required]
        Reduction on each side of each plot along the "row" dimension. The same logic applies as with rangelength and rangebuffer. If your plots are 
        5 feet wide but you have 0.6 feet (7.2 inches) of spacing from the end of one plot before the next plot begins, this corresponds to a 0.3-foot reduction 
        on each side of your plots. You would enter 0.3 for the rowbuffer parameter.
    utm_epsg : int [required]
        EPSG code for the coordinate system (e.g., 32614 for UTM zone 14N). Obtain from a site like https://mangomap.com/robertyoung/maps/69585/what-utm-zone-am-i-in-
    plot_id_colname : str [required]
        Column name in the CSV identifying the unique plot or plant ID.
    range_colname : str [required]
        Column name in the CSV indicating the "range" index (Y axis of the field).
    row_colname : str [required]
        Column name in the CSV indicating the "row" index (X axis of the field).
    ft_toggle : float [optional]
        Conversion factor for feet-to-meters. Use 1 if your data is in meters already.
    shift_rows_to_one : bool [optional] but highly recommended. 
        If True, shift the row_colname values so the smallest value becomes 1
        before optionally doing the remap that's available with the remap_rows argument below. 
        Defaults to True. For example, if your first row starts at 300, the first time you drag 
        and drop the newly created shapefile into QGIS, it will appear far off to the right of where 
        the A coodinate is in the AB line. Setting this argument to True will make it so that rows
        that begin with a number other than one, e.g. (300, 301, 302) will be remapped to (1, 2, 3).
    remap_rows : bool [optional]
        If True, remap the values in the row_colname column to consecutive integers
        while preserving sorted order. Defaults to True.
        If you have odd-numbered row indices, there will be unwanted space between plots. This
        parameter mitigates this issue by taking odd or even-numbered row names (such as 1, 3, 5)
        and remapping them to consecutive integers (1, 2, 3).
    
    Returns
    -------
    gpd.GeoDataFrame
        A GeoDataFrame containing the created polygons. Also writes out a shapefile.
    """

    # Read CSV
    df = pd.read_csv(csv_file)
    
    # Shift rows so that the samllest value becomes 1, if requested
    # Example: if row values are [331, 333, 335], min is 331, so subtract (331 - 1) = 330 => [1, 3, 5]
    if shift_rows_to_one:
        min_val = df[row_colname].min()  # e.g. 331
        shift_value = min_val - 1        # e.g. 330
        df[row_colname] = df[row_colname] - shift_value

    # 3. If requested, remap row_colname values to consecutive integers
    if remap_rows:
        unique_rows = sorted(df[row_colname].unique())
        row_map = {old_val: i + 1 for i, old_val in enumerate(unique_rows)}
        df[row_colname] = df[row_colname].map(row_map)

    # 4. Calculate angle of the AB line relative to the x-axis
    dx = B[0] - A[0]
    dy = B[1] - A[1]
    angle = math.degrees(math.atan2(dy, dx))

    polygons = []
    plot_ids = []

    # 5. Build polygons for each row in the DataFrame
    for _, row_data in df.iterrows():
        
        plot_id = row_data[plot_id_colname]
        range_number = row_data[range_colname]
        row_number = row_data[row_colname]

        # Calculate length and width (minus buffers)
        plot_length = (rangelength * ft_toggle) - 2 * (rangebuffer * ft_toggle)
        plot_width = (rowwidth * ft_toggle) - 2 * (rowbuffer * ft_toggle)

        # Calculate the center coordinates for this plot
        x_center = A[0] + (range_number * rangelength * ft_toggle)
        y_center = A[1] - (row_number * rowwidth * ft_toggle)

        # Define rectangle corners
        x1 = (x_center - plot_length / 2)
        x2 = (x_center + plot_length / 2)
        y1 = (y_center - plot_width / 2)
        y2 = (y_center + plot_width / 2)

        # Create and rotate rectangle
        polygon = Polygon([(x1, y1), (x2, y1), (x2, y2), (x1, y2), (x1, y1)])
        rotated_polygon = shapely.affinity.rotate(polygon, angle, origin=A)

        polygons.append(rotated_polygon)
        plot_ids.append(plot_id)

    # 6. Create a GeoDataFrame
    gdf = gpd.GeoDataFrame(
        {plot_id_colname: plot_ids, 'geometry': polygons},
        crs=f"EPSG:{utm_epsg}"
    )

    # 7. Write to shapefile and return the GeoDataFrame
    gdf.to_file(shapefile_output)
    return gdf


#### Example usage ####
# Path to CSV
csv_path = 'C:/Users/aaron.desalvio/Downloads/Rust_Data'
# Write the path where you'd like to save the shapefile
output_path = 'C:/Users/aaron.desalvio/Downloads/Rust_Data'

gdf = create_shapefile_from_csv(
    csv_file=os.path.join(csv_path, '2021_Rust_Scores_Range_Row_V2.csv'),
    shapefile_output=os.path.join(output_path, '20210716_Rust.shp'),
    A=(746701.776,3381975.148),
    B=(746832.897,3381851.032),
    ft_toggle=0.3048,
    rangelength=25,
    rangebuffer=1.5,
    rowwidth=5.04,
    rowbuffer=0.3,
    utm_epsg=32614,
    plot_id_colname='Range_Row',
    range_colname='Range',
    row_colname='Row.Odd',
    remap_rows=True,
    shift_rows_to_one=True
)