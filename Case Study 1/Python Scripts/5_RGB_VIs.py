import numpy as np
import os
import cv2
import concurrent.futures
import multiprocessing.reduction as reduction
import cloudpickle
import csv

reduction.ForkingPickler.dumps = cloudpickle.dumps

def get_image_file_paths(folder):
    """
    Instead of loading images, just list file paths and base filenames.
    """
    file_paths = []
    img_files = []
    valid_extensions = ('.png', '.jpg', '.jpeg', '.tif', '.tiff')
    for filename in os.listdir(folder):
        if filename.lower().endswith(valid_extensions):
            file_paths.append(os.path.join(folder, filename))
            img_files.append(os.path.splitext(filename)[0])
    return file_paths, img_files

def safe_mean_median(array):
    valid = ~np.isnan(array) & ~np.isinf(array)
    if np.any(valid):
        array = array[valid]
        return np.mean(array), np.median(array)
    else:
        return np.nan, np.nan

# Default scale parameter
scale = 255

# Define vegetation indices as lambda functions (RGB version)
base_vi_functions = {
    'BCC':      lambda scale, R, G, B: (B/scale) / ((R/scale) + (G/scale) + (B/scale)),
    'BGI':      lambda scale, R, G, B: (B/scale) / (G/scale),
    'BI':       lambda scale, R, G, B: np.sqrt(((R/scale)**2 + (G/scale)**2 + (B/scale)**2) / 3),
    'CIVE':     lambda scale, R, G, B: (0.441*(R/scale) - 0.811*(G/scale) + 0.385*(B/scale) + 18.78745),
    'COM1':     lambda scale, R, G, B: ((2*(G/scale) - (R/scale) - (B/scale)) 
                                 + (0.441*(R/scale) - 0.811*(G/scale) + 0.385*(B/scale) + 18.78745) 
                                 + (3*(G/scale) - 2.4*(R/scale) - (B/scale)) 
                                 + ((G/scale) / ((R/scale)**0.667 * (B/scale)**0.334))),
    'COM2':     lambda scale, R, G, B: (0.36*(2*(G/scale) - (R/scale) - (B/scale)) 
                                 + 0.47*(0.441*(R/scale) - 0.811*(G/scale) + 0.385*(B/scale) + 18.78745)
                                 + 0.17*((G/scale) / ((R/scale)**0.667 * (B/scale)**0.334))),
    'ExG':      lambda scale, R, G, B: (2*(G/scale) - (R/scale) - (B/scale)),
    'ExG2':     lambda scale, R, G, B: (2*(G/scale) - (R/scale) - (B/scale)) / ((R/scale) + (G/scale) + (B/scale)),
    'ExGR':     lambda scale, R, G, B: (3*(G/scale) - 2.4*(R/scale) - (B/scale)),
    'EXR':      lambda scale, R, G, B: (1.4*(R/scale) - (G/scale)),
    'GminusB':  lambda scale, R, G, B: ((G/scale) - (B/scale)),
    'GminusR':  lambda scale, R, G, B: ((G/scale) - (R/scale)),
    'GdivB':    lambda scale, R, G, B: (G/scale) / (B/scale),
    'GdivR':    lambda scale, R, G, B: (G/scale) / (R/scale),
    'GCC':      lambda scale, R, G, B: (G/scale) / ((R/scale) + (G/scale) + (B/scale)),
    'GLI':      lambda scale, R, G, B: (2*(G/scale) - (R/scale) - (B/scale)) / (2*(G/scale) + (R/scale) + (B/scale)),
    'MExG':     lambda scale, R, G, B: (1.262*(G/scale) - 0.884*(R/scale) - 0.311*(B/scale)),
    'MGVRI':    lambda scale, R, G, B: (((G/scale)**2 - (R/scale)**2) / ((G/scale)**2 + (R/scale)**2)),
    'MRCCbyAlper': lambda scale, R, G, B: ((R/scale)**3) / ((R/scale) + (G/scale) + (B/scale)),
    'MSRGR':    lambda scale, R, G, B: np.sqrt((G/scale) / (R/scale)),
    'NDI':      lambda scale, R, G, B: 128 * ((((G/scale) - (R/scale)) / ((G/scale) + (R/scale))) + 1),
    'NDRBI':    lambda scale, R, G, B: ((R/scale) - (B/scale)) / ((R/scale) + (G/scale)),
    'NGBDI':    lambda scale, R, G, B: ((G/scale) - (B/scale)) / ((G/scale) + (B/scale)),
    'NGRDI':    lambda scale, R, G, B: ((G/scale) - (R/scale)) / ((G/scale) + (R/scale)),
    'NRMBI':    lambda scale, R, G, B: ((R/scale) - (B/scale)) / (G/scale),
    'RCC':      lambda scale, R, G, B: (R/scale) / ((R/scale) + (G/scale) + (B/scale)),
    'RdivB':    lambda scale, R, G, B: (R/scale) / (B/scale),
    'RGBVI':    lambda scale, R, G, B: (((G/scale)**2 - (R/scale)*(B/scale)) / ((G/scale)**2 + (R/scale)*(B/scale))),
    'RminusB':  lambda scale, R, G, B: ((R/scale) - (B/scale)),
    'TGI':      lambda scale, R, G, B: ((G/scale) - (0.39*(R/scale) - 0.69*(B/scale))),
    'TNDGR':    lambda scale, R, G, B: np.sqrt((((G/scale) - (R/scale)) / ((G/scale) + (R/scale))) + 0.5),
    'VARI':     lambda scale, R, G, B: ((G/scale) - (R/scale)) / ((G/scale) + (R/scale) - (B/scale)),
    'VEG':      lambda scale, R, G, B: (G/scale) / ((R/scale)**0.667 * (B/scale)**0.334)
}

# --- Helper Function for Parallel Processing ---
def process_single_rgb_image(args):
    """
    Processes a single image from file, loads it on the fly, and returns its VI stats
    and (optionally) distribution rows.
    """
    file_path, fname, final_vi_dict, distributions, probability_steps, scale_param = args
    # Load image on the fly
    img = cv2.imread(file_path, cv2.IMREAD_COLOR)
    if img is None:
        return {'file_name': fname, 'total_pixels': 0, 'black_pixels': 0, 'non_black_pixels': 0}, []
    
    height, width, _ = img.shape
    total_pixels = height * width
    # Identify black pixels: a pixel is black if all channels are zero.
    black_pixel_mask = np.all(img == 0, axis=2)
    black_pixels = np.sum(black_pixel_mask)
    non_black_pixels = total_pixels - black_pixels

    row_data = {
        'file_name': fname,
        'total_pixels': total_pixels,
        'black_pixels': black_pixels,
        'non_black_pixels': non_black_pixels
    }
    distribution_rows = []

    if non_black_pixels > 0:
        # Extract channels as float (cv2 loads in BGR order)
        B = img[..., 0].astype(float)
        G = img[..., 1].astype(float)
        R = img[..., 2].astype(float)
        # Mask out black pixels
        mask_non_black = ~black_pixel_mask
        R = R[mask_non_black]
        G = G[mask_non_black]
        B = B[mask_non_black]

        for vi_name, vi_func in final_vi_dict.items():
            vi_values = vi_func(scale_param, R, G, B)
            # Clean up infinities or NaNs
            vi_values[np.isinf(vi_values)] = np.nan
            vi_values[np.isnan(vi_values)] = np.nan
            mean_val, median_val = safe_mean_median(vi_values)
            row_data[f'{vi_name}_mean'] = mean_val
            row_data[f'{vi_name}_median'] = median_val

            if distributions:
                valid_vals = vi_values[~np.isnan(vi_values)]
                if valid_vals.size == 0:
                    quantiles = np.full(probability_steps.shape, np.nan)
                else:
                    quantiles = np.quantile(valid_vals, probability_steps)
                for prob, quant in zip(probability_steps, quantiles):
                    distribution_rows.append({
                        'file_name': fname,
                        'vi_name': vi_name,
                        'probability': prob,
                        'quantile': quant
                    })
    else:
        for vi_name in final_vi_dict.keys():
            row_data[f'{vi_name}_mean'] = np.nan
            row_data[f'{vi_name}_median'] = np.nan

    return row_data, distribution_rows

# --- Parallel RGB VI Computation Function with Batch Processing and Streaming ---
def compute_rgb_vegetation_indices_parallel(
    image_folder: str,
    output_folder: str,
    output_csv_name: str,
    vi_functions: dict = None,
    replace_vi_dictionary: bool = False,
    distributions: bool = False,
    output_csv_name_distribution: str = None,
    max_workers: int = None,
    batch_size: int = 100,  # Process tasks in batches to reduce memory pressure
    scale_param: float = 255
):
    """
    Computes vegetation indices from images in image_folder using parallel processing,
    and streams results directly to CSV files in output_folder.
    Both the main VI stats and distribution data are written immediately to disk.
    Tasks are processed in batches to avoid accumulating too many future objects and to recycle worker processes.

    Parameters
    ----------
    image_folder : str [Required]
        Path to the folder containing input image files.
    output_folder : str [Required]
        Path to the directory where the output CSV files will be saved.
    output_csv_name : str [Required]
        Filename for the main vegetation index statistics CSV (mean and median values only; no distributional data).
    vi_functions : dict [Optional]
        Custom dictionary of vegetation index functions to extend or override the defaults.
    replace_vi_dictionary : bool [Optional]
        If True, the provided vi_functions will entirely replace the default functions. If False and VIs are
        provided in vi_functions, then these user-defined VIs will be appended to the default list.
    distributions : bool [Optional]
        If True, computes and streams quantile distribution data for each vegetation index.
    output_csv_name_distribution : str [Optional]
        Filename for the CSV containing distribution statistics (used if distributions is True).
    max_workers : int [Optional]
        Maximum number of worker processes to use for parallel processing.
    batch_size : int [Optional]
        Number of image tasks processed per batch to minimize memory overhead.
    scale_param : float [Optional]
        Scaling factor used to normalize pixel values (e.g., 255 for RGB images).

    Returns
    -------
    pd.DataFrame
        DataFrame containing the main vegetation index statistics reloaded from the output CSV.
        (Note: The results are primarily written to disk to maintain memory efficiency.)

    """
    # Decide which VI dictionary to use
    if replace_vi_dictionary:
        final_vi_dict = vi_functions if vi_functions is not None else base_vi_functions
    else:
        final_vi_dict = {**base_vi_functions, **(vi_functions or {})}
    
    # Prepare file paths and probability steps (if distributions are requested)
    file_paths, img_files = get_image_file_paths(image_folder)
    probability_steps = np.linspace(0.01, 1.00, 100) if distributions else None

    # Prepare CSV file for streaming main VI stats
    os.makedirs(output_folder, exist_ok=True)
    output_path = os.path.join(output_folder, output_csv_name)
    
    # Define main fieldnames based on final_vi_dict keys
    main_fieldnames = ["file_name", "total_pixels", "black_pixels", "non_black_pixels"]
    for vi in final_vi_dict.keys():
        main_fieldnames.extend([f'{vi}_mean', f'{vi}_median'])
    
    main_file = open(output_path, "w", newline="")
    main_csv_writer = csv.DictWriter(main_file, fieldnames=main_fieldnames)
    main_csv_writer.writeheader()
    
    # Prepare CSV file for streaming distribution rows if requested
    distribution_file = None
    distribution_csv_writer = None
    if distributions and output_csv_name_distribution is not None:
        output_path_dist = os.path.join(output_folder, output_csv_name_distribution)
        distribution_file = open(output_path_dist, "w", newline="")
        dist_fieldnames = ["file_name", "vi_name", "probability", "quantile"]
        distribution_csv_writer = csv.DictWriter(distribution_file, fieldnames=dist_fieldnames)
        distribution_csv_writer.writeheader()

    # Build all tasks
    tasks = [
        (file_path, fname, final_vi_dict, distributions, probability_steps, scale_param)
        for file_path, fname in zip(file_paths, img_files)
    ]
    
    # Process tasks in batches to reduce the number of futures and recycle workers
    for i in range(0, len(tasks), batch_size):
        batch_tasks = tasks[i:i+batch_size]
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            # Using executor.map minimizes the accumulation of Future objects
            for row_data, dist_rows in executor.map(process_single_rgb_image, batch_tasks):
                # Write the main VI stats row immediately
                main_csv_writer.writerow(row_data)
                main_file.flush()
                # Write distribution rows immediately if applicable
                if distribution_csv_writer is not None:
                    for row in dist_rows:
                        distribution_csv_writer.writerow(row)
                    distribution_file.flush()

    # Close CSV files
    main_file.close()
    print(f"Saved VI stats to: {output_path}")
    if distribution_file is not None:
        distribution_file.close()
        print(f"Saved VI distribution stats to: {output_path_dist}")

# To run the script using the base VI functions, highlight the code starting with the if__name__ portion down to the
# end of the first compute_rgb_vegetation_indices_parallel function and run it.
if __name__ == "__main__":
    # Path to soil-masked images
    image_folder = 'path/to/soil-masked/images'
    # Path where you'd like the vegetation indices to be saved
    output_folder = 'path/to/save/VIs'
    # Name of the vegetation index file
    output_csv_name = '2021_Rust_VIs.csv'

    # Determine the number of workers
    safe_num_workers = os.process_cpu_count() - 12  # adjust as needed

    # First run using base VI functions
    compute_rgb_vegetation_indices_parallel(
        image_folder=image_folder,
        output_folder=output_folder,
        output_csv_name=output_csv_name,
        vi_functions=None,
        replace_vi_dictionary=False,
        max_workers=safe_num_workers,
        batch_size=100,
        scale_param=255
    )
    
    # Example of defining a custom VI:
    custom_vi_dict = {
        "SimpleIndex": lambda R, G, B: (R/scale) - (G/scale)
    }


    # Example: merging custom VI with base functions
    compute_rgb_vegetation_indices_parallel(
        image_folder=image_folder,
        output_folder=output_folder,
        output_csv_name="Additional_VI.csv",
        vi_functions=custom_vi_dict,
        replace_vi_dictionary=False,
        max_workers=4,
        batch_size=100,
        scale_param=255
    )

    # Example: using only custom VI (replacing base functions)
    compute_rgb_vegetation_indices_parallel(
        image_folder=image_folder,
        output_folder=output_folder,
        output_csv_name="User_VI_Only.csv",
        vi_functions=custom_vi_dict,
        replace_vi_dictionary=True,
        max_workers=4,
        batch_size=100,
        scale_param=255
    )
