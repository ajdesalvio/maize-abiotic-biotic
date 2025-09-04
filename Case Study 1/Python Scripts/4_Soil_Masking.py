import numpy as np
import os
import cv2


def load_images_from_folder(folder):
    images = []
    img_files = []
    # Filter to image files only to avoid errors
    valid_extensions = ('.png', '.jpg', '.jpeg', '.tif', '.tiff')
    for filename in os.listdir(folder):
        if filename.lower().endswith(valid_extensions):
            img_path = os.path.join(folder, filename)
            img = cv2.imread(img_path, cv2.IMREAD_COLOR)
            if img is not None:
                images.append(img)
                # Store filename without extension
                img_files.append(os.path.splitext(filename)[0])
    return images, img_files


##### ExG Soil Removal Function #####
# Soil removal index calculation based on ExG (Credit: Travis Parker https://youtu.be/U2qJn7jHYDw?si=nrF9iz1t8wKxJjot)
def rem_soil_index(img, 
                   threshold=0.1, 
                   index_name='ExG', 
                   comparison='less',
                   custom_indices=None):
    """
    Compute a vegetation index for each pixel and remove pixels (set to np.nan) 
    based on a threshold. Also removes all-black pixels.

    Parameters
    ----------
    img : numpy.ndarray
        Input image in BGR format (cv2 imports images in BGR order)
    threshold : float
        Threshold used for filtering.
    index_name : str
        Which index function to use from 'custom_indices' if the user is supplying
        their own index.
    comparison : str
        'less' or 'greater' - direction of threshold check.
    custom_indices : dict
        A dictionary of lambdas {index_name: lambda R,G,B: ...}.
        If None, defaults to a dict containing only 'ExG'.

    Returns
    -------
    np.ndarray:
        A 2D array of the computed vegetation index with pixels that exceed or fall
        below the threshold OR problematic pixels (those with a divide by zero
        error) are set to np.nan.
    """

    # If no custom dictionary is provided, use a default dictionary with ExG
    if custom_indices is None:
        custom_indices = {
            'ExG': lambda R, G, B: 2*(G/(R+G+B)) - (R/(R+G+B)) - (B/(R+G+B))
        }
    
    # Check if the requested index_name is in the dictionary
    if index_name not in custom_indices:
        raise ValueError(f"Index type '{index_name}' not found in custom_indices. "
                         f"Available keys: {list(custom_indices.keys())}")

    # Extract individual channels (B, G, R) - note that this band ordering was determined
    # by using cv2 to import images, which defaults to BGR
    B = img[:, :, 0].astype(np.float32)
    G = img[:, :, 1].astype(np.float32)
    R = img[:, :, 2].astype(np.float32)

    # Calculate the vegetation index using the selected lambda
    remsoil_index = custom_indices[index_name](R, G, B)
    # Divide by zero instances will produce np.nan or np.inf automatically, which are removed below.

    # Identify all-black pixels (these need to be removed too)
    black_pixel_mask = (R == 0) & (G == 0) & (B == 0)

    # Apply threshold in the chosen direction
    if comparison == 'less':
        remove_mask = remsoil_index < threshold
    elif comparison == 'greater':
        remove_mask = remsoil_index > threshold
    else:
        raise ValueError("comparison must be 'less' or 'greater'.")

    # Mark black pixels OR threshold-failing pixels as np.nan
    remsoil_index[black_pixel_mask | remove_mask] = np.nan

    return remsoil_index


def create_binary_mask(remsoil_index):
    """
    Create a binary mask from the Hue index.
    Valid pixels are 1, and NA (numpy.nan) values are 0.

    Parameters
    ----------
    remsoil_index : numpy array
        Vegetation index with possible NA values.

    Returns
    -------
    numpy array
        Binary mask (1 for valid, 0 for invalid).
    """
    binary_mask = np.zeros_like(remsoil_index, dtype=np.uint8)
    valid_mask = ~np.isnan(remsoil_index)
    binary_mask[valid_mask] = 1
    return binary_mask


def process_and_save_masks(input_images: str,
                           masks_out: str,
                           masked_images_out: str,
                           threshold: float = 0.1,
                           index_name: str = 'ExG',
                           comparison: str = 'less',
                           custom_indices = None,
                           output_type: str = '.jpg'):
    """
    1) Loads images from 'input_images' folder.
    2) Computes vegetation index for each pixel within each image using rem_soil_index.
    3) Creates a binary mask and saves it to 'masks_out' folder.
    4) Applies the mask to the image and saves the masked result to 'masked_images_out'.

    Parameters
    ----------
    input_images : str [required]
        Path to folder containing images with soil present.
    masks_out : str [Required]
        Path to folder where binary masks will be saved.
    masked_images_out : str [Required]
        Path to folder where masked images will be saved.
    threshold : float [Required]
        Threshold for rem_soil_index. Default is 0.1 corresponding to the ExG index.
    index_name : str [Required]
        Name of the index to use from custom_indices. Default is 'ExG'.
    comparison : str [Required]
        'less' or 'greater'. Default is 'less', again corresponding to the ExG index.
    custom_indices : dict [Optional] 
        dict of lambdas {index_name: lambda R,G,B: ...}. 
        If None, will use a default dict with only ExG.
    output_type : str [Required]
        File extension for output images (e.g., '.jpg', '.tif'). Default is '.jpg'.
        
    Returns
    -------
    None
    """

    # Ensure output directories exist
    os.makedirs(masks_out, exist_ok=True)
    os.makedirs(masked_images_out, exist_ok=True)

    # Load images
    images, img_files = load_images_from_folder(input_images)

    # Process each image
    for i, img in enumerate(images):
        # Calculate index
        remsoil = rem_soil_index(
            img, 
            threshold=threshold, 
            index_name=index_name, 
            comparison=comparison,
            custom_indices=custom_indices
        )

        # Create binary mask
        binary_mask = create_binary_mask(remsoil)

        # Save the mask to a file
        mask_output_path = os.path.join(masks_out, f"{img_files[i]}_mask{output_type}")
        cv2.imwrite(mask_output_path, (binary_mask * 255).astype(np.uint8))

        # Apply mask to the original image
        masked_img = img.copy()
        for channel in range(3):
            masked_img[:, :, channel] *= binary_mask

        # Save the masked image
        masked_img_output_path = os.path.join(
            masked_images_out, f"{img_files[i]}_masked{output_type}"
        )
        cv2.imwrite(masked_img_output_path, masked_img)

        print(f"Processed and saved mask for {img_files[i]} at {mask_output_path}")

# Path to cleaned cropped images from the previous script
input_images = 'path/to/cleaned/cropped/images'
# Path where you'd like to save the binary masks
masks_out = 'path/to/save/masks'
# Path where you'd like the masked images to be saved
masked_images_out = 'path/to/masked/images'

# If specifying your own vegetation index, create it within a dictionary and pass this to the function:
# As an example, the HUE index from FIELDimageR is given here. In this case, the value of the threshold
# parameter should be set to 0.7 and the comparison parameter should be set to 'greater'
my_index = {'HUE': lambda R,G,B: np.arctan(2*(B-G-R)/30.5*(G-R))}

# Using the function with the ExG index
process_and_save_masks(input_images,
                       masks_out,
                       masked_images_out,
                       threshold=0.1,
                       index_name='ExG',
                       comparison='less',
                       custom_indices=None,
                       output_type='.tif')

# Using the function with the HUE index from FIELDimageR
process_and_save_masks(input_images,
                       masks_out,
                       masked_images_out,
                       threshold=0.7,
                       index_name='HUE',
                       comparison='greater',
                       custom_indices=my_index,
                       output_type='.tif')