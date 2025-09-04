import os
import cv2
import numpy as np

# Many of these functions are adapted from a script by Lei Mao, University of Chicago, 3/1/2018

def inside_rect(rect, num_cols, num_rows):
    """
    Determine if the four corners of a rotated rectangle (rect) lie
    fully within the image.
    """
    rect_center = rect[0]
    rect_center_x = rect_center[0]
    rect_center_y = rect_center[1]

    if (rect_center_x < 0) or (rect_center_x > num_cols):
        return False
    if (rect_center_y < 0) or (rect_center_y > num_rows):
        return False

    box = cv2.boxPoints(rect)
    x_max = int(np.max(box[:, 0]))
    x_min = int(np.min(box[:, 0]))
    y_max = int(np.max(box[:, 1]))
    y_min = int(np.min(box[:, 1]))

    if (x_max <= num_cols) and (x_min >= 0) and (y_max <= num_rows) and (y_min >= 0):
        return True
    else:
        return False


def rect_bbx(rect):
    """
    Return the bounding box (as an upright rect) of a rotated rectangle (rect).
    """
    box = cv2.boxPoints(rect)
    x_max = int(np.max(box[:, 0]))
    x_min = int(np.min(box[:, 0]))
    y_max = int(np.max(box[:, 1]))
    y_min = int(np.min(box[:, 1]))

    center = (int((x_min + x_max) // 2), int((y_min + y_max) // 2))
    width = int(x_max - x_min)
    height = int(y_max - y_min)
    angle = 0

    return (center, (width, height), angle)


def image_rotate_without_crop(mat, angle):
    """
    Rotate an image by 'angle' degrees without clipping its corners.
    """
    height, width = mat.shape[:2]
    image_center = (width / 2, height / 2)

    rotation_mat = cv2.getRotationMatrix2D(image_center, angle, 1.0)
    abs_cos = abs(rotation_mat[0, 0])
    abs_sin = abs(rotation_mat[0, 1])

    bound_w = int(height * abs_sin + width * abs_cos)
    bound_h = int(height * abs_cos + width * abs_sin)

    rotation_mat[0, 2] += bound_w / 2 - image_center[0]
    rotation_mat[1, 2] += bound_h / 2 - image_center[1]

    rotated_mat = cv2.warpAffine(mat, rotation_mat, (bound_w, bound_h))
    return rotated_mat


def crop_rectangle(image, rect):
    """
    Crop an upright rectangle from the image, given by 'rect'.
    """
    num_rows = image.shape[0]
    num_cols = image.shape[1]

    if not inside_rect(rect=rect, num_cols=num_cols, num_rows=num_rows):
        print("Proposed rectangle is not fully in the image.")
        return None

    rect_center_x = rect[0][0]
    rect_center_y = rect[0][1]
    rect_width = rect[1][0]
    rect_height = rect[1][1]

    return image[
        rect_center_y - rect_height // 2 : rect_center_y + rect_height - rect_height // 2,
        rect_center_x - rect_width // 2  : rect_center_x + rect_width - rect_width // 2
    ]


def crop_rotated_rectangle(image, rect):
    """
    Crop a rotated rectangle from the image, given by 'rect'.
    """
    num_rows = image.shape[0]
    num_cols = image.shape[1]

    if not inside_rect(rect=rect, num_cols=num_cols, num_rows=num_rows):
        print("Proposed rectangle is not fully in the image.")
        return None

    rotated_angle = rect[2]
    rect_bbx_upright = rect_bbx(rect=rect)

    # Crop the upright bounding box of the rotated rect
    rect_bbx_upright_image = crop_rectangle(image=image, rect=rect_bbx_upright)
    if rect_bbx_upright_image is None:
        return None

    # Rotate that bounding box subimage by the rotated rect angle
    rotated_rect_bbx_upright_image = image_rotate_without_crop(
        mat=rect_bbx_upright_image, angle=rotated_angle
    )

    rect_width = rect[1][0]
    rect_height = rect[1][1]

    crop_center = (
        rotated_rect_bbx_upright_image.shape[1] // 2,
        rotated_rect_bbx_upright_image.shape[0] // 2
    )

    return rotated_rect_bbx_upright_image[
        int(crop_center[1] - rect_height // 2) : int(crop_center[1] + (rect_height - rect_height // 2)),
        int(crop_center[0] - rect_width // 2)  : int(crop_center[0] + (rect_width - rect_width // 2))
    ]


def get_contour(img):
    """
    Find and return the sorted contours of an image (largest first).
    """
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    _, thres = cv2.threshold(gray, 1, 255, cv2.THRESH_BINARY)
    contours, _ = cv2.findContours(thres, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    contours = sorted(contours, key=cv2.contourArea, reverse=True)
    return contours


def get_cropped_rotated(img, rect):
    """
    Helper function to rotate and crop an image around a given rotated rect.
    """
    mult = 1.0
    box = cv2.boxPoints(rect)
    box = np.int32(box)

    Xs = [i[0] for i in box]
    Ys = [i[1] for i in box]
    x1 = min(Xs)
    x2 = max(Xs)
    y1 = min(Ys)
    y2 = max(Ys)

    angle = rect[2]
    # Adjust angle if height is actually larger than width
    if y2 - y1 < x2 - x1:
        angle -= 90

    size = (int(mult * (x2 - x1)), int(mult * (y2 - y1)))
    M = cv2.getRotationMatrix2D((size[0] / 2, size[1] / 2), angle, 1.0)
    cropped = cv2.warpAffine(img, M, size)
    return cropped


def crop_image_region_of_interest(img):
    """
    Given an image, add a small border, find the largest contour,
    compute the minimum bounding rectangle, and crop to that region.
    """
    # Add a small black border so contours go to the edge if needed
    img = cv2.copyMakeBorder(img, 1, 1, 1, 1, cv2.BORDER_CONSTANT, value=0)
    cnt = get_contour(img)
    if not cnt:
        return img  # No contours found, just return original

    cnt = cnt[0]  # Largest contour
    min_rect = cv2.minAreaRect(cnt)
    
    # Crop the rotated rectangle
    cropped = crop_rotated_rectangle(img, min_rect)
    return cropped if cropped is not None else img


def iter_images_from_folder(folder, valid_exts=('.png', '.jpg', '.jpeg',
                                               '.tif', '.tiff')):
    """
    Yield (image, filename_without_ext) pairs one at a time.
    Reading happens on‑the‑fly so only one image is ever in memory.
    """
    for entry in os.scandir(folder):             # faster than os.listdir＋join
        if entry.is_file() and entry.name.lower().endswith(valid_exts):
            img = cv2.imread(entry.path, cv2.IMREAD_COLOR)
            if img is not None:
                yield img, os.path.splitext(entry.name)[0]
            else:
                print(f"Could not read {entry.path}")


def crop_ROI(
    images_to_crop: str, 
    clean_images_out: str, 
    output_file_type: str = '.tif'
):
    """
    Accepts images that contain a rectangle or square with a region of interest but that 
    also contain excess black background, usually because these images were cropped from a GeoTIFF
    and they retained their spatial orientation.
      1. Reads all images in 'images_to_crop' folder.
      2. For each image, crops out the black background area.
      3. Saves the cleaned images to 'clean_images_out' with the specified 'output_type'.
      
    Parameters
    ----------
    images_to_crop : str [required]
        Path to the images that contain a region of interest within a black background that
        needs to be removed.
    clean_images_out : str [required]
        Path to the folder where cropped images will be saved.
    output_file_type: str [optional]
        Type of image that will be saved. Usually this will be .tif or .jpg. Defaults to .jpg.
        
    Returns
    -------
    None

    """

    # Loop over images, crop, and save (after making the output directory if it doesn't exist)
    os.makedirs(clean_images_out, exist_ok=True)
    
    for img, stem in iter_images_from_folder(images_to_crop):
        cleaned = crop_image_region_of_interest(img)
        save_path = os.path.join(clean_images_out, f"{stem}{output_file_type}")
        cv2.imwrite(save_path, cleaned)

    print(f"All cleaned images have been saved to {clean_images_out}")


# Path to cropped images (most of the time these will contain residual black background that we want to remove)
images_to_crop = 'path/to/cropped/images'
# Path where you'd like the cleaned images to be saved
clean_images_out = 'path/to/cleaned/cropped/images'
# I highly recommend working with TIFs, not JPGs. For downstream masking and vegetation index steps, JPGs have been known to show
# pixel blending, where a pixel that's supposed to be completely black (meaning values of 0 for the R, G, and B bands) will have small
# nonzero values. Since the scripts following this one use all-0 pixels as a way to determine which ones should be excluded from computations,
# this can seriously impact results.
output_file_type = '.tif'

crop_ROI(
    images_to_crop=images_to_crop,
    clean_images_out=clean_images_out,
    output_file_type=output_file_type)
