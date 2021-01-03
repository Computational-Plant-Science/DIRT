import numpy as np
import cv2


def simple_threshold(image: np.ndarray, threshold: int = 90) -> np.ndarray:
    if threshold < 0 or threshold > 255:
        raise ValueError(f"Threshold must be between 0 and 255")

    _, image = cv2.threshold(image, threshold, 255, cv2.THRESH_BINARY)
    return image


def adaptive_threshold_mean(image: np.ndarray) -> np.ndarray:
    image = cv2.adaptiveThreshold(image, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY, 81, 1)
    return image


def adaptive_threshold_gaussian(image: np.ndarray) -> np.ndarray:
    image = cv2.adaptiveThreshold(image, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY, 81, 1)
    return image


def otsu_threshold(image: np.ndarray) -> np.ndarray:
    image = cv2.createCLAHE(clipLimit=3).apply(image)
    _, image = cv2.threshold(image, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    return image


def binary_threshold(image: np.ndarray):
    if len(np.unique(image)) <= 2:
        print('Binary input detected, skipping thresholding')
        idx1 = np.where(image == np.unique(image)[0])
        idx2 = np.where(image == np.unique(image)[1])
        image[idx1] = False
        image[idx2] = True
    else:
        print('Greyscale input detected, applying binary threshold')
        image = otsu_threshold(image)

    # just a quick fix of the dilation function that caused the binary image to consist of 0 and 2
    idx1 = np.where(image == np.unique(image)[0])
    idx2 = np.where(image == np.unique(image)[1])
    image[idx1] = 0
    image[idx2] = 255

    w, h = np.shape(image)
    image[0, :] = 0
    image[:, 0] = 0
    image[w - 1, :] = 0
    image[:, h - 1] = 0

    return image
