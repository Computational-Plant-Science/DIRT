from PIL import Image, ImageFile

__all__ = ['fix_orientation']

# PIL's Error "Suspension not allowed here" work around:
# s. http://mail.python.org/pipermail/image-sig/1999-August/000816.html
ImageFile.MAXBLOCK = 1024 * 1024

# The EXIF tag that holds orientation data.
EXIF_ORIENTATION_TAG = 274

# Obviously the only ones to process are 3, 6 and 8.
# All are documented here for thoroughness.
ORIENTATIONS = {
    1: ("Normal", 0),
    2: ("Mirrored left-to-right", 0),
    3: ("Rotated 180 degrees", 180),
    4: ("Mirrored top-to-bottom", 0),
    5: ("Mirrored along top-left diagonal", 0),
    6: ("Rotated 90 degrees", -90),
    7: ("Mirrored along top-right diagonal", 0),
    8: ("Rotated 270 degrees", -270)
}


def fix_orientation(image_path, replace=False):
    # All credits go to Kyle Fox who wrote this EXIF orientation patch.
    # We just modified tiny pieces. https://github.com/kylefox
    image = Image.open(image_path)
    try:
        orientation = image._getexif()[EXIF_ORIENTATION_TAG]
    except (TypeError, AttributeError, KeyError):
        print("WARNING: Image file has no EXIF data")
        orientation = -1
        pass
    if orientation in [3, 6, 8]:
        degrees = ORIENTATIONS[orientation][1]
        image = image.rotate(degrees)
        if replace:
            image.save(image_path, quality=100)
        return image, degrees
    else:
        return image, 0
