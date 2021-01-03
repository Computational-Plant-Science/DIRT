import os
import subprocess

import PIL as pil
import imageio


def image_to_scratch(im, scratch_image_name):
    # Saves image in memory to scratch file.
    # The .bmp format will be read correctly by Tesseract
    im.save(scratch_image_name, dpi=(300, 300))


def retrieve_text(scratch_text_name_root):
    inf = open(scratch_text_name_root + '.txt')
    text = inf.read()
    inf.close()
    return text


def perform_cleanup(scratch_image_name, scratch_text_name_root):
    # Clean up temporary files from disk
    for name in (scratch_image_name, scratch_text_name_root + '.txt', "tesseract.log"):
        try:
            os.remove(name)
        except OSError:
            pass


tesseract_exe_name = '/opt/local/bin/tesseract'  # Name of executable to be called at command line
scratch_image_name = "temp.bmp"  # This file must be .bmp or other Tesseract-compatible format
scratch_text_name_root = os.curdir + "/temp"  # Leave out the .txt extension
cleanup_scratch_flag = True  # Temporary files cleaned up after OCR operation


def set_scratch_text_name_root(path, text):
    global scratch_image_name
    scratch_image_name = text + ".bmp"  # This file must be .bmp or other Tesseract-compatible format
    global scratch_text_name_root
    scratch_text_name_root = path  # Leave out the .txt extension


def call_tesseract(input_filename, output_filename):
    # Calls external tesseract.exe on input file (restrictions on types),
    # outputting output_filename+'txt'
    args = [tesseract_exe_name, input_filename, output_filename]
    subprocess.call(args)


def image_to_string(im, cleanup=cleanup_scratch_flag):
    # Converts im to file, applies tesseract, and fetches resulting text.
    # If cleanup=True, delete scratch files after operation.
    print("------------------")
    print(scratch_text_name_root)
    print(scratch_image_name)
    try:
        image_to_scratch(im, scratch_text_name_root + scratch_image_name)
        call_tesseract(scratch_image_name, scratch_text_name_root)
        text = retrieve_text(scratch_text_name_root)
    finally:
        if cleanup:
            perform_cleanup(scratch_image_name, scratch_text_name_root)
    return text


def check_for_errors(logfile="tesseract.log"):
    inf = open(logfile)
    text = inf.read()
    inf.close()
    # all error conditions result in "Error" somewhere in logfile
    if text.find("Error") != -1:
        raise RuntimeError(text)


def image_file_to_string(filename, cleanup=cleanup_scratch_flag, graceful_errors=True):
    # Applies tesseract to filename; or, if image is incompatible and graceful_errors=True,
    # converts to compatible format and then applies tesseract.  Fetches resulting text.
    # If cleanup=True, delete scratch files after operation.
    try:
        try:
            call_tesseract(filename, scratch_text_name_root)
            text = retrieve_text(scratch_text_name_root)
        except Exception:
            if graceful_errors:
                im = pil.Image.open(filename)
                text = image_to_string(im, cleanup)
            else:
                raise
    finally:
        if cleanup:
            perform_cleanup(scratch_image_name, scratch_text_name_root)
    return text


def find_text(image, path='temp'):
    imageio.imwrite(path + 'temp.bmp', image)
    set_scratch_text_name_root(path, 'temp.bmp')
    text = image_file_to_string(path + 'temp.bmp', cleanup=cleanup_scratch_flag, graceful_errors=True)
    text = text.translate(None, ",!.;:'{}[]-=()*&^%$#@!~`<>?/|\_+")
    text = ''.join(c for c in text if (c.isalnum() or ' ' or ','))
    text = ' '.join(text.split())
    print(f"Found text: {text}")
    return text


def find_barcode(image, path='temp'):
    imageio.imwrite(path + 'temp.bmp', image)
    args = ['/usr/local/bin/zbarimg', '-q', path + 'temp.bmp']
    code = subprocess.Popen(args, stdout=subprocess.PIPE).communicate()[0]
    print(f"Found barcode: {code}")
    return code
