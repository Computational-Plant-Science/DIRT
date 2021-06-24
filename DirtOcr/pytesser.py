"""OCR in Python using the Tesseract engine from Google
http://code.google.com/p/pytesser/
by Michael J.T. O'Kelly
V 0.0.1, 3/10/07"""

# from PIL import Image
# import subprocess
# import os
# import util
# import errors
#
# tesseract_exe_name = '/opt/local/bin/tesseract' # Name of executable to be called at command line
# scratch_image_name = "temp.bmp" # This file must be .bmp or other Tesseract-compatible format
# scratch_text_name_root = os.curdir+"/temp" # Leave out the .txt extension
# cleanup_scratch_flag = True  # Temporary files cleaned up after OCR operation
#
# def set_scratch_text_name_root(path,text):
# 	global scratch_image_name
# 	scratch_image_name = text+".bmp" # This file must be .bmp or other Tesseract-compatible format
# 	global scratch_text_name_root
# 	scratch_text_name_root = path # Leave out the .txt extension
#
# def call_tesseract(input_filename, output_filename):
# 	"""Calls external tesseract.exe on input file (restrictions on types),
# 	outputting output_filename+'txt'"""
# 	args = [tesseract_exe_name, input_filename, output_filename]
# 	proc = subprocess.call(args)
# 	#retcode = proc.wait()
# 	#if retcode!=0:
# 		#errors.check_for_errors()
#
# def image_to_string(im, cleanup = cleanup_scratch_flag):
# 	"""Converts im to file, applies tesseract, and fetches resulting text.
# 	If cleanup=True, delete scratch files after operation."""
# 	print( "------------------")
# 	print( scratch_text_name_root)
# 	print(scratch_image_name)
# 	try:
# 		util.image_to_scratch(im, scratch_text_name_root+scratch_image_name)
# 		call_tesseract(scratch_image_name, scratch_text_name_root)
# 		text = util.retrieve_text(scratch_text_name_root)
# 	finally:
# 		if cleanup:
# 			util.perform_cleanup(scratch_image_name, scratch_text_name_root)
# 	return text
#
# def image_file_to_string(filename, cleanup = cleanup_scratch_flag, graceful_errors=True):
# 	"""Applies tesseract to filename; or, if image is incompatible and graceful_errors=True,
# 	converts to compatible format and then applies tesseract.  Fetches resulting text.
# 	If cleanup=True, delete scratch files after operation."""
# 	try:
# 		try:
# 			call_tesseract(filename, scratch_text_name_root)
# 			text = util.retrieve_text(scratch_text_name_root)
# 		except errors.Tesser_General_Exception:
# 			if graceful_errors:
# 				im = Image.open(filename)
# 				text = image_to_string(im, cleanup)
# 			else:
# 				raise
# 	finally:
# 		if cleanup:
# 			util.perform_cleanup(scratch_image_name, scratch_text_name_root)
# 	return text
#
#
# if __name__=='__main__':
# 	im = Image.open('phototest.tif')
# 	text = image_to_string(im)
# 	print(text)
# 	try:
# 		text = image_file_to_string('fnord.tif', graceful_errors=False)
# 	except errors.Tesser_General_Exception as value:
# 		print("fnord.tif is incompatible filetype.  Try graceful_errors=True")
# 		print(value)
# 	text = image_file_to_string('fnord.tif', graceful_errors=True)
# 	print("fnord.tif contents:", text)
# 	text = image_file_to_string('fonts_test.png', graceful_errors=True)
# 	print(text)


