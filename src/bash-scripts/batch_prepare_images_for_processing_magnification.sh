# this is the array of displacements
grad_x_array=(0.50 1.00 1.50 2.00 2.50 3.00 3.50 4.00 4.50 5.00)
#grad_x_array=(5.00)

echo 'grad_x:' ${grad_x_array[*]}

# this is the array of focal lengths
magnification_array=(100 125 150 175 200 225 250 275 300)
echo 'magnification: ' ${magnification_array[*]}

for grad_x in ${grad_x_array[*]}
do
	for magnification in ${magnification_array[*]}
	do
		echo $grad_x $magnification
		# this is the top write directory
		top_write_directory=/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/magnification/grad_x=$grad_x/$magnification
		./create_sample_processing_directory_02.sh $top_write_directory
		
		# this is the source file path
		source_filepath=$top_write_directory
		# this is the destination to copy the files
		destination=$top_write_directory/processing/reordered-images
		# this is the image base name
		image_base_name=magnification
		./copy_images_to_single_folder_04.sh $source_filepath $destination $image_base_name
		
		# this is the source file path
		source_filepath=$top_write_directory/processing/reordered-images/
		# this is the destination where the cropped images will be saved
		destination=$top_write_directory/processing/cropped-images/
		# these are the min and max coordiantes to crop
		xmin=416
		xmax=608
		ymin=416
		ymax=608
		python ../pre-processing-codes/crop_images_02.py $source_filepath $destination $xmin $xmax $ymin $ymax
	done
done
