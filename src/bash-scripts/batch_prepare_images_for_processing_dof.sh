# this is the array of displacements
grad_x_array=(5.00)
echo 'grad_x:' ${grad_x_array[*]}

# this is the array of f numbers
f_number_array=(02 04 05 08 11 16 22 32 64)
echo "f_numbers:" ${f_number_array[*]}

for grad_x in ${grad_x_array[*]}
do
	for f_number in ${f_number_array[*]}
	do
		echo $grad_x $f_number

		# this is the top write directory
		top_write_directory=/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/data/images/dof/grad_x=$grad_x/f_number=$f_number
		
		./create_sample_processing_directory_02.sh $top_write_directory
		
		# this is the source file path
		source_filepath=$top_write_directory
		# this is the destination to copy the files
		destination=$top_write_directory/processing/reordered-images
		# this is the image base name
		image_base_name=dof
		./copy_images_to_single_folder_04.sh $source_filepath $destination $image_base_name
		
		# this is the source file path
		source_filepath=$top_write_directory/processing/reordered-images/
		# this is the destination where the cropped images will be saved
		destination=$top_write_directory/processing/cropped-images/
		# these are the min and max coordiantes to crop
		xmin=256
		xmax=768
		ymin=256
		ymax=768
		python ../pre-processing-codes/crop_images_02.py $source_filepath $destination $xmin $xmax $ymin $ymax


	done
done
