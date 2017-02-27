# this is the array of displacements
grad_x_array=(0.50 1.00 1.50 2.00 2.50 3.00 3.50 4.00 4.50 5.00)
echo 'grad_x:' ${grad_x_array[*]}

# # this is the array of seeding densities
# seeding_densities=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20)
# echo 'seeding_densities: ' ${seeding_densities[*]}

# top read directory
top_write_directory=/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/data/images/dotsize


for grad_x in ${grad_x_array[*]}
do

	# get list of folders for this displacement
	folder_names=$(find $top_write_directory/grad_x=$grad_x/ -maxdepth 1 -type d -printf '%P\n' | sort --version-sort)
	echo ${folder_names[*]}
	
	for folder in ${folder_names[*]}
	do
		echo $grad_x $folder

		./create_sample_processing_directory_02.sh $top_write_directory/grad_x=$grad_x/$folder

		# this is the source file path
		source_filepath=$top_write_directory/grad_x=$grad_x/$folder
		# this is the destination to copy the files
		destination=$top_write_directory/grad_x=$grad_x/$folder/processing/reordered-images
		# this is the image base name
		image_base_name=dotsize
		./copy_images_to_single_folder_04.sh $source_filepath $destination $image_base_name

		# this is the source file path
		source_filepath=$top_write_directory/grad_x=$grad_x/$folder/processing/reordered-images/
		# this is the destination where the cropped images will be saved
		destination=$top_write_directory/grad_x=$grad_x/$folder/processing/cropped-images/
		# these are the min and max coordiantes to crop
		xmin=256
		xmax=768
		ymin=256
		ymax=768
		python ../pre-processing-codes/crop_images_02.py $source_filepath $destination $xmin $xmax $ymin $ymax

	done
done
