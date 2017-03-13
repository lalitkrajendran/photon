# this is the array of displacements
grad_x_array=(5.00)
echo 'grad_x:' ${grad_x_array[*]}

for grad_x in ${grad_x_array[*]}
do
	echo "grad_x", $grad_x
	top_directory=/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/data/images/focusing/stray-light-study/grad_x=$grad_x
  echo 'top_directory' $top_directory
	./create_sample_processing_directory_02.sh $top_directory

	# find the names of the folders and store them in an array
	folder_names=($(find $top_directory/ -maxdepth 1 -type d -name "perturbation*" -printf '%P\n' | sort --version-sort))

	# this is the destination to copy the files
	destination=$top_directory/processing/reordered-images
	# remove destination folder if it exists
	rm -r $destination
	mkdir -p $destination
	
	for folder_name in ${folder_names[*]}
	do
		echo $folder_name

		# this is the source file path
		source_filepath=$top_directory/$folder_name
		# this is the image base name
		image_base_name=focusing
		./copy_images_to_single_folder_05.sh $source_filepath $destination $image_base_name
	done
done
