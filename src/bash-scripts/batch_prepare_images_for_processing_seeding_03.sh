# top_directory=/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/data/images/seeding
top_directory=/home/shannon/c/aether/Projects/BOS/error-analysis/analysis/data/images/seeding
# find the number of folders in this path
num_folders=$(find $top_directory/* -maxdepth 1 -type d -name "delta_x*"| wc -l)
echo "num_folders", $num_folders

num_folders=15

# find the names of the folders and store them in an array
folder_names=($(find $top_directory/ -maxdepth 1 -type d -name "delta_x*" -printf '%P\n' | sort --version-sort))
echo ${folder_names[*]}

for i in $(seq 10 $((num_folders - 1)))
do
	echo 'i=' $i
	folder_name=${folder_names[i]}
	echo 'folder name', $folder_name

	sub_folder_names=$(find $top_directory/$folder_name/ -maxdepth 1 -type d -printf '%P\n' | sort --version-sort)

	echo 'subfolder names', ${sub_folder_names[*]}

	for sub_folder_name in ${sub_folder_names[*]}
	do
		# create a processing directory
		echo 'sub folder_name', $sub_folder_name
		./create_sample_processing_directory_02.sh $top_directory/$folder_name/$sub_folder_name

		# copy images to single folder
		# this is the source file path
		source_filepath=$top_directory/$folder_name/$sub_folder_name
		# this is the destination to copy the files
		destination=$top_directory/$folder_name/$sub_folder_name/processing/reordered-images
		# this is the image base name
		image_base_name=seeding
		./copy_images_to_single_folder_04.sh $source_filepath $destination $image_base_name

		# this is the source file path
		source_filepath=$top_directory/$folder_name/$sub_folder_name/processing/reordered-images/
		# this is the destination where the cropped images will be saved
		destination=$top_directory/$folder_name/$sub_folder_name/processing/cropped-images/
		# these are the min and max coordiantes to crop
		xmin=256
		xmax=768
		ymin=256
		ymax=768
		python ../pre-processing-codes/crop_images_02.py $source_filepath $destination $xmin $xmax $ymin $ymax

	done
done
