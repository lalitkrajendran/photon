# this script copies all the image pairs stored in different folders, 
# to a single folder and renames them in the process

# # this is the grad_x value
# grad_x=$1
#
# # this is the seeding density
# density=$2

# this is the array of displacements
grad_x_array=(0.50 1.00 1.50 2.00 2.50 3.00 3.50 4.00 4.50 5.00)
echo 'grad_x:' ${grad_x[*]}

# this is the array of seeding densities
seeding_densities_array=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20)
echo 'seeding_densities: ' ${seeding_densities[*]}

for grad_x in ${grad_x_array[*]}
do
	for density in ${seeding_densities_array[*]}
	do
		echo $grad_x $density

		# this is the top write directory
		top_write_directory=/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/seeding/grad_x=$grad_x/$density
		
		# this is the source file path
		filepath=$top_write_directory
		
		echo 'filepath', $filepath


		# find the number of folders in this path
		num_folders=$(find $filepath/* -maxdepth 0 -type d | wc -l)

		# subtract by one because last folder is the processing folder
		num_folders=$(($num_folders-1))

		# print out the number of folders
		echo "number of folders: " $num_folders

		# find the names of the folders in the path
		# folder_names=$(find $filepath/* -maxdepth 0 -type d | sort --version-sort)
		folder_names=$(seq 1 $num_folders)
		echo $folder_names

		# iterate through the folder names
		# for dir in $folder_names
		for name in $folder_names
		do
		  # print out the folder names
		  dir=$filepath/$name
		  echo $dir
		  
		  # delete light ray positions and directions
		  rm -rf $dir/light-ray-positions
		  rm -rf $dir/light-ray-directions
		done
	done
done








