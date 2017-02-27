# this script copies all the image pairs stored in different folders, 
# to a single folder and renames them in the process

# get source filepath from user
filepath=$1

# get destination from user
destination=$2

# get image base name from user
image_base_name=$3_

# # if the destination directory does not exist, create one; it it does, then delete it
# rm -r $destination
# mkdir $destination

echo 'filepath', $filepath
# find the number of folders in this path
num_folders=$(find $filepath/ -maxdepth 1 -type d | wc -l)

# subtract by one because last folder is the processing folder
num_folders=$(($num_folders-1))

# print out the number of folders
echo "number of folders: " $num_folders

# find the names of the folders in the path
folder_names=$(find $filepath/ -maxdepth 1 -type d -printf '%P\n' | sort --version-sort)
echo ${folder_names[*]}

# iterate through the folder names
for name in $folder_names
do
	# print out the folder names
	dir=$filepath/$name
	echo $dir

	counter_var=0

	# set the new filename under which the image will be saved
	# this is the first part of the name
	new_filename_1=$image_base_name
	# echo "new_filename_1", $new_filename_1

	# # this is the second part of the name
	# new_filename_2="$(basename $filepath)"_
	# # echo "new_filename_2", $new_filename_2

	# find the number of files already in the destination
	num_previous_files=$(find $destination/ -name *.tif | wc -l)
	echo "num_previous_files", $num_previous_files
	
	new_filename_3=$((num_previous_files + 1))
	# this is the full name of the file including the path
	# new_filename="$destination/$new_filename_1$new_filename_2$new_filename_3.tif"
	new_filename="$destination/$new_filename_1$_$(printf %04d $new_filename_3).tif"

	# print out new filename
	echo "new_filename", $new_filename

	# copy file to the destination but under the new file name
	cp $(find $dir/ -name *1.tif) $new_filename
done





