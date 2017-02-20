# this script deletes old vectors

# this is the grad_x value
grad_x=$1

# this is the seeding density
density=$2

# this is the base name for the images
image_base_name=seeding_
# this is the destination where the images will be copied to
destination=/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/seeding/grad_x=$grad_x/$density/processing/reordered-images

# if the destination directory does not exist, create one; it it does, then delete it
rm -r $destination
mkdir $destination

# this is the source filepath
filepath=/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/seeding/grad_x=$grad_x/$density

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

