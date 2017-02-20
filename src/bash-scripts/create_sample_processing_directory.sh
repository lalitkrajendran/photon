# this script creates all the subfolders for a processing case

# specify top directory
top_write_directory=/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/seeding/processing
#top_write_directory=/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing


# set case name
#case_name=50x50-f16-grad_x3.0
case_name=$1

# create new directory with case name
destination=$top_write_directory/$case_name
echo "creating a folder at " $destination

mkdir -p $destination

# create subfolders within the the new directory

# create directory to store images
mkdir -p $destination/reordered-images

# create directory to store cropped images
mkdir -p $destination/cropped-images

# create directory to store plots
mkdir -p $destination/plots

# create directory to store results
mkdir -p $destination/results

# create subfolder to store particle-id
mkdir -p $destination/results/particle-id

# create subfolder to store particle-size
mkdir -p $destination/results/particle-size

# create subfolder to store vectors
mkdir -p $destination/results/vectors
