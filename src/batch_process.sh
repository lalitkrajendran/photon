# set parameter
density_grad=12

# generate case name
case_name=150x150-f16-disp$density_grad
echo "case name: " $case_name

# set image directory
image_dir=/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/$case_name
echo "image directory:" $image_dir

# create simulation parameters
echo "******************* creating simulation parameters ************************"
python batch_create_bos_params.py $density_grad

# create density gradient data
echo "******************* creating density gradient data ************************"
python createNRRD.py $density_grad

# make folder to store images
mkdir -p $image_dir

# render images
echo "*************************** generating images *****************************"
python batch_run_bos_simulation.py

# create processing folder
echo "******************** creating processing directory ************************"
cd ./bash-scripts/
./create_sample_processing_directory.sh $case_name

# copy images to a single folder
echo "***************** copying images to a single folder ***********************"
./copy_images_to_single_folder.sh $case_name

# crop images
cd ../
echo "************************** cropping images ********************************"
python crop_images.py $case_name

