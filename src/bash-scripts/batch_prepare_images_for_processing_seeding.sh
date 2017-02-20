# this is the array of displacements
grad_x=(0.50 1.00 1.50 2.00 2.50 3.00 3.50 4.00 4.50 5.00)
echo 'grad_x:' ${grad_x[*]}

# this is the array of seeding densities
seeding_densities=(1 3 5 7 9 11 13 15 17 19)
echo 'seeding_densities: ' ${seeding_densities[*]}

for i in ${grad_x[*]}
do
	for j in ${seeding_densities[*]}
	do
		echo $i $j
		./create_sample_processing_directory_02_seeding.sh $i $j
		./copy_images_to_single_folder_03_seeding.sh $i $j
		python ../post-processing-codes/crop_images.py $i $j
	done
done
