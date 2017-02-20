# this script deletes old vectors


# this is the array of displacements
grad_x_array=(0.50 1.00 1.50 2.00 2.50 3.00 3.50 4.00 4.50 5.00)

echo 'grad_x:' ${grad_x_array[*]}

# this is the array of seeding densities
seeding_densities_array=(2 4 6 8 10 12 14 16 18 20)
echo 'seeding_densities: ' ${seeding_densities_array[*]}

for grad_x in ${grad_x_array[*]}
do
	for density in ${seeding_densities_array[*]}
	do
		echo $grad_x $density
		# this is the source filepath
		filepath=/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/seeding/grad_x=$grad_x/$density

		# show files before deleting
		ls $filepath/processing/results/vectors/*.mat
		
		# remove the mat files in the vectors folder
		rm $filepath/processing/results/vectors/*.mat

		# show files after deleting
		ls $filepath/processing/results/vectors/*.mat
				
	done
done

