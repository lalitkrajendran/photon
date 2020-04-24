################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CU_SRCS += \
../parallel_ray_tracing.cu 

OBJS += \
./parallel_ray_tracing.o 

CU_DEPS += \
./parallel_ray_tracing.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-10.1/bin/nvcc -I/usr/include -I/scratch/shannon/a/lrajendr/usr/include/ -I/scratch/shannon/a/lrajendr/usr/include/teem/ -I/scratch/shannon/c/aether/Projects/BOS/error-analysis/analysis/src/photon/cuda_codes -I/scratch/shannon/c/aether/Projects/BOS/image-generation/analysis/src/cuda-practice/cubic-interpolation/CubicInterpolationCUDA/examples/cuda5_fix -I/scratch/shannon/c/aether/Projects/BOS/image-generation/analysis/src/cuda-practice/cubic-interpolation/CubicInterpolationCUDA/code -I/usr/lib64 -I/usr/local/lib -I../ -G -g -O0 -Xcompiler -fPIC -gencode arch=compute_30,code=sm_30  -odir "." -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-10.1/bin/nvcc -I/usr/include -I/scratch/shannon/a/lrajendr/usr/include/ -I/scratch/shannon/a/lrajendr/usr/include/teem/ -I/scratch/shannon/c/aether/Projects/BOS/error-analysis/analysis/src/photon/cuda_codes -I/scratch/shannon/c/aether/Projects/BOS/image-generation/analysis/src/cuda-practice/cubic-interpolation/CubicInterpolationCUDA/examples/cuda5_fix -I/scratch/shannon/c/aether/Projects/BOS/image-generation/analysis/src/cuda-practice/cubic-interpolation/CubicInterpolationCUDA/code -I/usr/lib64 -I/usr/local/lib -I../ -G -g -O0 -Xcompiler -fPIC --compile --relocatable-device-code=false -gencode arch=compute_30,code=compute_30 -gencode arch=compute_30,code=sm_30  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


