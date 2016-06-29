################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CU_SRCS += \
../parallel_ray_tracing.cu 

CU_DEPS += \
./parallel_ray_tracing.d 

OBJS += \
./parallel_ray_tracing.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/opt/cuda/7.0/bin/nvcc -I/home/barracuda/a/lrajendr/usr/include -I/usr/include -I/usr/lib64 -I/usr/local/lib -I../ -I/home/barracuda/a/lrajendr/usr/lib -I/home/barracuda/a/lrajendr/usr/ -G -g -O0 -maxrregcount 32 -Xcompiler -fPIC -gencode arch=compute_20,code=sm_20  -odir "." -M -o "$(@:%.o=%.d)" "$<"
	/opt/cuda/7.0/bin/nvcc -I/home/barracuda/a/lrajendr/usr/include -I/usr/include -I/usr/lib64 -I/usr/local/lib -I../ -I/home/barracuda/a/lrajendr/usr/lib -I/home/barracuda/a/lrajendr/usr/ -G -g -O0 -maxrregcount 32 -Xcompiler -fPIC --compile --relocatable-device-code=false -gencode arch=compute_20,code=compute_20 -gencode arch=compute_20,code=sm_20  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


