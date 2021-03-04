################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include subdir.mk
-include .metadata/.plugins/org.eclipse.cdt.make.core/subdir.mk

OS_SUFFIX := $(subst Linux,linux,$(subst Darwin/x86_64,darwin,$(shell uname -s)/$(shell uname -m)))

-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(CU_DEPS)),)
-include $(CU_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

-include ../makefile.defs

# OpenGL+GLUT OS-specific define
ifeq ($(shell uname -s),Darwin)
GLUT_LIBS := -Xlinker -framework -Xlinker GLUT -Xlinker -framework -Xlinker OpenGL
else
GLUT_LIBS := -lGL -lGLU -lglut 
endif


# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: libparallel_ray_tracing.so

# Tool invocations
libparallel_ray_tracing.so: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: NVCC Linker'
	# /usr/local/cuda-10.1/bin/nvcc --cudart static -L/usr/include/ -L/scratch/shannon/a/lrajendr/usr/lib/ -L/scratch/shannon/a/lrajendr/usr/bin/ -L/scratch/shannon/a/lrajendr/usr/lib64 -L/scratch/shannon/a/lrajendr/usr/include/teem/ -L/scratch/shannon/a/lrajendr/usr/include/ -L/usr/local/include/ -L/usr/lib64/ -shared --relocatable-device-code=false -gencode arch=compute_30,code=compute_30 -gencode arch=compute_30,code=sm_30 -link -o  "libparallel_ray_tracing.so" $(OBJS) $(USER_OBJS) $(LIBS)
	/usr/local/cuda-10.1/bin/nvcc --cudart static -L/usr/include/ -L/scratch/shannon/a/lrajendr/usr/lib/ -L/scratch/shannon/a/lrajendr/usr/bin/ -L/scratch/shannon/a/lrajendr/usr/lib64 -L teem/ -L/scratch/shannon/a/lrajendr/usr/include/ -L/usr/local/include/ -L/usr/lib64/ -shared --relocatable-device-code=false -gencode arch=compute_30,code=compute_30 -gencode arch=compute_30,code=sm_30 -link -o  "libparallel_ray_tracing.so" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(LIBRARIES)$(CC_DEPS)$(C++_DEPS)$(C_UPPER_DEPS)$(CXX_DEPS)$(OBJS)$(CU_DEPS)$(CPP_DEPS)$(C_DEPS) libparallel_ray_tracing.so
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets
