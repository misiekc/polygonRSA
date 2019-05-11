################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Config.cpp \
../Main.cpp \
../Packing.cpp \
../PackingGenerator.cpp \
../Parameters.cpp \
../ProgramArguments.cpp \
../RND.cpp \
../Surface.cpp \
../ThreadLocalRND.cpp \
../Timer.cpp \
../Utils.cpp \
../Voxel.cpp \
../VoxelList.cpp 

OBJS += \
./Config.o \
./Main.o \
./Packing.o \
./PackingGenerator.o \
./Parameters.o \
./ProgramArguments.o \
./RND.o \
./Surface.o \
./ThreadLocalRND.o \
./Timer.o \
./Utils.o \
./Voxel.o \
./VoxelList.o 

CPP_DEPS += \
./Config.d \
./Main.d \
./Packing.d \
./PackingGenerator.d \
./Parameters.d \
./ProgramArguments.d \
./RND.d \
./Surface.d \
./ThreadLocalRND.d \
./Timer.d \
./Utils.d \
./Voxel.d \
./VoxelList.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


