################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../surfaces/NBoxFBC.cpp \
../surfaces/NBoxPBC.cpp 

OBJS += \
./surfaces/NBoxFBC.o \
./surfaces/NBoxPBC.o 

CPP_DEPS += \
./surfaces/NBoxFBC.d \
./surfaces/NBoxPBC.d 


# Each subdirectory must supply rules for building sources it contributes
surfaces/%.o: ../surfaces/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


