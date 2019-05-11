################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../shape/Shape.cpp \
../shape/ShapeFactory.cpp 

OBJS += \
./shape/Shape.o \
./shape/ShapeFactory.o 

CPP_DEPS += \
./shape/Shape.d \
./shape/ShapeFactory.d 


# Each subdirectory must supply rules for building sources it contributes
shape/%.o: ../shape/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


