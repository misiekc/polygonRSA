################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../shape/shapes/polygon/Polygon.cpp \
../shape/shapes/polygon/Triangle.cpp 

OBJS += \
./shape/shapes/polygon/Polygon.o \
./shape/shapes/polygon/Triangle.o 

CPP_DEPS += \
./shape/shapes/polygon/Polygon.d \
./shape/shapes/polygon/Triangle.d 


# Each subdirectory must supply rules for building sources it contributes
shape/shapes/polygon/%.o: ../shape/shapes/polygon/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


