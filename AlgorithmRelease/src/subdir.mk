################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/binary_functions.cpp \
../src/binary_regression.cpp \
../src/projections.cpp 

OBJS += \
./src/binary_functions.o \
./src/binary_regression.o \
./src/projections.o 

CPP_DEPS += \
./src/binary_functions.d \
./src/binary_regression.d \
./src/projections.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++  -std=c++0x -I/Users/yaara/Documents/cppworkspace/projections/lib/sparsehash/include -O3 -Wall -c -fmessage-length=0 -D NDEBUG -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


