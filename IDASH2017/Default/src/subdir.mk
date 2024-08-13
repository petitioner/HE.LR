################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/MyIDASH2017.cpp \
../src/MyMethods.cpp \
../src/MyTools.cpp 

OBJS += \
./src/MyIDASH2017.o \
./src/MyMethods.o \
./src/MyTools.o 

CPP_DEPS += \
./src/MyIDASH2017.d \
./src/MyMethods.d \
./src/MyTools.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler [CNNinference/Default/src/subdir.mk]'
	g++ -I$(INCLUDE_DIR) -I$(HEAAN_INCLUDE_DIR) -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


