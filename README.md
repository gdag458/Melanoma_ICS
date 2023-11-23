# Melanoma_ICS
ICS model code and figure data

The melanoma ICS model is a customized version of SPPARKS (developed and distributed by Sandia National Laboratories). Here one can find the customized cpp/h files to build the ICS model with SPPARKS in linux. To build the ICS SPPARKS executable, first download SPPARKS from https://sjplimp.github.io//download.html. Decompress the "spparks.tar.gz" file resulting in a directory of the form "spparks-DDMon_YY" (eg. spparks-6Sep23). Rename this directory as you like, from here we will refer to it as "ICS_Model". Navigate into the directory to the "src" directory.

Download the ICS_customized_SPPARKS_code directory and place it in the ICS_Model directory. Copy all .cpp and .h into the src directory to put into place the customized code. From within the ICS_Model directory:

cp ICS_customized_SPPARKS_code/app* src/

Make sure you have GCC and OpenMPI modules installed and loaded. To load them, use the command:

module load GCC/7.3.0-2.30 OpenMPI/3.1.1

Navigate into the src directory and build the ICS_model executable using:

cd src; make clean-all; make redsky

Move the executable "spk_redsky", the example initialization file, 16BL_responder_SPPARKS_init, as well as the SPPARKS input file, in.full_trained_model, into a new directory within "ICS_Model" called "simulation_execution":

cd ..; mkdir simulation_execution; cp src/spk_redsky simulation_execution; cp ICS_customized_SPPARKS_code/in.full_trained_model simulation_execution; cp ICS_customized_SPPARKS_code/16BL_responder_SPPARKS_init simulation_execution; cd simulation_execution

Note that the initialization file may be changed out with other slide initialization files which can be found in the data files for Figure S6. When changing the initialization file, make sure to also change its reference in the input file (in.full_trained_model) following SPPARKS syntax. Once in the simulation_execution directory, the ICS model can be run with the command:

mpirun -np 1 ./spk_redsky < in.full_trained_model

