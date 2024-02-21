# Prepare the docker environment for EQquasi
# Use the latest official Ubuntu Docker image as the base
FROM ubuntu:latest

# Install OpenMPI and MUMPS
#RUN apt-get update && apt-get install -y \
#    openmpi-bin \
#    libmumps-dev

# Add a new user 'eq'
RUN useradd -ms /bin/bash eq

# Set 'eq' as the default user
USER eq

# Set the working directory
WORKDIR /home/eq
