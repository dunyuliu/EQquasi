# Use the latest official Ubuntu Docker image
FROM ubuntu:latest

# Install OpenMPI
RUN apt-get update && apt-get install -y openmpi-bin

# Add a new user 'eq'
RUN useradd -ms /bin/bash eq

# Set 'eq' as the default user
USER eq

# Set the working directory
WORKDIR /home/eq
