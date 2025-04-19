# Use the smallest Linux distribution
FROM ubuntu

RUN cat /proc/cpuinfo

# Update the package list
RUN apt-get update

# Install dependencies for building the application
RUN apt-get install -y \
    build-essential gcc autotools-dev automake autoconf

# Set the working directory
WORKDIR /app

# Copy application files into the container
COPY . /app

# Generate the configure script
RUN autoconf

# Generate the Makefile
RUN automake

# Configure the build environment
RUN ./configure --build=aarch64-unknown-linux-gnu

# Build and install the application
RUN make && make install

# Set the default command to run the application
CMD ["./dwsolver --version"]
