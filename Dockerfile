# Use the smallest Linux distribution
FROM ubuntu

RUN cat /proc/cpuinfo

# Update the package list
RUN apt-get update

# Install dependencies for building the application
RUN apt-get install -y \
    build-essential gcc autotools-dev automake autoconf

# Update config.guess and config.sub to support modern architectures
RUN apt-get install -y wget && \
    wget -O config.guess https://git.savannah.gnu.org/cgit/config.git/plain/config.guess && \
    wget -O config.sub https://git.savannah.gnu.org/cgit/config.git/plain/config.sub && \
    chmod +x config.guess config.sub && \
    mv config.guess config.sub /usr/share/automake-1.16/

# Set the working directory
WORKDIR /app

# Copy application files into the container
COPY . /app

# Update config.guess and config.sub in the project directory
RUN cp /usr/share/automake-1.16/config.guess /app/config.guess && \
    cp /usr/share/automake-1.16/config.sub /app/config.sub

# Generate the configure script
RUN autoreconf

# Generate the Makefile
RUN automake

# Clean up any stale object files before building
RUN make clean || true

# Configure the build environment with an explicit compiler
RUN ./configure --build=$(dpkg-architecture -qDEB_BUILD_GNU_TYPE) CC=gcc

# Build and install the application
RUN make && make install

# Set the default command to run the application
CMD ["./dwsolver --version"]
