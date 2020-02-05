# Use an official Python runtime as a parent image
FROM rocker/shiny-verse:3.5.0

# Set the working directory to /app
WORKDIR /app

# Copy the current directory contents into the container at /app
COPY . /app

#PROXY settings
ENV http_proxy http://www-cache:3128
ENV https_proxy https://www-cache:3128

# Install any needed packages specified in requirements.txt
RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    pkg-config \
    libnlopt-dev \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libgsl-dev \
    libssl-dev \
    libssh2-1-dev \
    libxml2-dev \
    openssl \
    libpng-dev
    
RUN Rscript installation_script.R

# Make port 80 available to the world outside this container
EXPOSE 80

# Define environment variable
ENV NAME World

# Run app.py when the container launches
#CMD ["Rscript", "runApp.R"]
