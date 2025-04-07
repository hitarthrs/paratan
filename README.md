# paratan

**Paratan** is a lightweight OpenMC wrapper designed for rapid parametric modeling of tandem mirror fusion devices. It simplifies geometry construction, material assignment, and source specification, enabling easy setup and iteration of neutronics simulations. Built for early-stage design and analysis, **Paratan** helps researchers and engineers explore design spaces efficiently.

**🐳 Setting Up the Environment with Docker**

This project includes a Dockerfile to ensure a consistent and reproducible development environment across different systems.
🔧 Step 1: Build the Docker Image

Run the following command from the root directory of the repository (where the Dockerfile is located):

    sudo docker build -t openmc-paratan .

This builds a Docker image named openmc-paratan using the instructions provided in the Dockerfile.

**🚀 Step 2: Run the Docker Container**

After the image is built, you can launch the container with:

    sudo docker run -v /local_paratan_directory/:/local_workspace -it openmc-paratan

The -v flag mounts the project directory into the container so that changes you make on your host system are reflected inside the container.

The -it flag runs the container in interactive mode with a terminal.

    
