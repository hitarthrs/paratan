# Use the OpenMC DAGMC base image
FROM openmc/openmc:develop-dagmc

# ________________________________________ SETUP DAGMC ________________________________________ #
ENV PATH $PATH:/root/DAGMC/bin

# ________________________________________ INSTALL CORE UTILITIES ________________________________________ #
RUN apt-get update && apt-get install -y \
    nano \
    gmsh \
    paraview

# ________________________________________ INSTALL PYTHON LIBRARIES ________________________________________ #
RUN pip install --no-cache-dir \
    pyyaml \
    h5py \
    matplotlib \
    jupyter notebook \
    openmc_plot \
    ipython \
    netCDF4 \
    vtk \
    progress \
    openmc_data_downloader

# ________________________________________ INSTALL ENDF/B-VIII.0 DATA ________________________________________ #
RUN wget -O endf.xz https://anl.box.com/shared/static/uhbxlrx7hvxqw27psymfbhi7bx7s6u6a.xz \
    && tar -xf endf.xz \
    && rm -f endf.xz
RUN echo 'export OPENMC_CROSS_SECTIONS=/endfb-viii.0-hdf5/cross_sections.xml' >> ~/.bashrc

# ________________________________________ INSTALL PYNE WITHOUT SUDO ________________________________________ #
RUN wget https://raw.githubusercontent.com/pyne/install_scripts/main/ubuntu.sh \
    && sed -e 's/\<sudo\>//g' ubuntu.sh > ubuntuNoSudo.sh \
    && rm -f ubuntu.sh \
    && chmod +x ./ubuntuNoSudo.sh \
    && /bin/bash -c './ubuntuNoSudo.sh'

# ________________________________________ FINAL SETUP ________________________________________ #
ENV PATH $PATH:"/root/.local/bin:${PATH}"
