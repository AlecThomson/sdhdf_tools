FROM ubuntu:18.04 AS spython-base
ENV SDHDF_RUNTIME=/usr/include/runtime
ENV PATH="/Processing_Software/bin:$PATH"
RUN echo "Updating apt repositories"
RUN apt update && apt upgrade -y
RUN DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata
RUN apt install -y  \
    git wget \
    liberfa1 libhdf5-100 libhdf5-dev liberfa-dev \
    pgplot5 \
    tcsh \
    autoconf \
    libtool \
    libcfitsio5 \
    libcfitsio-dev \
    screen \
    gnuplot \
    libx11-dev \
    make \
    cmake \
    build-essential \
    pkg-config \
    casacore-data casacore-dev \
    libblas-dev liblapack-dev \
    python3 \
    libboost-date-time-dev \
    libboost-filesystem-dev \
    libboost-system-dev \
    libboost-test-dev \
    libcfitsio-dev \
    libfftw3-dev \
    libgsl-dev \
    libgtkmm-3.0-dev \
    liblua5.3-dev \
    libpng-dev  \
    libpython3-dev \
    libxml2-dev \
    apt autoremove -y
RUN apt clean -y
RUN echo "Installing calceph"
RUN cd / \
    && wget https://www.imcce.fr/content/medias/recherche/equipes/asd/calceph/calceph-3.4.7.tar.gz \
    && tar -zxvf calceph-3.4.7.tar.gz  \
    && cd calceph-3.4.7/ \
    && apt install -y gfortran \
    && ./configure F77=gfortran --prefix=/usr/include/calceph-3.4.7 \
    && make \
    && make install \
    && cd / \
    && rm -rf calceph-3.4.7.tar.gz calceph-3.4.7/
RUN echo "Installing sdhdfProc"
RUN cd / \
    && git clone https://bitbucket.csiro.au/scm/cpda/sdhdf_tools.git \
    && cd sdhdf_tools/c/sdhdfProc \
    && ./bootstrap  \
    && ./configure --prefix=/ CFLAGS="-I/usr/include/calceph-3.4.7/include -I/usr/include/hdf5/serial/" LDFLAGS="-L/usr/include/calceph-3.4.7/lib -L/usr/lib/x86_64-linux-gnu/hdf5/serial"  \
    && make \
    && make install \
    && cd / \
    && mv sdhdf_tools/c/sdhdfProc/runtime /usr/include/ \
    && rm -rf sdhdf_tools/
RUN echo "Installing aoflagger"
RUN cd / \
    && git clone https://gitlab.com/aroffringa/aoflagger.git \
    && cd aoflagger/ \
    && mkdir build \
    && cd build/ \
    && cmake -DPORTABLE=ON ..  \
    && make -j 8 \
    && make install \
    && cd / \
    && rm -rf aoflagger/
CMD exec "$@"