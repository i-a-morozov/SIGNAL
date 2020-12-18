# docker build -t signal .
# docker run --rm -it --net=host --name signal -e DISPLAY -v ~/.Xauthority:/root/.Xauthority --entrypoint=/bin/bash signal

FROM ubuntu

ENV DEBIAN_FRONTEND noninteractive

RUN apt update
RUN apt -y upgrade

RUN apt -y install gcc gcc-10
RUN apt -y install gfortran gfortran-10
RUN apt -y install g++
RUN apt -y install make
RUN apt -y install cmake
RUN apt -y install libatlas-base-dev
RUN apt -y install libomp-dev
RUN apt -y install git
RUN apt -y install libfftw3-dev
RUN apt -y install libgsl-dev

RUN git clone https://github.com/Reference-LAPACK/lapack.git
RUN cd lapack && mkdir build && cd build && cmake -DCMAKE_INSTALL_LIBDIR=$HOME/.local/lapack .. && cmake --build . -j --target install

RUN git clone https://github.com/opencollab/arpack-ng.git
RUN cd arpack-ng && mkdir build && cd build && cmake -DCMAKE_INSTALL_LIBDIR=$HOME/.local/arpack .. -DBUILD_SHARED_LIBS=OFF && make && make install

RUN git clone https://github.com/epics-base/epics-base.git
RUN cd epics-base && make

ENV EPICS_BASE /epics-base
ENV EPICS_HOST_ARCH linux-x86_64
ENV PATH ${EPICS_BASE}/bin/${EPICS_HOST_ARCH}:${PATH}

RUN mkdir signal
COPY src /signal/src
COPY examples /signal/examples
COPY epics /signal/epics
RUN cd /signal/src && make
RUN cd $HOME/.local && mkdir signal && cd signal && cp -rf /signal/src/libsignal.a . && cp -rf /signal/src/mod .

RUN cp -rf $HOME/.local/lapack/* /usr/lib
RUN cp -rf $HOME/.local/arpack/* /usr/lib
RUN cp -rf $HOME/.local/signal/* /usr/lib

ENV DISPLAY :0
ENV NO_AT_BRIDGE 1