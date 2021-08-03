FROM python:3.8-slim-buster as build

# Install unzip and texlive. The latter one is needed for Python PyX.

RUN apt-get update
RUN apt-get install -y build-essential
RUN apt-get install -y unzip
RUN apt-get clean && rm -rf /var/lib/apt/lists/*

# Prepare LibGAL distribution for build

RUN mkdir /opt/src
COPY ./LibGAL/libgal-0.6.0.zip /opt/src
WORKDIR /opt/src
RUN unzip -q ./libgal-0.6.0.zip

# Build LibGAL library

WORKDIR /opt/src/libgal-0.6.0
RUN ./configure --prefix=/opt/libgal-0.6.0
RUN make clean
RUN make
RUN make install

# Make psgp4gm module written in Cython

RUN pip install cython

RUN mkdir /opt/satpass
RUN mkdir /opt/satpass/psgp4gm
COPY ./psgp4gm/psgp4gm.pyx /opt/satpass/psgp4gm
COPY ./psgp4gm/build2.sh /opt/satpass/psgp4gm
WORKDIR /opt/satpass/psgp4gm
RUN chmod +x ./build2.sh
RUN ./build2.sh

# Start new image

FROM python:3.8-slim-buster as satpass

COPY --from=build /opt/libgal-0.6.0/ /opt/libgal-0.6.0 

# Create an appropriate directory structure and copy over SatPass scripts.

RUN mkdir /opt/satpass
RUN mkdir /opt/satpass/psgp4gm
RUN mkdir /opt/satpass/passes
RUN mkdir /opt/satpass/passes/data
RUN mkdir /opt/satpass/satools
COPY ./psgp4gm/psgp4gm.pyx /opt/satpass/psgp4gm
COPY ./psgp4gm/build2.sh /opt/satpass/psgp4gm
RUN chmod +x opt/satpass/psgp4gm/build2.sh
COPY --from=build /opt/satpass/psgp4gm/psgp4gm.so /opt/satpass/psgp4gm
COPY ./passes/*.py /opt/satpass/passes
COPY ./passes/data/satellites.conf /opt/satpass/passes/data
COPY ./satools/*.py /opt/satpass/satools

# Install TeX needed by the Python PyX module

RUN apt-get update
RUN apt-get install -y texlive
RUN apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Python modules used by SatPass

RUN /usr/local/bin/python -m pip install --upgrade pip
RUN pip install pyyaml
RUN pip install geographiclib
RUN pip install pytz
RUN pip install timezonefinder
RUN pip install pyx

# Some final touches ...

WORKDIR /opt/satpass/passes
RUN echo "\e[B": history-search-forward > ~/.inputrc
RUN echo "\e[A": history-search-backward >> ~/.inputrc
RUN echo "PS1='\w\$ '" >> ~/.profile
RUN echo "export INPUTRC=~/.inputrc" >> ~/.profile
RUN echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/libgal-0.6.0/lib" >> ~/.bashrc
