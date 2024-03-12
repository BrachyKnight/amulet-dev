# FROM mcr.microsoft.com/vscode/devcontainers/base:ubuntu-18.04 AS base

# ENV ROOT_PACKAGE_NAME=root_v6.24.00.Linux-ubuntu18-x86_64-gcc7.5.tar.gz
FROM mcr.microsoft.com/vscode/devcontainers/base:ubuntu-20.04 AS base

ENV ROOT_PACKAGE_NAME=root_v6.30.04.Linux-ubuntu20.04-x86_64-gcc9.4.tar.gz
ENV ROOT_PACKAGE_URL=https://root.cern/download/${ROOT_PACKAGE_NAME}
ENV USER=vscode

USER ${USER}
# Needed for initial configuration
RUN mkdir -p /home/vscode/.config/vscode-dev-containers

FROM base AS update

RUN sudo apt update && sudo apt upgrade -y

FROM update AS dependencies
RUN sudo apt install -y dpkg-dev cmake g++ gcc binutils libx11-dev libxpm-dev \
    libxft-dev libxext-dev libssl-dev

FROM dependencies AS optional
RUN sudo apt install -y gfortran libpcre3-dev \
    xlibmesa-glu-dev libglew-dev libftgl-dev \
    libmysqlclient-dev libfftw3-dev libcfitsio-dev \
    graphviz-dev libavahi-compat-libdnssd-dev \
    libldap2-dev libxml2-dev libkrb5-dev \
    libgsl0-dev qtwebengine5-dev

FROM optional AS python
RUN sudo apt install -y python3 python3-pip python3-dev python3-numpy \
    libxml2-dev libxslt1-dev && \
    pip3 install lxml

FROM python AS root
WORKDIR /home/${USER}
RUN wget ${ROOT_PACKAGE_URL} && \
    tar -xzvf ${ROOT_PACKAGE_NAME} && \
    rm -f ${ROOT_PACKAGE_NAME} && \
    echo "alias python='python3'" >> $HOME/.bashrc && \
    echo "PATH=$PATH:~/root/bin" >> $HOME/.bashrc && \
    echo "source ~/root/bin/thisroot.sh" >> $HOME/.bashrc && \
    sudo chown -R ${USER}:${USER} ~/root && \
    sudo chmod -R 755 ~/root


CMD [ "sleep", "infinity" ]