FROM tiagopeixoto/graph-tool

LABEL maintainer="Wes Bonelli"

# RUN apt-get update && \
#     echo 'deb http://downloads.skewed.de/apt/xenial xenial universe' | tee -a  /etc/apt/sources.list && \
#     echo 'deb-src http://downloads.skewed.de/apt/xenial xenial universe' | tee -a  /etc/apt/sources.list && \
COPY . /opt/DIRT

RUN pacman -S --noconfirm gcc git python-pip && \
    cd /opt/DIRT && \
    sed -i 's#/usr/local/bin/zbarimg#/usr/bin/zbarimg#' /opt/DIRT/DirtOcr/__init__.py && \
    pip install --upgrade pip && \
    pip install -r /opt/DIRT/requirements.txt

ENV LC_ALL=C
ENV DISPLAY=:1

CMD python /opt/DIRT/main.py "$@"
