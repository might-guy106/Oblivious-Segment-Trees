FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt update && apt install -y wget git build-essential net-tools iproute2 iperf iputils-ping libbsd-dev libboost-all-dev numactl time python3 python3-pip
RUN pip3 install pandas matplotlib seaborn rich
WORKDIR /root
COPY . prac
WORKDIR /root/prac
RUN make clean
RUN make -j `nproc`
