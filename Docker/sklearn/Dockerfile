FROM ubuntu:22.04 AS system


# system container dependencies
RUN apt-get update 
RUN DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y install tzdata
RUN apt-get install -q -y --no-install-recommends \
    libglu1-mesa-dev \
    libgomp1 \
    libopenslide-dev \
    python3 \
    python3-pip \
    python3-venv \
    git \
    unzip
RUN apt-get clean 
RUN rm -rf /var/cache/apt/archives /var/lib/apt/lists/*

# Manage to use the virtual environment
RUN python3 -m venv /.venv
ENV VIRTUAL_ENV="/.venv"
RUN export VIRTUAL_ENV
ENV _OLD_VIRTUAL_PATH="$PATH"
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
RUN export PATH

# Install python dependencies
RUN pip install SimpleITK numpy scipy scikit-learn scikit-image
RUN pip install m2aia

# provide app data
FROM system AS application
COPY _sklearn.py _sklearn.py
ENTRYPOINT [ "python", "_sklearn.py" ]
