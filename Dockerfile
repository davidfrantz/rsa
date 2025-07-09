FROM davidfrantz/base:latest AS builder

# disable interactive frontends
ENV DEBIAN_FRONTEND=noninteractive 

# Environment variables
ENV SOURCE_DIR=$HOME/src/rsa
ENV INSTALL_DIR=$HOME/bin

# Copy src to SOURCE_DIR
RUN mkdir -p $SOURCE_DIR
WORKDIR $SOURCE_DIR
COPY --chown=docker:docker . .

# Build, install
RUN echo "building rsa tools" && \
  make && \
  make install

FROM davidfrantz/base:latest AS final

COPY --chown=docker:docker --from=builder $HOME/bin $HOME/bin

WORKDIR /home/docker
