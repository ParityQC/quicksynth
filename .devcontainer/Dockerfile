# Copyright (c) 2026 Parity Quantum Computing GmbH
# SPDX-License-Identifier: MIT

# Shared dev image for quicksynth: Python 3.12, uv and a C++ toolchain. It is
# built as quicksynth-dev-base; each developer's .devcontainer/local/Dockerfile
# adds a personal layer on top (README.md in this folder).
#
# The official python image is used rather than mcr.microsoft.com/devcontainers
# because the latter ships a yarn apt repository whose signing key has rotated,
# which breaks `apt-get update` inside the image.
FROM python:3.12-bookworm

ARG DEBIAN_FRONTEND=noninteractive
ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=1000

# C++ toolchain: build, debug, format.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        g++ \
        clang-format \
        cmake \
        ninja-build \
        gdb \
        git \
        curl \
        less \
        sudo \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Non-root user with passwordless sudo, matching the usual devcontainer setup.
# On Linux hosts the dev container tooling changes its UID to the host user's.
RUN groupadd --gid ${USER_GID} ${USERNAME} \
    && useradd --uid ${USER_UID} --gid ${USER_GID} -m -s /bin/bash ${USERNAME} \
    && echo "${USERNAME} ALL=(root) NOPASSWD:ALL" > /etc/sudoers.d/${USERNAME} \
    && chmod 0440 /etc/sudoers.d/${USERNAME}

# uv drives the Python environment, at the version the lock file was written with.
COPY --from=ghcr.io/astral-sh/uv:0.12.17 /uv /uvx /usr/local/bin/

# The Python environment is per checkout (UV_PROJECT_ENVIRONMENT, set in
# devcontainer.json), so it is not created here.
ENV UV_LINK_MODE=copy \
    UV_COMPILE_BYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    CMAKE_GENERATOR=Ninja

RUN mkdir -p /home/${USERNAME}/.cache/uv \
    && chown -R ${USER_UID}:${USER_GID} /home/${USERNAME}/.cache

# A login shell re-reads /etc/profile, which resets PATH; put the checkout's
# environment (DEV_VENV, set in devcontainer.json) back in front.
RUN printf '%s\n' \
        'if [ -n "${DEV_VENV:-}" ]; then export PATH="$DEV_VENV/bin:$PATH"; fi' \
        > /etc/profile.d/10-quicksynth-venv.sh \
    && chmod 0644 /etc/profile.d/10-quicksynth-venv.sh

USER ${USERNAME}
