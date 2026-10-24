# ChannelBridge development image.
#
# Three kinds of dependency, handled three different ways, because they have
# three different licences and sizes:
#
#   * The Python package's own dependencies (numpy, scipy, matplotlib,
#     pillow, PyYAML, pytest, ruff) -- pip, from pyproject.toml, at
#     postCreateCommand so the lockstep with the repo is visible.
#   * PyMOL -- the distribution's open-source build. The pose and figure
#     stages shell out to it; they do not import it, so any PyMOL 2.x/3.x
#     works. Not pip-installable in a form worth depending on.
#   * CAVER -- NOT in this image. It is GPL-3 software from Masaryk
#     University and we do not redistribute it; see .devcontainer/README.md.
#     Only a JRE is installed here, so a downloaded CAVER runs as-is.
FROM mcr.microsoft.com/devcontainers/python:3.11-bookworm

# pymol here is the open-source build; it pulls a large Qt/GL stack, which is
# why this is a Dockerfile rather than a features-only devcontainer.json.
# default-jre-headless is for CAVER: it is pure Java and needs no display.
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && apt-get install -y --no-install-recommends \
        pymol \
        default-jre-headless \
        git \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Where a downloaded CAVER is expected. Created empty so the bind/copy target
# exists even before anyone puts CAVER there, and owned by the non-root user.
RUN mkdir -p /workspaces/ChannelBridge/vendor \
    && chown -R vscode:vscode /workspaces/ChannelBridge/vendor

USER vscode
