#!/bin/bash

# Install apt and R packages if they are not yet installed.

set -e

# a function to install apt packages only if they are not installed
function apt_install() {
    if ! dpkg -s "$@" >/dev/null 2>&1; then
        if [ "$(find /var/lib/apt/lists/* | wc -l)" = "0" ]; then
            apt-get update
        fi
        apt-get install -y --no-install-recommends "$@"
    fi
}

# To add additional apt packages,
# for example, to install tmux and vim, do the following:
apt_install \
  tmux      \
  vim

# radian
apt_install python3-pip
pip3 install -U radian

# To add additional R packages,
# for example, to install tidyverse and dplyr
# do the following:
install2.r --error --skipinstalled -n -1 \
  tidyverse \
  dplyr
