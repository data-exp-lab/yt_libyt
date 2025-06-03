#!/bin/bash
# This script downloads test data for yt_libyt
# Run this script inside the tests/data directory: sh get-test-data.sh

# gamer
mkdir -p gamer
curl -o "gamer/Plummer.tar.gz" "https://yt-project.org/data/Plummer.tar.gz"
tar -xvzf "gamer/Plummer.tar.gz" -C gamer

# enzo
mkdir -p enzo
curl -o "enzo/IsolatedGalaxy.tar.gz" "https://yt-project.org/data/IsolatedGalaxy.tar.gz"
tar -xvzf "enzo/IsolatedGalaxy.tar.gz" -C enzo
curl -o "enzo/EnzoKelvinHelmholtz.tar.gz" "https://yt-project.org/data/EnzoKelvinHelmholtz.tar.gz"
tar -xvzf "enzo/EnzoKelvinHelmholtz.tar.gz" -C enzo
