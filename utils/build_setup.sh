#!/usr/bin/env bash

if [["$TRAVIS_OS_NAME" == "linux"]]; then
    sudo add-apt-repository --yes ppa:ubuntu-sdk-team/ppa
    sudo apt-get update -qq
    sudo apt-get install -qq x11-apps libxt6 libxrender1 libxext6 libgl1-mesa-glx qtbase5-dev qtdeclarative5-dev
    sudo apt-get install -y --no-install-recommends texlive-fonts-recommended texlive-latex-extra texlive-fonts-extra dvipng texlive-latex-recommended
    export GKSwstype=svg
elif  [[ "$TRAVIS_OS_NAME" == "osx" ]];
    brew update
    brew install qt5
    brew cask install basictex
fi