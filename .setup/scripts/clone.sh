#!/bin/bash

git clone git@github.com:jeksterslab/simStateSpace.git
rm -rf "$PWD.git"
mv simStateSpace/.git "$PWD"
rm -rf simStateSpace
