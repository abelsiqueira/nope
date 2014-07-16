#!/bin/bash

g=$(git ls-files .. | xargs grep -n "\s\+$")
if [ ! -z "$g" ]; then
  echo -e "\e[0;31mThere are files with trailing whitespaces\e[0m"
  echo -e "\e[0;33m$g\e[0m"
  exit 1
fi
