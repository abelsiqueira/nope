#!/bin/bash

g=$(git ls-files .. | xargs grep "\s\+$")
[ -z "$g" ] || (echo $g && exit 1)
