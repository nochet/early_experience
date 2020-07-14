#!/usr/bin/sh

for f in ./*; do mv "$f" "${f%-*.JPG}.JPG" ; done
