#!/usr/bin/sh

a=1
for i in *.JPG; do
  new=$(printf "%03d.JPG" "$a") #03 pad to length of 3
  mv -i -- "$i" "$i_$new"
  let a=a+1
done
