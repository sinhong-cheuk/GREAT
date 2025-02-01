#!/bin/bash



alphas=(0.1 0.25)
for alpha in "${alphas[@]}"
do
    java -cp .:/home/zhuoxh/waiting_room-master/fastutil-7.2.0.jar Main "$alpha" 100000 &
    wait
done