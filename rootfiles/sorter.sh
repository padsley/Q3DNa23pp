#!/bin/bash

for COUNT in {1..120}

do
    echo $COUNT
    if [ $COUNT -lt 10 ]
    then
	FILE="run00"$COUNT".root"
    elif [ $COUNT -lt 100 ]
    then
	FILE="run0"$COUNT".root"
    else
	FILE="run"$COUNT".root"
    fi
    echo $FILE
    if [ $COUNT -lt 10 ]
    then
	FILE2="run00"$COUNT"_an.root"
    elif [ $COUNT -lt 100 ]
    then
	FILE2="run0"$COUNT"_an.root"
    else
	FILE2="run"$COUNT"_an.root"
    fi
    echo $FILE2
    if [ -e $FILE ]
    then
	if [ ! -e $FILE2 ]
	then
	    echo $FILE
	    /home/mwilliams/Documents/MunichQ3DNa23pp/MunichQ3DSort/adsleyse $FILE
	fi
    fi
done
