#!/bin/bash

cp revisedtests/TestsProg2/$1 app.c
make app
mpirun -np $2 --hostfile hostfile ./app

