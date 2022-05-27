#!/bin/sh
if [ "$1" = "start" ]; then
    python3 ion_watchdog.py >> watchdog.log &
else
    PID=$(ps -ef |grep ion_watchdog |awk 'NR==1{print $2}')
    if [ ! -z "$PID" ];then kill $PID
    fi
fi