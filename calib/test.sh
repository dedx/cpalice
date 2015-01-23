#!/bin/sh                                                               
RUNLIST="177501 177497 177477 177496"
touch myfooB.txt
TAILLIST="002 003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022 023 024 025 026 027 028 029 030"
for RUN in $RUNLIST; do
for TAIL in $TAILLIST; do
  echo $RUN
      NUMS=10
      while [ $NUMS -lt 50 ]; do
          echo "/scratch2/alice/data/2012/LHC12b/000$RUN/ESDs/pass1/12000$RUN$TAIL.$NUMS/AliESDs.root" >> myfooB.txt
          NUMS=$(($NUMS+1))
      done
done
done

