#!/bin/sh
odir="/home/daminton/git/swiftest/examples/symba_test_restart/original"
rdir="/home/daminton/git/swiftest/examples/symba_test_restart/restart"
rsync -va $odir/dump_pl1.bin $odir/dump_tp1.bin $odir/bin.dat $odir/particle.dat $rdir/
grep TSTOP $odir/dump_param1.dat | sed 's/1\.0/10\.0/' > $rdir/param.in
grep -v TSTOP $odir/dump_param1.dat >> $rdir/param.in
