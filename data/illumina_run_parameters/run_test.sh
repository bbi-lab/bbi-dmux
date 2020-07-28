#!/bin/bash

PRG="../../bin/read_run_info.py"

echo "MiSeq unpaired read"
RUN_XML="MiSeq/read_single/200529_M00511_0140_000000000-CR4KR"
$PRG $RUN_XML
echo

echo "MiSeq paired read"
RUN_XML="MiSeq/read_pair/190423_M00296_0105_000000000-CCNJC"
$PRG $RUN_XML
echo

echo "MiSeq paired read"
RUN_XML="MiSeq/read_pair/200720_M00511_0150_000000000-J32RF"
$PRG $RUN_XML
echo

echo "HiSeq3000"
RUN_XML="HiSeq3000/130531_SN743_0384_Ac1vt1acxx"
$PRG $RUN_XML
echo

echo "HiSeq3000"
RUN_XML="HiSeq3000/180905_D00584_0284_AHLFNYBCX2"
$PRG $RUN_XML
echo

echo "NextSeq"
RUN_XML="NextSeq/200106_NS500488_0971_AHW5CYBGXC"
$PRG $RUN_XML
echo

echo "NextSeq2000"
RUN_XML="NextSeq2000/200722_VH00123_5_AAAFWVYM5"
$PRG $RUN_XML
echo

