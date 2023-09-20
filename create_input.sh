#!/bin/bash

dia_inicial=20110101
dia_final=20120731
expt=('expt_00.2')
rodada=('004j_002')
legenda=('FREE RUN')
inc_tempo=1

for ((i=0;i<${#expt[@]};i++)) do
    current=${dia_inicial}
    while [ ${current} -le ${dia_final} ]; do
           echo ${current}
           
           rm input_parameters.txt_shell_002 > /dev/null 2>&1
           echo ${expt[${i}]} >> input_parameters.txt_shell_002
           echo ${rodada[${i}]} >> input_parameters.txt_shell_002
           echo ${legenda[${i}]} >> input_parameters.txt_shell_002
           echo `date -u -d "${current}" +%d/%m/%Y` >> input_parameters.txt_shell_002
           last_day=`date -u -d "${current} +1 month -1 days" +%d/%m/%Y`
           echo ${last_day} >> input_parameters.txt_shell_002
           echo ${inc_tempo} >> input_parameters.txt_shell_002
           echo >> input_parameters.txt_shell_002

           python calc_RMSD_ROMS_HYCOM.py input_parameters.txt_shell_002
           wait

           current=`date -u -d "${current} +1 month" +%Y%m%d`
    done
done
