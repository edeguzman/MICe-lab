export PATH=/OGS/bin/linux-x64:/axiom2/projects/software/arch/linux-xenial-xerus/bin/:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/micehome/matthijs/64-bit/bin:/micehome/matthijs/scripts:/axiom2/projects/software/arch/linux-xenial-xerus/src/mrtrix3/release/bin
export PYTHONPATH=/axiom2/projects/software/arch/linux-xenial-xerus/python
export LD_LIBRARY_PATH=/axiom2/projects/software/arch/linux-xenial-xerus/lib/:/axiom2/projects/software/arch/linux-xenial-xerus/lib/InsightToolkit/:/axiom2/projects/software/arch/linux-xenial-xerus/src/mrtrix3/release/lib
export PERL5LIB=/axiom2/projects/software/arch/linux-xenial-xerus/perl/
 
# first find all log files that need to be processed
logfiles=`find /hpf/largeprojects/MICe/HPF_CT_2_MINC/  -path  /hpf/largeprojects/MICe/HPF_CT_2_MINC/MINC  -prune -o -name '*.log' -print`;
 
datenow=`date`
echo
echo
echo
echo "*****************************************************"
echo "*********************************************************************************************"
echo "*****************************************************"
echo
echo
echo
echo $datenow
 
for log in $logfiles;
  do echo
  dir=`dirname $log`;
  tiffdir=`dirname $log`;
  dirbase=`basename $dir`;
  dir_minus_dirbase=`dirname $dir`
  sleep 2
  while [ $dir_minus_dirbase != /hpf/largeprojects/MICe/HPF_CT_2_MINC ];
    do cur_dirbase=$dirbase;
       next_dirbase=`basename $dir_minus_dirbase`
       #concatenate directory stucture information
       dirbase=${next_dirbase}_${dirbase}
       #update the directory name minus the latest dirbase
       dir_minus_dirbase=`dirname $dir_minus_dirbase`
       sleep 2
    done
  logbase=`basename $log .log`;
  lockfile=${log}__work_in_progess__;
  main_MINC_output=/hpf/largeprojects/MICe/HPF_CT_2_MINC/MINC/${dirbase}/${logbase}.mnc;
  main_MINC_output_dir=/hpf/largeprojects/MICe/HPF_CT_2_MINC/MINC/${dirbase}/;
  move_to_tiff_dir=/hpf/largeprojects/MICe/HPF_CT_2_MINC/MINC/processed_tiff_directories/${dirbase}/;
  user=`ls -ld $log | awk '{print $3}'`; 
  # create a lock-file in order to run a single conversion at the time:
  if [ ! -f ${lockfile} ]; then
    su $user -c 'touch '$lockfile'';
    su $user -c 'echo Going to convert files for log file: '$log' for user: '$user' creating MINC file: '$main_MINC_output'';
    su $user -c 'echo COMMAND: mkdir -p '$main_MINC_output_dir'';
    su $user -c 'mkdir -p '$main_MINC_output_dir'';
    su $user -c 'echo COMMAND: convert_CT_image.py --check-file-transfer-complete '$log'  '$main_MINC_output'';
    su $user -c '. /hpf/largeprojects/MICe/HPF_CT_2_MINC/MINC/setup_paths_for_conversion ; convert_CT_image.py --check-file-transfer-complete '$log'  '$main_MINC_output'';
    if [ -f $main_MINC_output ] ; then
      su $user -c 'echo        MINC FILE EXISTS  ***************    COMMAND: mv '${tiffdir}' '${move_to_tiff_dir}'';
      su $user -c 'mv '${tiffdir}' '${move_to_tiff_dir}''; fi;
  fi;
done

