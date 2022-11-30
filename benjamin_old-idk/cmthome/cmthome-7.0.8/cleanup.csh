# echo "cleanup cmt_standalone v0 in /home/benjamin/cmthome/cmthome-7.0.8"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtcmt_standalonetempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtcmt_standalonetempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt cleanup -csh -pack=cmt_standalone -version=v0 -path=/home/benjamin/cmthome/cmthome-7.0.8  $* >${cmtcmt_standalonetempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt cleanup -csh -pack=cmt_standalone -version=v0 -path=/home/benjamin/cmthome/cmthome-7.0.8  $* >${cmtcmt_standalonetempfile}"
  set cmtcleanupstatus=2
  /bin/rm -f ${cmtcmt_standalonetempfile}
  unset cmtcmt_standalonetempfile
  exit $cmtcleanupstatus
endif
set cmtcleanupstatus=0
source ${cmtcmt_standalonetempfile}
if ( $status != 0 ) then
  set cmtcleanupstatus=2
endif
/bin/rm -f ${cmtcmt_standalonetempfile}
unset cmtcmt_standalonetempfile
exit $cmtcleanupstatus

