# echo "cleanup XiXibarAlg XiXibarAlg-00-00-03 in /home/benjamin"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtXiXibarAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtXiXibarAlgtempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt cleanup -csh -pack=XiXibarAlg -version=XiXibarAlg-00-00-03 -path=/home/benjamin  $* >${cmtXiXibarAlgtempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt cleanup -csh -pack=XiXibarAlg -version=XiXibarAlg-00-00-03 -path=/home/benjamin  $* >${cmtXiXibarAlgtempfile}"
  set cmtcleanupstatus=2
  /bin/rm -f ${cmtXiXibarAlgtempfile}
  unset cmtXiXibarAlgtempfile
  exit $cmtcleanupstatus
endif
set cmtcleanupstatus=0
source ${cmtXiXibarAlgtempfile}
if ( $status != 0 ) then
  set cmtcleanupstatus=2
endif
/bin/rm -f ${cmtXiXibarAlgtempfile}
unset cmtXiXibarAlgtempfile
exit $cmtcleanupstatus

