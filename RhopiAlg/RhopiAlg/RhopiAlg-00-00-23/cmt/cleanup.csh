# echo "cleanup RhopiAlg RhopiAlg-00-00-23 in /home/benjamin/workarea"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtRhopiAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtRhopiAlgtempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt cleanup -csh -pack=RhopiAlg -version=RhopiAlg-00-00-23 -path=/home/benjamin/workarea  $* >${cmtRhopiAlgtempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt cleanup -csh -pack=RhopiAlg -version=RhopiAlg-00-00-23 -path=/home/benjamin/workarea  $* >${cmtRhopiAlgtempfile}"
  set cmtcleanupstatus=2
  /bin/rm -f ${cmtRhopiAlgtempfile}
  unset cmtRhopiAlgtempfile
  exit $cmtcleanupstatus
endif
set cmtcleanupstatus=0
source ${cmtRhopiAlgtempfile}
if ( $status != 0 ) then
  set cmtcleanupstatus=2
endif
/bin/rm -f ${cmtRhopiAlgtempfile}
unset cmtRhopiAlgtempfile
exit $cmtcleanupstatus

