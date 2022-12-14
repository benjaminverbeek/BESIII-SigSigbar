# echo "setup RhopiAlg RhopiAlg-00-00-23 in /home/benjamin/workarea"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtRhopiAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtRhopiAlgtempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt setup -csh -pack=RhopiAlg -version=RhopiAlg-00-00-23 -path=/home/benjamin/workarea  -no_cleanup $* >${cmtRhopiAlgtempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt setup -csh -pack=RhopiAlg -version=RhopiAlg-00-00-23 -path=/home/benjamin/workarea  -no_cleanup $* >${cmtRhopiAlgtempfile}"
  set cmtsetupstatus=2
  /bin/rm -f ${cmtRhopiAlgtempfile}
  unset cmtRhopiAlgtempfile
  exit $cmtsetupstatus
endif
set cmtsetupstatus=0
source ${cmtRhopiAlgtempfile}
if ( $status != 0 ) then
  set cmtsetupstatus=2
endif
/bin/rm -f ${cmtRhopiAlgtempfile}
unset cmtRhopiAlgtempfile
exit $cmtsetupstatus

