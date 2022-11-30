# echo "setup SigmaSigmabarAlg SigmaSigmabarAlg-00-00-01 in /home/benjamin/workarea"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtSigmaSigmabarAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtSigmaSigmabarAlgtempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt setup -csh -pack=SigmaSigmabarAlg -version=SigmaSigmabarAlg-00-00-01 -path=/home/benjamin/workarea  -no_cleanup $* >${cmtSigmaSigmabarAlgtempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt setup -csh -pack=SigmaSigmabarAlg -version=SigmaSigmabarAlg-00-00-01 -path=/home/benjamin/workarea  -no_cleanup $* >${cmtSigmaSigmabarAlgtempfile}"
  set cmtsetupstatus=2
  /bin/rm -f ${cmtSigmaSigmabarAlgtempfile}
  unset cmtSigmaSigmabarAlgtempfile
  exit $cmtsetupstatus
endif
set cmtsetupstatus=0
source ${cmtSigmaSigmabarAlgtempfile}
if ( $status != 0 ) then
  set cmtsetupstatus=2
endif
/bin/rm -f ${cmtSigmaSigmabarAlgtempfile}
unset cmtSigmaSigmabarAlgtempfile
exit $cmtsetupstatus

