# echo "cleanup SigmaSigmabarAlg SigmaSigmabarAlg-00-00-01 in /home/benjamin/workarea"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtSigmaSigmabarAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtSigmaSigmabarAlgtempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt cleanup -sh -pack=SigmaSigmabarAlg -version=SigmaSigmabarAlg-00-00-01 -path=/home/benjamin/workarea  $* >${cmtSigmaSigmabarAlgtempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt cleanup -sh -pack=SigmaSigmabarAlg -version=SigmaSigmabarAlg-00-00-01 -path=/home/benjamin/workarea  $* >${cmtSigmaSigmabarAlgtempfile}"
  cmtcleanupstatus=2
  /bin/rm -f ${cmtSigmaSigmabarAlgtempfile}
  unset cmtSigmaSigmabarAlgtempfile
  return $cmtcleanupstatus
fi
cmtcleanupstatus=0
. ${cmtSigmaSigmabarAlgtempfile}
if test $? != 0 ; then
  cmtcleanupstatus=2
fi
/bin/rm -f ${cmtSigmaSigmabarAlgtempfile}
unset cmtSigmaSigmabarAlgtempfile
return $cmtcleanupstatus

