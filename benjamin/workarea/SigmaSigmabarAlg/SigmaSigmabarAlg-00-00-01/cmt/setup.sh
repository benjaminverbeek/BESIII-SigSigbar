# echo "setup SigmaSigmabarAlg SigmaSigmabarAlg-00-00-01 in /home/benjamin/workarea"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtSigmaSigmabarAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtSigmaSigmabarAlgtempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt setup -sh -pack=SigmaSigmabarAlg -version=SigmaSigmabarAlg-00-00-01 -path=/home/benjamin/workarea  -no_cleanup $* >${cmtSigmaSigmabarAlgtempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt setup -sh -pack=SigmaSigmabarAlg -version=SigmaSigmabarAlg-00-00-01 -path=/home/benjamin/workarea  -no_cleanup $* >${cmtSigmaSigmabarAlgtempfile}"
  cmtsetupstatus=2
  /bin/rm -f ${cmtSigmaSigmabarAlgtempfile}
  unset cmtSigmaSigmabarAlgtempfile
  return $cmtsetupstatus
fi
cmtsetupstatus=0
. ${cmtSigmaSigmabarAlgtempfile}
if test $? != 0 ; then
  cmtsetupstatus=2
fi
/bin/rm -f ${cmtSigmaSigmabarAlgtempfile}
unset cmtSigmaSigmabarAlgtempfile
return $cmtsetupstatus

