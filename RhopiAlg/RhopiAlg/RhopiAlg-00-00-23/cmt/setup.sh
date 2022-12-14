# echo "setup RhopiAlg RhopiAlg-00-00-23 in /home/benjamin/workarea"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtRhopiAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtRhopiAlgtempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt setup -sh -pack=RhopiAlg -version=RhopiAlg-00-00-23 -path=/home/benjamin/workarea  -no_cleanup $* >${cmtRhopiAlgtempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt setup -sh -pack=RhopiAlg -version=RhopiAlg-00-00-23 -path=/home/benjamin/workarea  -no_cleanup $* >${cmtRhopiAlgtempfile}"
  cmtsetupstatus=2
  /bin/rm -f ${cmtRhopiAlgtempfile}
  unset cmtRhopiAlgtempfile
  return $cmtsetupstatus
fi
cmtsetupstatus=0
. ${cmtRhopiAlgtempfile}
if test $? != 0 ; then
  cmtsetupstatus=2
fi
/bin/rm -f ${cmtRhopiAlgtempfile}
unset cmtRhopiAlgtempfile
return $cmtsetupstatus

