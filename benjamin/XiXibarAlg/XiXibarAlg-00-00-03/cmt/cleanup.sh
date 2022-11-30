# echo "cleanup XiXibarAlg XiXibarAlg-00-00-03 in /home/benjamin"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtXiXibarAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtXiXibarAlgtempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt cleanup -sh -pack=XiXibarAlg -version=XiXibarAlg-00-00-03 -path=/home/benjamin  $* >${cmtXiXibarAlgtempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt cleanup -sh -pack=XiXibarAlg -version=XiXibarAlg-00-00-03 -path=/home/benjamin  $* >${cmtXiXibarAlgtempfile}"
  cmtcleanupstatus=2
  /bin/rm -f ${cmtXiXibarAlgtempfile}
  unset cmtXiXibarAlgtempfile
  return $cmtcleanupstatus
fi
cmtcleanupstatus=0
. ${cmtXiXibarAlgtempfile}
if test $? != 0 ; then
  cmtcleanupstatus=2
fi
/bin/rm -f ${cmtXiXibarAlgtempfile}
unset cmtXiXibarAlgtempfile
return $cmtcleanupstatus

