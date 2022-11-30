# echo "setup XiXibarAlg XiXibarAlg-00-00-03 in /home/benjamin"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtXiXibarAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtXiXibarAlgtempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt setup -sh -pack=XiXibarAlg -version=XiXibarAlg-00-00-03 -path=/home/benjamin  -no_cleanup $* >${cmtXiXibarAlgtempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt setup -sh -pack=XiXibarAlg -version=XiXibarAlg-00-00-03 -path=/home/benjamin  -no_cleanup $* >${cmtXiXibarAlgtempfile}"
  cmtsetupstatus=2
  /bin/rm -f ${cmtXiXibarAlgtempfile}
  unset cmtXiXibarAlgtempfile
  return $cmtsetupstatus
fi
cmtsetupstatus=0
. ${cmtXiXibarAlgtempfile}
if test $? != 0 ; then
  cmtsetupstatus=2
fi
/bin/rm -f ${cmtXiXibarAlgtempfile}
unset cmtXiXibarAlgtempfile
return $cmtsetupstatus

