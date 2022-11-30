# echo "setup XiXibarAlg XiXibarAlg-00-00-03 in /home/benjamin"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtXiXibarAlgtempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtXiXibarAlgtempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt setup -csh -pack=XiXibarAlg -version=XiXibarAlg-00-00-03 -path=/home/benjamin  -no_cleanup $* >${cmtXiXibarAlgtempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt setup -csh -pack=XiXibarAlg -version=XiXibarAlg-00-00-03 -path=/home/benjamin  -no_cleanup $* >${cmtXiXibarAlgtempfile}"
  set cmtsetupstatus=2
  /bin/rm -f ${cmtXiXibarAlgtempfile}
  unset cmtXiXibarAlgtempfile
  exit $cmtsetupstatus
endif
set cmtsetupstatus=0
source ${cmtXiXibarAlgtempfile}
if ( $status != 0 ) then
  set cmtsetupstatus=2
endif
/bin/rm -f ${cmtXiXibarAlgtempfile}
unset cmtXiXibarAlgtempfile
exit $cmtsetupstatus

