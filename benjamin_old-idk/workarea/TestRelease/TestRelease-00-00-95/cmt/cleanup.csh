# echo "cleanup TestRelease TestRelease-00-00-95 in /home/benjamin/workarea"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtTestReleasetempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtTestReleasetempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt cleanup -csh -pack=TestRelease -version=TestRelease-00-00-95 -path=/home/benjamin/workarea  $* >${cmtTestReleasetempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt cleanup -csh -pack=TestRelease -version=TestRelease-00-00-95 -path=/home/benjamin/workarea  $* >${cmtTestReleasetempfile}"
  set cmtcleanupstatus=2
  /bin/rm -f ${cmtTestReleasetempfile}
  unset cmtTestReleasetempfile
  exit $cmtcleanupstatus
endif
set cmtcleanupstatus=0
source ${cmtTestReleasetempfile}
if ( $status != 0 ) then
  set cmtcleanupstatus=2
endif
/bin/rm -f ${cmtTestReleasetempfile}
unset cmtTestReleasetempfile
exit $cmtcleanupstatus

