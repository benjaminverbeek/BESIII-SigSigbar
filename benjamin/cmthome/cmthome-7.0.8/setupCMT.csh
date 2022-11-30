setenv DEVROOT /cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/SLC6
setenv EXTERNALLIBDIR $DEVROOT/ExternalLib
setenv CONTRIBDIR $DEVROOT/contrib

if ( $?LD_LIBRARY_PATH == 0 ) then
    setenv LD_LIBRARY_PATH ""
endif
# setup-gcc
setenv PATH $CONTRIBDIR/gcc/4.6.3/x86_64-slc6/bin:$PATH
setenv LD_LIBRARY_PATH $CONTRIBDIR/gcc/4.6.3/x86_64-slc6/lib:$LD_LIBRARY_PATH
setenv LD_LIBRARY_PATH $CONTRIBDIR/gcc/4.6.3/x86_64-slc6/lib64:$LD_LIBRARY_PATH
# setup-cmt
source $CONTRIBDIR/CMT/v1r25/mgr/setup.csh

