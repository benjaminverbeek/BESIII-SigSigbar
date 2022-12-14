#-- start of make_header -----------------

#====================================
#  Document RhopiAlg_check_install_runtime
#
#   Generated Fri Nov 18 11:25:57 2022  by benjamin
#
#====================================

include ${CMTROOT}/src/Makefile.core

ifdef tag
CMTEXTRATAGS = $(tag)
else
tag       = $(CMTCONFIG)
endif

cmt_RhopiAlg_check_install_runtime_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_RhopiAlg_check_install_runtime_has_target_tag

tags      = $(tag),$(CMTEXTRATAGS),target_RhopiAlg_check_install_runtime

RhopiAlg_tag = $(tag)

#cmt_local_tagfile_RhopiAlg_check_install_runtime = $(RhopiAlg_tag)_RhopiAlg_check_install_runtime.make
cmt_local_tagfile_RhopiAlg_check_install_runtime = $(bin)$(RhopiAlg_tag)_RhopiAlg_check_install_runtime.make

else

tags      = $(tag),$(CMTEXTRATAGS)

RhopiAlg_tag = $(tag)

#cmt_local_tagfile_RhopiAlg_check_install_runtime = $(RhopiAlg_tag).make
cmt_local_tagfile_RhopiAlg_check_install_runtime = $(bin)$(RhopiAlg_tag).make

endif

include $(cmt_local_tagfile_RhopiAlg_check_install_runtime)
#-include $(cmt_local_tagfile_RhopiAlg_check_install_runtime)

ifdef cmt_RhopiAlg_check_install_runtime_has_target_tag

cmt_final_setup_RhopiAlg_check_install_runtime = $(bin)setup_RhopiAlg_check_install_runtime.make
cmt_dependencies_in_RhopiAlg_check_install_runtime = $(bin)dependencies_RhopiAlg_check_install_runtime.in
#cmt_final_setup_RhopiAlg_check_install_runtime = $(bin)RhopiAlg_RhopiAlg_check_install_runtimesetup.make
cmt_local_RhopiAlg_check_install_runtime_makefile = $(bin)RhopiAlg_check_install_runtime.make

else

cmt_final_setup_RhopiAlg_check_install_runtime = $(bin)setup.make
cmt_dependencies_in_RhopiAlg_check_install_runtime = $(bin)dependencies.in
#cmt_final_setup_RhopiAlg_check_install_runtime = $(bin)RhopiAlgsetup.make
cmt_local_RhopiAlg_check_install_runtime_makefile = $(bin)RhopiAlg_check_install_runtime.make

endif

#cmt_final_setup = $(bin)setup.make
#cmt_final_setup = $(bin)RhopiAlgsetup.make

#RhopiAlg_check_install_runtime :: ;

dirs ::
	@if test ! -r requirements ; then echo "No requirements file" ; fi; \
	  if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi

javadirs ::
	@if test ! -d $(javabin) ; then $(mkdir) -p $(javabin) ; fi

srcdirs ::
	@if test ! -d $(src) ; then $(mkdir) -p $(src) ; fi

help ::
	$(echo) 'RhopiAlg_check_install_runtime'

binobj = 
ifdef STRUCTURED_OUTPUT
binobj = RhopiAlg_check_install_runtime/
#RhopiAlg_check_install_runtime::
#	@if test ! -d $(bin)$(binobj) ; then $(mkdir) -p $(bin)$(binobj) ; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)$(binobj)
endif

${CMTROOT}/src/Makefile.core : ;
ifdef use_requirements
$(use_requirements) : ;
endif

#-- end of make_header ------------------
#-- start of cmt_action_runner_header ---------------

ifdef ONCE
RhopiAlg_check_install_runtime_once = 1
endif

ifdef RhopiAlg_check_install_runtime_once

RhopiAlg_check_install_runtimeactionstamp = $(bin)RhopiAlg_check_install_runtime.actionstamp
#RhopiAlg_check_install_runtimeactionstamp = RhopiAlg_check_install_runtime.actionstamp

RhopiAlg_check_install_runtime :: $(RhopiAlg_check_install_runtimeactionstamp)
	$(echo) "RhopiAlg_check_install_runtime ok"
#	@echo RhopiAlg_check_install_runtime ok

#$(RhopiAlg_check_install_runtimeactionstamp) :: $(RhopiAlg_check_install_runtime_dependencies)
$(RhopiAlg_check_install_runtimeactionstamp) ::
	$(silent) /cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.8/BesPolicy/BesPolicy-01-05-05/cmt/bes_check_installations.sh -files= -s=../share *.txt   -installdir=/home/benjamin/workarea/InstallArea/share
	$(silent) cat /dev/null > $(RhopiAlg_check_install_runtimeactionstamp)
#	@echo ok > $(RhopiAlg_check_install_runtimeactionstamp)

RhopiAlg_check_install_runtimeclean ::
	$(cleanup_silent) /bin/rm -f $(RhopiAlg_check_install_runtimeactionstamp)

else

#RhopiAlg_check_install_runtime :: $(RhopiAlg_check_install_runtime_dependencies)
RhopiAlg_check_install_runtime ::
	$(silent) /cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.8/BesPolicy/BesPolicy-01-05-05/cmt/bes_check_installations.sh -files= -s=../share *.txt   -installdir=/home/benjamin/workarea/InstallArea/share

endif

install ::
uninstall ::

#-- end of cmt_action_runner_header -----------------
#-- start of cleanup_header --------------

clean :: RhopiAlg_check_install_runtimeclean ;
#	@cd .

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(RhopiAlg_check_install_runtime.make) $@: No rule for such target" >&2
else
.DEFAULT::
	$(error PEDANTIC: $@: No rule for such target)
endif

RhopiAlg_check_install_runtimeclean ::
#-- end of cleanup_header ---------------
