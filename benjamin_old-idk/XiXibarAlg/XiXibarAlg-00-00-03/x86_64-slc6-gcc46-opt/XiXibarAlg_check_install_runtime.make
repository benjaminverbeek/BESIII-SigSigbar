#-- start of make_header -----------------

#====================================
#  Document XiXibarAlg_check_install_runtime
#
#   Generated Thu Nov 17 09:41:49 2022  by benjamin
#
#====================================

include ${CMTROOT}/src/Makefile.core

ifdef tag
CMTEXTRATAGS = $(tag)
else
tag       = $(CMTCONFIG)
endif

cmt_XiXibarAlg_check_install_runtime_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_XiXibarAlg_check_install_runtime_has_target_tag

tags      = $(tag),$(CMTEXTRATAGS),target_XiXibarAlg_check_install_runtime

XiXibarAlg_tag = $(tag)

#cmt_local_tagfile_XiXibarAlg_check_install_runtime = $(XiXibarAlg_tag)_XiXibarAlg_check_install_runtime.make
cmt_local_tagfile_XiXibarAlg_check_install_runtime = $(bin)$(XiXibarAlg_tag)_XiXibarAlg_check_install_runtime.make

else

tags      = $(tag),$(CMTEXTRATAGS)

XiXibarAlg_tag = $(tag)

#cmt_local_tagfile_XiXibarAlg_check_install_runtime = $(XiXibarAlg_tag).make
cmt_local_tagfile_XiXibarAlg_check_install_runtime = $(bin)$(XiXibarAlg_tag).make

endif

include $(cmt_local_tagfile_XiXibarAlg_check_install_runtime)
#-include $(cmt_local_tagfile_XiXibarAlg_check_install_runtime)

ifdef cmt_XiXibarAlg_check_install_runtime_has_target_tag

cmt_final_setup_XiXibarAlg_check_install_runtime = $(bin)setup_XiXibarAlg_check_install_runtime.make
cmt_dependencies_in_XiXibarAlg_check_install_runtime = $(bin)dependencies_XiXibarAlg_check_install_runtime.in
#cmt_final_setup_XiXibarAlg_check_install_runtime = $(bin)XiXibarAlg_XiXibarAlg_check_install_runtimesetup.make
cmt_local_XiXibarAlg_check_install_runtime_makefile = $(bin)XiXibarAlg_check_install_runtime.make

else

cmt_final_setup_XiXibarAlg_check_install_runtime = $(bin)setup.make
cmt_dependencies_in_XiXibarAlg_check_install_runtime = $(bin)dependencies.in
#cmt_final_setup_XiXibarAlg_check_install_runtime = $(bin)XiXibarAlgsetup.make
cmt_local_XiXibarAlg_check_install_runtime_makefile = $(bin)XiXibarAlg_check_install_runtime.make

endif

#cmt_final_setup = $(bin)setup.make
#cmt_final_setup = $(bin)XiXibarAlgsetup.make

#XiXibarAlg_check_install_runtime :: ;

dirs ::
	@if test ! -r requirements ; then echo "No requirements file" ; fi; \
	  if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi

javadirs ::
	@if test ! -d $(javabin) ; then $(mkdir) -p $(javabin) ; fi

srcdirs ::
	@if test ! -d $(src) ; then $(mkdir) -p $(src) ; fi

help ::
	$(echo) 'XiXibarAlg_check_install_runtime'

binobj = 
ifdef STRUCTURED_OUTPUT
binobj = XiXibarAlg_check_install_runtime/
#XiXibarAlg_check_install_runtime::
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
XiXibarAlg_check_install_runtime_once = 1
endif

ifdef XiXibarAlg_check_install_runtime_once

XiXibarAlg_check_install_runtimeactionstamp = $(bin)XiXibarAlg_check_install_runtime.actionstamp
#XiXibarAlg_check_install_runtimeactionstamp = XiXibarAlg_check_install_runtime.actionstamp

XiXibarAlg_check_install_runtime :: $(XiXibarAlg_check_install_runtimeactionstamp)
	$(echo) "XiXibarAlg_check_install_runtime ok"
#	@echo XiXibarAlg_check_install_runtime ok

#$(XiXibarAlg_check_install_runtimeactionstamp) :: $(XiXibarAlg_check_install_runtime_dependencies)
$(XiXibarAlg_check_install_runtimeactionstamp) ::
	$(silent) /cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.8/BesPolicy/BesPolicy-01-05-05/cmt/bes_check_installations.sh -files= -s=../share *.txt   -installdir=/share
	$(silent) cat /dev/null > $(XiXibarAlg_check_install_runtimeactionstamp)
#	@echo ok > $(XiXibarAlg_check_install_runtimeactionstamp)

XiXibarAlg_check_install_runtimeclean ::
	$(cleanup_silent) /bin/rm -f $(XiXibarAlg_check_install_runtimeactionstamp)

else

#XiXibarAlg_check_install_runtime :: $(XiXibarAlg_check_install_runtime_dependencies)
XiXibarAlg_check_install_runtime ::
	$(silent) /cvmfs/bes3.ihep.ac.cn/bes3sw/Boss/7.0.8/BesPolicy/BesPolicy-01-05-05/cmt/bes_check_installations.sh -files= -s=../share *.txt   -installdir=/share

endif

install ::
uninstall ::

#-- end of cmt_action_runner_header -----------------
#-- start of cleanup_header --------------

clean :: XiXibarAlg_check_install_runtimeclean ;
#	@cd .

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(XiXibarAlg_check_install_runtime.make) $@: No rule for such target" >&2
else
.DEFAULT::
	$(error PEDANTIC: $@: No rule for such target)
endif

XiXibarAlg_check_install_runtimeclean ::
#-- end of cleanup_header ---------------
