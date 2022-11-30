#-- start of make_header -----------------

#====================================
#  Document config
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

cmt_config_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_config_has_target_tag

tags      = $(tag),$(CMTEXTRATAGS),target_config

XiXibarAlg_tag = $(tag)

#cmt_local_tagfile_config = $(XiXibarAlg_tag)_config.make
cmt_local_tagfile_config = $(bin)$(XiXibarAlg_tag)_config.make

else

tags      = $(tag),$(CMTEXTRATAGS)

XiXibarAlg_tag = $(tag)

#cmt_local_tagfile_config = $(XiXibarAlg_tag).make
cmt_local_tagfile_config = $(bin)$(XiXibarAlg_tag).make

endif

include $(cmt_local_tagfile_config)
#-include $(cmt_local_tagfile_config)

ifdef cmt_config_has_target_tag

cmt_final_setup_config = $(bin)setup_config.make
cmt_dependencies_in_config = $(bin)dependencies_config.in
#cmt_final_setup_config = $(bin)XiXibarAlg_configsetup.make
cmt_local_config_makefile = $(bin)config.make

else

cmt_final_setup_config = $(bin)setup.make
cmt_dependencies_in_config = $(bin)dependencies.in
#cmt_final_setup_config = $(bin)XiXibarAlgsetup.make
cmt_local_config_makefile = $(bin)config.make

endif

#cmt_final_setup = $(bin)setup.make
#cmt_final_setup = $(bin)XiXibarAlgsetup.make

#config :: ;

dirs ::
	@if test ! -r requirements ; then echo "No requirements file" ; fi; \
	  if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi

javadirs ::
	@if test ! -d $(javabin) ; then $(mkdir) -p $(javabin) ; fi

srcdirs ::
	@if test ! -d $(src) ; then $(mkdir) -p $(src) ; fi

help ::
	$(echo) 'config'

binobj = 
ifdef STRUCTURED_OUTPUT
binobj = config/
#config::
#	@if test ! -d $(bin)$(binobj) ; then $(mkdir) -p $(bin)$(binobj) ; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)$(binobj)
endif

${CMTROOT}/src/Makefile.core : ;
ifdef use_requirements
$(use_requirements) : ;
endif

#-- end of make_header ------------------

config :: ../XiXibarAlg/config.h
	@/bin/echo "------> config.h ok"

../XiXibarAlg/config.h :: ../XiXibarAlg/config.h.in
	@if test -f ../XiXibarAlg/config.h.in; then \
	  cd ../XiXibarAlg; \
	  $(config_command) config.h.in config.h; \
        fi


#-- start of cleanup_header --------------

clean :: configclean ;
#	@cd .

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(config.make) $@: No rule for such target" >&2
else
.DEFAULT::
	$(error PEDANTIC: $@: No rule for such target)
endif

configclean ::
#-- end of cleanup_header ---------------
