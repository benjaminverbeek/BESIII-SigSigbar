#-- start of make_header -----------------

#====================================
#  Library RhopiAlg
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

cmt_RhopiAlg_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_RhopiAlg_has_target_tag

tags      = $(tag),$(CMTEXTRATAGS),target_RhopiAlg

RhopiAlg_tag = $(tag)

#cmt_local_tagfile_RhopiAlg = $(RhopiAlg_tag)_RhopiAlg.make
cmt_local_tagfile_RhopiAlg = $(bin)$(RhopiAlg_tag)_RhopiAlg.make

else

tags      = $(tag),$(CMTEXTRATAGS)

RhopiAlg_tag = $(tag)

#cmt_local_tagfile_RhopiAlg = $(RhopiAlg_tag).make
cmt_local_tagfile_RhopiAlg = $(bin)$(RhopiAlg_tag).make

endif

include $(cmt_local_tagfile_RhopiAlg)
#-include $(cmt_local_tagfile_RhopiAlg)

ifdef cmt_RhopiAlg_has_target_tag

cmt_final_setup_RhopiAlg = $(bin)setup_RhopiAlg.make
cmt_dependencies_in_RhopiAlg = $(bin)dependencies_RhopiAlg.in
#cmt_final_setup_RhopiAlg = $(bin)RhopiAlg_RhopiAlgsetup.make
cmt_local_RhopiAlg_makefile = $(bin)RhopiAlg.make

else

cmt_final_setup_RhopiAlg = $(bin)setup.make
cmt_dependencies_in_RhopiAlg = $(bin)dependencies.in
#cmt_final_setup_RhopiAlg = $(bin)RhopiAlgsetup.make
cmt_local_RhopiAlg_makefile = $(bin)RhopiAlg.make

endif

#cmt_final_setup = $(bin)setup.make
#cmt_final_setup = $(bin)RhopiAlgsetup.make

#RhopiAlg :: ;

dirs ::
	@if test ! -r requirements ; then echo "No requirements file" ; fi; \
	  if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi

javadirs ::
	@if test ! -d $(javabin) ; then $(mkdir) -p $(javabin) ; fi

srcdirs ::
	@if test ! -d $(src) ; then $(mkdir) -p $(src) ; fi

help ::
	$(echo) 'RhopiAlg'

binobj = 
ifdef STRUCTURED_OUTPUT
binobj = RhopiAlg/
#RhopiAlg::
#	@if test ! -d $(bin)$(binobj) ; then $(mkdir) -p $(bin)$(binobj) ; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)$(binobj)
endif

${CMTROOT}/src/Makefile.core : ;
ifdef use_requirements
$(use_requirements) : ;
endif

#-- end of make_header ------------------
#-- start of libary_header ---------------

RhopiAlglibname   = $(bin)$(library_prefix)RhopiAlg$(library_suffix)
RhopiAlglib       = $(RhopiAlglibname).a
RhopiAlgstamp     = $(bin)RhopiAlg.stamp
RhopiAlgshstamp   = $(bin)RhopiAlg.shstamp

RhopiAlg :: dirs  RhopiAlgLIB
	$(echo) "RhopiAlg ok"

#-- end of libary_header ----------------

RhopiAlgLIB :: $(RhopiAlglib) $(RhopiAlgshstamp)
	@/bin/echo "------> RhopiAlg : library ok"

$(RhopiAlglib) :: $(bin)Rhopi.o $(bin)Rhopi_load.o $(bin)Rhopi_entries.o
	$(lib_echo) library
	$(lib_silent) cd $(bin); \
	  $(ar) $(RhopiAlglib) $?
	$(lib_silent) $(ranlib) $(RhopiAlglib)
	$(lib_silent) cat /dev/null >$(RhopiAlgstamp)

#------------------------------------------------------------------
#  Future improvement? to empty the object files after
#  storing in the library
#
##	  for f in $?; do \
##	    rm $${f}; touch $${f}; \
##	  done
#------------------------------------------------------------------

$(RhopiAlglibname).$(shlibsuffix) :: $(RhopiAlglib) $(RhopiAlgstamps)
	$(lib_silent) cd $(bin); QUIET=$(QUIET); $(make_shlib) "$(tags)" RhopiAlg $(RhopiAlg_shlibflags)

$(RhopiAlgshstamp) :: $(RhopiAlglibname).$(shlibsuffix)
	@if test -f $(RhopiAlglibname).$(shlibsuffix) ; then cat /dev/null >$(RhopiAlgshstamp) ; fi

RhopiAlgclean ::
	$(cleanup_echo) objects
	$(cleanup_silent) cd $(bin); /bin/rm -f $(bin)Rhopi.o $(bin)Rhopi_load.o $(bin)Rhopi_entries.o

#-----------------------------------------------------------------
#
#  New section for automatic installation
#
#-----------------------------------------------------------------

ifeq ($(INSTALLAREA),)
installarea = $(CMTINSTALLAREA)
else
ifeq ($(findstring `,$(INSTALLAREA)),`)
installarea = $(shell $(subst `,, $(INSTALLAREA)))
else
installarea = $(INSTALLAREA)
endif
endif

install_dir = ${installarea}/${CMTCONFIG}/lib
RhopiAlginstallname = $(library_prefix)RhopiAlg$(library_suffix).$(shlibsuffix)

RhopiAlg :: RhopiAlginstall

install :: RhopiAlginstall

RhopiAlginstall :: $(install_dir)/$(RhopiAlginstallname)
	@if test ! "${installarea}" = ""; then\
	  echo "installation done"; \
	fi

$(install_dir)/$(RhopiAlginstallname) :: $(bin)$(RhopiAlginstallname)
	@if test ! "${installarea}" = ""; then \
	  cd $(bin); \
	  if test ! "$(install_dir)" = ""; then \
	    if test ! -d "$(install_dir)"; then \
	      mkdir -p $(install_dir); \
	    fi ; \
	    if test -d "$(install_dir)"; then \
	      echo "Installing library $(RhopiAlginstallname) into $(install_dir)"; \
	      if test -e $(install_dir)/$(RhopiAlginstallname); then \
	        $(cmt_uninstall_area_command) $(install_dir)/$(RhopiAlginstallname); \
	        $(cmt_uninstall_area_command) $(install_dir)/$(RhopiAlginstallname).cmtref; \
	      fi; \
	      $(cmt_install_area_command) `pwd`/$(RhopiAlginstallname) $(install_dir)/$(RhopiAlginstallname); \
	      echo `pwd`/$(RhopiAlginstallname) >$(install_dir)/$(RhopiAlginstallname).cmtref; \
	    fi \
          else \
	    echo "Cannot install library $(RhopiAlginstallname), no installation directory specified"; \
	  fi; \
	fi

RhopiAlgclean :: RhopiAlguninstall

uninstall :: RhopiAlguninstall

RhopiAlguninstall ::
	@if test ! "${installarea}" = ""; then \
	  cd $(bin); \
	  if test ! "$(install_dir)" = ""; then \
	    if test -d "$(install_dir)"; then \
	      echo "Removing installed library $(RhopiAlginstallname) from $(install_dir)"; \
	      $(cmt_uninstall_area_command) $(install_dir)/$(RhopiAlginstallname); \
	      $(cmt_uninstall_area_command) $(install_dir)/$(RhopiAlginstallname).cmtref; \
	    fi \
          else \
	    echo "Cannot uninstall library $(RhopiAlginstallname), no installation directory specified"; \
	  fi; \
	fi




#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),RhopiAlgclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)Rhopi.d

$(bin)$(binobj)Rhopi.d :

$(bin)$(binobj)Rhopi.o : $(cmt_final_setup_RhopiAlg)

$(bin)$(binobj)Rhopi.o : $(src)Rhopi.cxx
	$(cpp_echo) $(src)Rhopi.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(RhopiAlg_pp_cppflags) $(lib_RhopiAlg_pp_cppflags) $(Rhopi_pp_cppflags) $(use_cppflags) $(RhopiAlg_cppflags) $(lib_RhopiAlg_cppflags) $(Rhopi_cppflags) $(Rhopi_cxx_cppflags)  $(src)Rhopi.cxx
endif
endif

else
$(bin)RhopiAlg_dependencies.make : $(Rhopi_cxx_dependencies)

$(bin)RhopiAlg_dependencies.make : $(src)Rhopi.cxx

$(bin)$(binobj)Rhopi.o : $(Rhopi_cxx_dependencies)
	$(cpp_echo) $(src)Rhopi.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(RhopiAlg_pp_cppflags) $(lib_RhopiAlg_pp_cppflags) $(Rhopi_pp_cppflags) $(use_cppflags) $(RhopiAlg_cppflags) $(lib_RhopiAlg_cppflags) $(Rhopi_cppflags) $(Rhopi_cxx_cppflags)  $(src)Rhopi.cxx

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),RhopiAlgclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)Rhopi_load.d

$(bin)$(binobj)Rhopi_load.d :

$(bin)$(binobj)Rhopi_load.o : $(cmt_final_setup_RhopiAlg)

$(bin)$(binobj)Rhopi_load.o : $(src)components/Rhopi_load.cxx
	$(cpp_echo) $(src)components/Rhopi_load.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(RhopiAlg_pp_cppflags) $(lib_RhopiAlg_pp_cppflags) $(Rhopi_load_pp_cppflags) $(use_cppflags) $(RhopiAlg_cppflags) $(lib_RhopiAlg_cppflags) $(Rhopi_load_cppflags) $(Rhopi_load_cxx_cppflags) -I../src/components $(src)components/Rhopi_load.cxx
endif
endif

else
$(bin)RhopiAlg_dependencies.make : $(Rhopi_load_cxx_dependencies)

$(bin)RhopiAlg_dependencies.make : $(src)components/Rhopi_load.cxx

$(bin)$(binobj)Rhopi_load.o : $(Rhopi_load_cxx_dependencies)
	$(cpp_echo) $(src)components/Rhopi_load.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(RhopiAlg_pp_cppflags) $(lib_RhopiAlg_pp_cppflags) $(Rhopi_load_pp_cppflags) $(use_cppflags) $(RhopiAlg_cppflags) $(lib_RhopiAlg_cppflags) $(Rhopi_load_cppflags) $(Rhopi_load_cxx_cppflags) -I../src/components $(src)components/Rhopi_load.cxx

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),RhopiAlgclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)Rhopi_entries.d

$(bin)$(binobj)Rhopi_entries.d :

$(bin)$(binobj)Rhopi_entries.o : $(cmt_final_setup_RhopiAlg)

$(bin)$(binobj)Rhopi_entries.o : $(src)components/Rhopi_entries.cxx
	$(cpp_echo) $(src)components/Rhopi_entries.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(RhopiAlg_pp_cppflags) $(lib_RhopiAlg_pp_cppflags) $(Rhopi_entries_pp_cppflags) $(use_cppflags) $(RhopiAlg_cppflags) $(lib_RhopiAlg_cppflags) $(Rhopi_entries_cppflags) $(Rhopi_entries_cxx_cppflags) -I../src/components $(src)components/Rhopi_entries.cxx
endif
endif

else
$(bin)RhopiAlg_dependencies.make : $(Rhopi_entries_cxx_dependencies)

$(bin)RhopiAlg_dependencies.make : $(src)components/Rhopi_entries.cxx

$(bin)$(binobj)Rhopi_entries.o : $(Rhopi_entries_cxx_dependencies)
	$(cpp_echo) $(src)components/Rhopi_entries.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(RhopiAlg_pp_cppflags) $(lib_RhopiAlg_pp_cppflags) $(Rhopi_entries_pp_cppflags) $(use_cppflags) $(RhopiAlg_cppflags) $(lib_RhopiAlg_cppflags) $(Rhopi_entries_cppflags) $(Rhopi_entries_cxx_cppflags) -I../src/components $(src)components/Rhopi_entries.cxx

endif

#-- end of cpp_library ------------------
#-- start of cleanup_header --------------

clean :: RhopiAlgclean ;
#	@cd .

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(RhopiAlg.make) $@: No rule for such target" >&2
else
.DEFAULT::
	$(error PEDANTIC: $@: No rule for such target)
endif

RhopiAlgclean ::
#-- end of cleanup_header ---------------
#-- start of cleanup_library -------------
	$(cleanup_echo) library RhopiAlg
	-$(cleanup_silent) cd $(bin); /bin/rm -f $(library_prefix)RhopiAlg$(library_suffix).a $(library_prefix)RhopiAlg$(library_suffix).s? RhopiAlg.stamp RhopiAlg.shstamp
#-- end of cleanup_library ---------------
