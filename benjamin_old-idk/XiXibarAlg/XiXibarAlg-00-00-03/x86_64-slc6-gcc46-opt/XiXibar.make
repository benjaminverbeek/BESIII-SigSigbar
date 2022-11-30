#-- start of make_header -----------------

#====================================
#  Library XiXibar
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

cmt_XiXibar_has_no_target_tag = 1

#--------------------------------------------------------

ifdef cmt_XiXibar_has_target_tag

tags      = $(tag),$(CMTEXTRATAGS),target_XiXibar

XiXibarAlg_tag = $(tag)

#cmt_local_tagfile_XiXibar = $(XiXibarAlg_tag)_XiXibar.make
cmt_local_tagfile_XiXibar = $(bin)$(XiXibarAlg_tag)_XiXibar.make

else

tags      = $(tag),$(CMTEXTRATAGS)

XiXibarAlg_tag = $(tag)

#cmt_local_tagfile_XiXibar = $(XiXibarAlg_tag).make
cmt_local_tagfile_XiXibar = $(bin)$(XiXibarAlg_tag).make

endif

include $(cmt_local_tagfile_XiXibar)
#-include $(cmt_local_tagfile_XiXibar)

ifdef cmt_XiXibar_has_target_tag

cmt_final_setup_XiXibar = $(bin)setup_XiXibar.make
cmt_dependencies_in_XiXibar = $(bin)dependencies_XiXibar.in
#cmt_final_setup_XiXibar = $(bin)XiXibarAlg_XiXibarsetup.make
cmt_local_XiXibar_makefile = $(bin)XiXibar.make

else

cmt_final_setup_XiXibar = $(bin)setup.make
cmt_dependencies_in_XiXibar = $(bin)dependencies.in
#cmt_final_setup_XiXibar = $(bin)XiXibarAlgsetup.make
cmt_local_XiXibar_makefile = $(bin)XiXibar.make

endif

#cmt_final_setup = $(bin)setup.make
#cmt_final_setup = $(bin)XiXibarAlgsetup.make

#XiXibar :: ;

dirs ::
	@if test ! -r requirements ; then echo "No requirements file" ; fi; \
	  if test ! -d $(bin) ; then $(mkdir) -p $(bin) ; fi

javadirs ::
	@if test ! -d $(javabin) ; then $(mkdir) -p $(javabin) ; fi

srcdirs ::
	@if test ! -d $(src) ; then $(mkdir) -p $(src) ; fi

help ::
	$(echo) 'XiXibar'

binobj = 
ifdef STRUCTURED_OUTPUT
binobj = XiXibar/
#XiXibar::
#	@if test ! -d $(bin)$(binobj) ; then $(mkdir) -p $(bin)$(binobj) ; fi
#	$(echo) "STRUCTURED_OUTPUT="$(bin)$(binobj)
endif

${CMTROOT}/src/Makefile.core : ;
ifdef use_requirements
$(use_requirements) : ;
endif

#-- end of make_header ------------------
#-- start of libary_header ---------------

XiXibarlibname   = $(bin)$(library_prefix)XiXibar$(library_suffix)
XiXibarlib       = $(XiXibarlibname).a
XiXibarstamp     = $(bin)XiXibar.stamp
XiXibarshstamp   = $(bin)XiXibar.shstamp

XiXibar :: dirs  XiXibarLIB
	$(echo) "XiXibar ok"

#-- end of libary_header ----------------

XiXibarLIB :: $(XiXibarlib) $(XiXibarshstamp)
	@/bin/echo "------> XiXibar : library ok"

$(XiXibarlib) :: $(bin)XiXibar.o $(bin)LL_entries.o $(bin)LL_load.o
	$(lib_echo) library
	$(lib_silent) cd $(bin); \
	  $(ar) $(XiXibarlib) $?
	$(lib_silent) $(ranlib) $(XiXibarlib)
	$(lib_silent) cat /dev/null >$(XiXibarstamp)

#------------------------------------------------------------------
#  Future improvement? to empty the object files after
#  storing in the library
#
##	  for f in $?; do \
##	    rm $${f}; touch $${f}; \
##	  done
#------------------------------------------------------------------

$(XiXibarlibname).$(shlibsuffix) :: $(XiXibarlib) $(XiXibarstamps)
	$(lib_silent) cd $(bin); QUIET=$(QUIET); $(make_shlib) "$(tags)" XiXibar $(XiXibar_shlibflags)

$(XiXibarshstamp) :: $(XiXibarlibname).$(shlibsuffix)
	@if test -f $(XiXibarlibname).$(shlibsuffix) ; then cat /dev/null >$(XiXibarshstamp) ; fi

XiXibarclean ::
	$(cleanup_echo) objects
	$(cleanup_silent) cd $(bin); /bin/rm -f $(bin)XiXibar.o $(bin)LL_entries.o $(bin)LL_load.o

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
XiXibarinstallname = $(library_prefix)XiXibar$(library_suffix).$(shlibsuffix)

XiXibar :: XiXibarinstall

install :: XiXibarinstall

XiXibarinstall :: $(install_dir)/$(XiXibarinstallname)
	@if test ! "${installarea}" = ""; then\
	  echo "installation done"; \
	fi

$(install_dir)/$(XiXibarinstallname) :: $(bin)$(XiXibarinstallname)
	@if test ! "${installarea}" = ""; then \
	  cd $(bin); \
	  if test ! "$(install_dir)" = ""; then \
	    if test ! -d "$(install_dir)"; then \
	      mkdir -p $(install_dir); \
	    fi ; \
	    if test -d "$(install_dir)"; then \
	      echo "Installing library $(XiXibarinstallname) into $(install_dir)"; \
	      if test -e $(install_dir)/$(XiXibarinstallname); then \
	        $(cmt_uninstall_area_command) $(install_dir)/$(XiXibarinstallname); \
	        $(cmt_uninstall_area_command) $(install_dir)/$(XiXibarinstallname).cmtref; \
	      fi; \
	      $(cmt_install_area_command) `pwd`/$(XiXibarinstallname) $(install_dir)/$(XiXibarinstallname); \
	      echo `pwd`/$(XiXibarinstallname) >$(install_dir)/$(XiXibarinstallname).cmtref; \
	    fi \
          else \
	    echo "Cannot install library $(XiXibarinstallname), no installation directory specified"; \
	  fi; \
	fi

XiXibarclean :: XiXibaruninstall

uninstall :: XiXibaruninstall

XiXibaruninstall ::
	@if test ! "${installarea}" = ""; then \
	  cd $(bin); \
	  if test ! "$(install_dir)" = ""; then \
	    if test -d "$(install_dir)"; then \
	      echo "Removing installed library $(XiXibarinstallname) from $(install_dir)"; \
	      $(cmt_uninstall_area_command) $(install_dir)/$(XiXibarinstallname); \
	      $(cmt_uninstall_area_command) $(install_dir)/$(XiXibarinstallname).cmtref; \
	    fi \
          else \
	    echo "Cannot uninstall library $(XiXibarinstallname), no installation directory specified"; \
	  fi; \
	fi




#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),XiXibarclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)XiXibar.d

$(bin)$(binobj)XiXibar.d :

$(bin)$(binobj)XiXibar.o : $(cmt_final_setup_XiXibar)

$(bin)$(binobj)XiXibar.o : $(src)XiXibar.cxx
	$(cpp_echo) $(src)XiXibar.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(XiXibar_pp_cppflags) $(lib_XiXibar_pp_cppflags) $(XiXibar_pp_cppflags) $(use_cppflags) $(XiXibar_cppflags) $(lib_XiXibar_cppflags) $(XiXibar_cppflags) $(XiXibar_cxx_cppflags)  $(src)XiXibar.cxx
endif
endif

else
$(bin)XiXibar_dependencies.make : $(XiXibar_cxx_dependencies)

$(bin)XiXibar_dependencies.make : $(src)XiXibar.cxx

$(bin)$(binobj)XiXibar.o : $(XiXibar_cxx_dependencies)
	$(cpp_echo) $(src)XiXibar.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(XiXibar_pp_cppflags) $(lib_XiXibar_pp_cppflags) $(XiXibar_pp_cppflags) $(use_cppflags) $(XiXibar_cppflags) $(lib_XiXibar_cppflags) $(XiXibar_cppflags) $(XiXibar_cxx_cppflags)  $(src)XiXibar.cxx

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),XiXibarclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)LL_entries.d

$(bin)$(binobj)LL_entries.d :

$(bin)$(binobj)LL_entries.o : $(cmt_final_setup_XiXibar)

$(bin)$(binobj)LL_entries.o : $(src)components/LL_entries.cxx
	$(cpp_echo) $(src)components/LL_entries.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(XiXibar_pp_cppflags) $(lib_XiXibar_pp_cppflags) $(LL_entries_pp_cppflags) $(use_cppflags) $(XiXibar_cppflags) $(lib_XiXibar_cppflags) $(LL_entries_cppflags) $(LL_entries_cxx_cppflags) -I../src/components $(src)components/LL_entries.cxx
endif
endif

else
$(bin)XiXibar_dependencies.make : $(LL_entries_cxx_dependencies)

$(bin)XiXibar_dependencies.make : $(src)components/LL_entries.cxx

$(bin)$(binobj)LL_entries.o : $(LL_entries_cxx_dependencies)
	$(cpp_echo) $(src)components/LL_entries.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(XiXibar_pp_cppflags) $(lib_XiXibar_pp_cppflags) $(LL_entries_pp_cppflags) $(use_cppflags) $(XiXibar_cppflags) $(lib_XiXibar_cppflags) $(LL_entries_cppflags) $(LL_entries_cxx_cppflags) -I../src/components $(src)components/LL_entries.cxx

endif

#-- end of cpp_library ------------------
#-- start of cpp_library -----------------

ifneq (-MMD -MP -MF $*.d -MQ $@,)

ifneq ($(MAKECMDGOALS),XiXibarclean)
ifneq ($(MAKECMDGOALS),uninstall)
-include $(bin)$(binobj)LL_load.d

$(bin)$(binobj)LL_load.d :

$(bin)$(binobj)LL_load.o : $(cmt_final_setup_XiXibar)

$(bin)$(binobj)LL_load.o : $(src)components/LL_load.cxx
	$(cpp_echo) $(src)components/LL_load.cxx
	$(cpp_silent) $(cppcomp) -MMD -MP -MF $*.d -MQ $@ -o $@ $(use_pp_cppflags) $(XiXibar_pp_cppflags) $(lib_XiXibar_pp_cppflags) $(LL_load_pp_cppflags) $(use_cppflags) $(XiXibar_cppflags) $(lib_XiXibar_cppflags) $(LL_load_cppflags) $(LL_load_cxx_cppflags) -I../src/components $(src)components/LL_load.cxx
endif
endif

else
$(bin)XiXibar_dependencies.make : $(LL_load_cxx_dependencies)

$(bin)XiXibar_dependencies.make : $(src)components/LL_load.cxx

$(bin)$(binobj)LL_load.o : $(LL_load_cxx_dependencies)
	$(cpp_echo) $(src)components/LL_load.cxx
	$(cpp_silent) $(cppcomp) -o $@ $(use_pp_cppflags) $(XiXibar_pp_cppflags) $(lib_XiXibar_pp_cppflags) $(LL_load_pp_cppflags) $(use_cppflags) $(XiXibar_cppflags) $(lib_XiXibar_cppflags) $(LL_load_cppflags) $(LL_load_cxx_cppflags) -I../src/components $(src)components/LL_load.cxx

endif

#-- end of cpp_library ------------------
#-- start of cleanup_header --------------

clean :: XiXibarclean ;
#	@cd .

ifndef PEDANTIC
.DEFAULT::
	$(echo) "(XiXibar.make) $@: No rule for such target" >&2
else
.DEFAULT::
	$(error PEDANTIC: $@: No rule for such target)
endif

XiXibarclean ::
#-- end of cleanup_header ---------------
#-- start of cleanup_library -------------
	$(cleanup_echo) library XiXibar
	-$(cleanup_silent) cd $(bin); /bin/rm -f $(library_prefix)XiXibar$(library_suffix).a $(library_prefix)XiXibar$(library_suffix).s? XiXibar.stamp XiXibar.shstamp
#-- end of cleanup_library ---------------
