
noinst_PROGRAMS = vec2add vec2sub vec2dotp vec2scalp vec2norm vec2vecp \
vec2pred vec3add vec3sub vec3dotp vec3scalp vec3vecp vec3norm vec3pred \
mat3add mat4add prtree2pt prtree2nftrav prtree2nftrave prtree3pt \
prtree3nftrav prtree3nftrave prtree2hilb prtree3hilb

AM_CFLAGS = -I$(top_srcdir)/build
LDADD = $(top_srcdir)/build/vsg/libvsgd.la


$(srcdir)/package.m4: $(top_srcdir)/configure.ac
	{                                      \
	echo '# Signature of the current package.'; \
	echo 'm4_define([AT_PACKAGE_NAME],      [@PACKAGE_NAME@])'; \
	echo 'm4_define([AT_PACKAGE_TARNAME],   [@PACKAGE_TARNAME@])'; \
	echo 'm4_define([AT_PACKAGE_VERSION],   [@PACKAGE_VERSION@])'; \
	echo 'm4_define([AT_PACKAGE_STRING],    [@PACKAGE_STRING@])'; \
	echo 'm4_define([AT_PACKAGE_BUGREPORT], [@PACKAGE_BUGREPORT@])'; \
	} >$(srcdir)/package.m4

testscripts = vec2ops.at vec3ops.at mat3ops.at mat4ops.at prtree2.at \
prtree3.at

EXTRA_DIST = testsuite.at $(testscripts) testsuite hilbert2-lvl1.ref \
hilbert2-lvl2.ref hilbert2-lvl3.ref hilbert2-lvl4.ref hilbert2-lvl5.ref \
hilbert3-lvl1.ref hilbert3-lvl2.ref hilbert3-lvl3.ref hilbert3-lvl4.ref


TESTSUITE = $(srcdir)/testsuite
check-local: atconfig atlocal $(TESTSUITE)
	$(SHELL) $(TESTSUITE)

AUTOTEST = $(AUTOM4TE) --language=autotest
$(TESTSUITE): $(srcdir)/testsuite.at $(srcdir)/package.m4 $(testscripts) 
	$(AUTOTEST) -I $(srcdir) $@.at -o $@.tmp
	mv $@.tmp $@
