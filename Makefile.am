AUTOMAKE_OPTIONS = check-news
ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = volregrid

EXTRA_DIST = $(m4_files) \
	epm-header.in

m4_files = m4/mni_REQUIRE_LIB.m4 m4/mni_REQUIRE_MNILIBS.m4	 \
	m4/mni_REQUIRE_OPENINVENTOR.m4				 \
	m4/mni_cxx_have_koenig_lookup.m4 m4/smr_CGAL_MAKEFILE.m4 \
	m4/smr_OPTIONAL_LIB.m4 m4/smr_REQUIRED_LIB.m4		 \
	m4/smr_WITH_BUILD_PATH.m4

volregrid_SOURCES = volregrid.c \
   arb_path_io.h arb_path_io.c \
   minc_support.h minc_support.c \
   lex.l
