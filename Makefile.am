ACLOCAL_AMFLAGS = -I m4

SUBDIRS = include doc po

AM_CPPFLAGS = -Iinclude
AM_CPPFLAGS += -D__STRICT_ANSI__
AM_CPPFLAGS += -DNDEBUG

WARNFLAGS = -Wall
WARNFLAGS += -Wshadow
WARNFLAGS += -pedantic-errors -Wpointer-arith -Waggregate-return
WARNFLAGS += -Wlong-long -Wredundant-decls -Wcast-qual -Wcast-align
WARNFLAGS += -Wmissing-prototypes -Wstrict-prototypes -Wmissing-declarations

AM_CFLAGS = $(WARNFLAGS) $(OPENMP_CXXFLAGS) -std=c11
AM_CXXFLAGS = $(WARNFLAGS) -std=c++11

localedir = $(datadir)/locale
DEFS = -DLOCALEDIR=\"$(localedir)\" @DEFS@

lib_LTLIBRARIES = libminidjvu-mod.la
noinst_LTLIBRARIES = libminidjvu-mod-settings.la

libminidjvu_mod_la_SOURCES = src/matcher/no_mdjvu.h src/matcher/bitmaps.h	\
 src/matcher/common.h src/djvu/bs.h src/jb2/jb2coder.h			\
 src/jb2/bmpcoder.h src/jb2/zp.h src/jb2/jb2const.h			\
 src/base/mdjvucfg.h src/matcher/cuts.c src/matcher/patterns.c		\
 src/matcher/frames.c src/matcher/bitmaps.c src/alg/nosubst.c		\
 src/alg/erosion.c src/alg/smooth.c src/alg/delegate.c			\
 src/alg/classify.c src/alg/render.c src/alg/clean.c			\
 src/alg/adjust_y.c src/alg/blitsort.c src/alg/split.c			\
 src/alg/average.c src/alg/compress.c src/djvu/djvuload.c		\
 src/djvu/djvusave.c src/djvu/djvuinfo.c src/djvu/iff.c			\
 src/image-io/tiffload.c src/image-io/tiff.c src/image-io/pbm.c		\
 src/image-io/tiffsave.c src/image-io/bmp.c src/jb2/proto.c		\
 src/base/3graymap.c src/base/2io.c src/base/5image.c			\
 src/base/4bitmap.c src/base/version.c src/base/6string.c		\
 src/base/1error.c src/base/0porting.c src/djvu/djvudir.cpp		\
 src/djvu/bs.cpp src/jb2/jb2coder.cpp src/jb2/bmpcoder.cpp		\
 src/jb2/jb2load.cpp src/jb2/zp.cpp src/jb2/jb2save.cpp

libminidjvu_mod_settings_la_SOURCES = \
 tools/settings-reader/Arrays.cpp tools/settings-reader/Arrays.h 			\
 tools/settings-reader/atomic.cpp tools/settings-reader/atomic.h 			\
 tools/settings-reader/ByteStream.cpp tools/settings-reader/ByteStream.h 		\
 tools/settings-reader/DjVuGlobal.cpp tools/settings-reader/DjVuGlobal.h		\
 tools/settings-reader/GContainer.cpp tools/settings-reader/GContainer.h		\
 tools/settings-reader/GException.cpp tools/settings-reader/GException.h		\
 tools/settings-reader/GOS.cpp tools/settings-reader/GOS.h				\
 tools/settings-reader/GSmartPointer.cpp tools/settings-reader/GSmartPointer.h		\
 tools/settings-reader/GString.cpp tools/settings-reader/GString.h			\
 tools/settings-reader/GUnicode.cpp 							\
 tools/settings-reader/GURL.cpp tools/settings-reader/GURL.h				\
 tools/settings-reader/ParsingByteStream.cpp tools/settings-reader/ParsingBytestream.h 	\
 tools/settings-reader/UnicodeByteStream.cpp tools/settings-reader/UnicodeByteStream.h	\
 tools/settings-reader/SettingsReader.cpp tools/settings-reader/SettingsReader.h	\
 tools/settings-reader/DjbzSetting.cpp

 
bin_PROGRAMS = minidjvu-mod

minidjvu_mod_SOURCES = tools/minidjvu-mod.c

minidjvu_mod_LDADD = libminidjvu-mod.la libminidjvu-mod-settings.la


minidjvu-mod.pc:
	echo 'prefix=$(prefix)'			>  $@
	echo 'exec_prefix=$(exec_prefix)'	>> $@
	echo 'includedir=$(includedir)'		>> $@
	echo 'libdir=$(libdir)'			>> $@
	echo					>> $@
	echo 'Name: $(PACKAGE_TARNAME)'		>> $@
	echo 'Description: $(PACKAGE_NAME)'	>> $@
	echo 'URL: $(PACKAGE_URL)'		>> $@
	echo 'Version: $(PACKAGE_VERSION)'	>> $@
	echo 'Cflags: -I$(includedir)'		>> $@
	echo 'Libs: -L$(libdir) -lminidjvu-mod'	>> $@
	echo 'Libs.private: $(LDFLAGS) $(LIBS)'	>> $@

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = minidjvu-mod.pc

MOSTLYCLEANFILES = $(pkgconfig_DATA)
