#include "SettingsReader.h"
#include "ParsingByteStream.h"
#include "GSmartPointer.h"
#include "GURL.h"
#include "GString.h"
#include <iostream>
#include <memory>

inline bool is_valid_classifier(int id) { return id > 0 && id < 4; }

SettingsReader::SettingsReader(const char* settings_file, AppOptions *appOptions):
    m_bs(nullptr), m_appOptions(appOptions)
{
    const GURL url = GURL::Filename::UTF8(GNativeString(settings_file));
    GP<ByteStream> bs = ByteStream::create(url, "r");

    //NOTE: this require a change in ParsingByteStream to keep bs from destruction
    m_bs = new GP<ParsingByteStream>(ParsingByteStream::create(bs));
}

SettingsReader::~SettingsReader()
{
    if (m_bs) {
        delete m_bs;
    }
}

bool
SettingsReader::readAllOptions()
{
    if (!m_bs || !m_appOptions) {
        return false;
    }

    ParsingByteStream& pbs = **m_bs;

    int c = pbs.get_spaces(true);
    if (c == EOF) {
        std::cerr << "Empty file or can't be read\n";
        return false;
    }

    if (c!='(') {
        std::cerr << "Syntax error in settings file: expecting '('\n";
        return false;
    }

    int app_options_found = 0;
    int input_files_found = 0;

    while (c != EOF) {

        GUTF8String token = pbs.get_token();
        if (token == "options") {
            app_options_found = 1;
            if ( input_files_found ) {
                std::cerr << "Error: \"options\" list must be declared before \"input-files\" lists.\n";
                return false;
            }
            if (!readAppOptions())
                return false;
        } else if (token == "input-files") {
            input_files_found = 1;
            if (!readInputFiles(&m_appOptions->file_list, false))
                return false;
        } else if (token == "djbz") {
            if ( !app_options_found || !input_files_found ) {
                std::cerr << "Error: \"djbz\" sttings lists must be declared after \"options\" and \"input-files\" lists.\n";
                return false;
            }
            struct DjbzOptions* djbz = djbz_setting_create(m_appOptions->default_djbz_options);
            if (!readDjbzOptions(djbz, false)) {
                djbz_setting_destroy(djbz);
                return false;
            }
            djbz_list_add_option(&m_appOptions->djbz_list, djbz);
        } else {
            std::cerr << "Syntax error in settings file: expecting '(djbz'\n";
            return false;
        }
        c = pbs.get_spaces(true);
    }

    return true;
}

/* readAppOptions()
(options             # application options and defaults

 (default-djbz       # default djbz settings
   averaging     0   # default averaging (off)
   aggression    100 # default agression level (100)
   classifier    3   # default classifier (max compression)
   erosion       0   # default erosion (disabled)
   no-prototypes 0   # default prototypes usage (on)
   xtension      djbz # default djbz id extension ("djbz")
 )

 (default-image      # default image options

   #dpi           300 # if set, use this dpi value for encoding all images
                     # except those that have personal dpi option set.
                     # if not set, use dpi of source image of each page.

   smooth       0    # default smoothing image before processing (off)
   clean        0    # default cleaning image after processing (off)
   erosion      0    # default erosion image after processing (off)
 )


 indirect       0    # save indirect djvu (multifile) (off)
 #lossy          1    # if set, turns off or on following options:
                     # default-djbz::erosion, default-djbz::averaging
                     # default-image::smooth, default-image::clean

 match
 pages-per-dict 10   # automatically assign pages that aren't referred
                     # in any djbz blocks to the new djbz dictionaries.
                     # New dictionaries contain 10 (default) pages or less.

 report         0    # report progress to stdout
 #threads-max   2    # if set, use max N threads for processing (each thread
                     # process one block of pages. One djbz is a one block).
                     # By default, if CPU have C cores:
                     # if C > 2 then N = C-1, otherwise N = 1
 verbose        1    # print verbose log to stdout
 warnings       1    # print libtiff warnings to stdout
)
*/

bool
SettingsReader::readAppOptions()
{
    ParsingByteStream& pbs = **m_bs;

    GUTF8String token; GUTF8String val; int i;
    while (true) {
        token = pbs.get_token(true).downcase();
        i = token.search(')');
        if (i > 0) {
            for (int j = token.length()-1; j >= i; j--) {
                pbs.unget(token[j]);
            }
            token = token.substr(0, i).downcase();
        }

        if (!token || token[0] == ')') {
            break;
        }

        if (token[0] == '(') {
            //nested list
            const int len = token.length();
            token = (len > 1) ? token.substr(1, len-1)          // "(file"
                              : pbs.get_token(true).downcase(); // "( file"

            if (token == "default-djbz") {
                if (!readDjbzOptions(m_appOptions->default_djbz_options, true)) {
                    return false;
                }
                continue;
            } else if (token == "default-image") {
                if (!readImageOptions(m_appOptions->default_image_options)) {
                    return false;
                }
                continue;
            } else {
                std::cerr << "ERROR: only \"default-djbz\" list could be nested in \"options\" list, but " <<
                             token.getUTF82Native().getbuf() << " found!\n";
                return false;
            }
        }



        if (token == "indirect") {
            if (!readValInt("indirect", m_appOptions->indirect)) return false;
        } else if (token == "lossy") {
            int dummy = 0;
            if (!readValInt("lossy", dummy, 0, 2)) return false;
            if (dummy == 2) {
                m_appOptions->Match = 1;
                dummy = 1;
            }
            m_appOptions->match =
                    m_appOptions->default_djbz_options->erosion =
                    m_appOptions->default_djbz_options->averaging =
                    m_appOptions->default_image_options->smooth =
                    m_appOptions->default_image_options->clean = dummy;
        } else if (token == "match") {
            if (!readValInt("match", m_appOptions->match, 0, 2)) return false;
            if (m_appOptions->match == 2) {
                m_appOptions->match = m_appOptions->Match = 1;
            }
        } else if (token == "pages-per-dict") {
            if (!readValInt("pages_per_dict", m_appOptions->pages_per_dict, 0, 1e6)) return false;
        } else if (token == "report") {
            if (!readValInt("report", m_appOptions->report)) return false;
        } else
#ifdef _OPENMP
            if (token == "threads-max") {
                if (!readValInt("threads-max", m_appOptions->max_threads, 1, 1e6)) return false;
            } else
#endif
                if (token == "verbose") {
                    if (!readValInt("verbose", m_appOptions->verbose)) return false;
                } else if (token == "xtension") {
                    if (!readValStr("xtension", val)) return false;
                    GNativeString str = val.getUTF82Native();
                    app_options_set_djbz_suffix(m_appOptions, str.getbuf());
                } else if (token == "warnings") {
                    if (!readValInt("warnings", m_appOptions->warnings)) return false;
                } else {
                    std::cerr << "Wrong options token: \"" << token.getbuf() << "\"\n";
                }
    }
    return true;
}

/*
 in (options ...):
  (default-image     # default image options

   #dpi 300          # if set, use this dpi value for encoding all images
                     # except those that have personal dpi option set.
                     # if not set, use dpi of source image of each page.

   smooth       0    # default smoothing image before processing (off)
   clean        0    # default cleaning image while processing (off)
   erosion      0    # default erosion image after processing (off)
                     # this applies to both images and djbz dictionaries.
  )

 or in (file ...):

  (image             # image options

   #dpi 300          # if set, use this dpi value for encoding this image
                     # if not set use default dpi value if it's set
                     # if default value isn't det, use dpi of source image

   smooth       0    # smoothing image before processing (off)
   clean        0    # cleaning image while processing (off)
   erosion      0    # erosion image after processing (off)
  )
 */

bool
SettingsReader::readImageOptions(struct ImageOptions* opts)
{
    ParsingByteStream& pbs = **m_bs;

    GUTF8String token;
    GUTF8String val;

    while (true) {
        token = pbs.get_token(true).downcase();
        int i = token.search(')');
        if (i > 0) {
            for (int j = token.length()-1; j >= i; j--) {
                pbs.unget(token[j]);
            }
            token = token.substr(0, i).downcase();
        }

        if (!token || token[0] == ')') {
            break;
        }

        if (token == "smooth") {
            if (!readValInt("smooth", opts->smooth)) return false;
        } else if (token == "clean") {
            if (!readValInt("clean", opts->clean)) return false;
        } else if (token == "erosion") {
            if (!readValInt("erosion", opts->erosion)) return false;
        } else if (token == "dpi") {
            if (!readValInt("dpi", opts->dpi, 20, 2000)) return false;
            opts->dpi_specified = 1;
        } else {
            std::cerr << "Wrong image token: \"" << token.getbuf() << "\"\n";
        }
    }


    return true;
}

/* readInputFiles()
(input-files       # A list of input image files
                   # the order is the same as the
                   # the order of pages in document.

                   # Multipage tiff's are expanded to
                   # a set of single page tiffs.

 path/file1        # full filename of the 1st image
                   # will use default image options

 (file             # 2nd image which overwrite default options

   path/file2      # full filename of the 2nd image
   (image          # image settings of the 2nd image
                   # that overrides default settings
     smooth   0
     clean    0
     # etc. as desctibed above
   )
   page       0    # if file is multipage, use only page 0
   page-start 0    # if file is multipage, use pages from 0 to page-end
   page-end   3    # if file is multipage, use pages from page-start to 3
  )

 # etc. for other files. Just write their filename if default settings is fine
 # or include filename in (file ...) list to use page-specific settings.
)
*/

bool
SettingsReader::readFile(struct FileList* file_list, bool ref_only)
{
    ParsingByteStream& pbs = **m_bs;

    GUTF8String orig_token;
    GUTF8String token;
    GUTF8String val;
    GUTF8String filename;

    auto opts = std::unique_ptr< struct ImageOptions, void(*)(void*)>{ nullptr, free };

    int multipage = 0;
    int page_start = 0;
    int page_end = 0;

    while (true) {
        orig_token = pbs.get_token(true);
        token = orig_token.downcase();
        int i = token.search(')');
        if (i > 0) {
            for (int j = token.length()-1; j >= i; j--) {
                pbs.unget(token[j]);
            }
            token = token.substr(0, i).downcase();
        }

        if (!token || token[0] == ')') {
            break;
        }

        if (token[0] == '(') {
            //nested list
            const int len = token.length();
            token = (len > 1) ? token.substr(1, len-1)          // "(file"
                              : pbs.get_token(true).downcase(); // "( file"

            if (token == "image") {

                opts.reset(image_options_create(m_appOptions->default_image_options));
                if (!readImageOptions(opts.get())) {
                    return false;
                }
                continue;
            } else {
                std::cerr << "ERROR: only \"image\" list could be nested in \"file\", but " <<
                             token.getUTF82Native().getbuf() << " found!\n";
                return false;
            }
        }


        if (token == "page") {
            if (!readValInt("page", page_end, 0, 1e6)) return false;
            multipage = 1;
            page_start = page_end;
        } else if (token == "page-start") {
            if (!readValInt("page-start", page_start, 0, 1e6)) return false;
            multipage = 1;
        } else if (token == "page-end") {
            if (!readValInt("page-end", page_end, 0, 1e6)) return false;
            multipage = 1;
        } else {
            filename = orig_token;
        }
    }

    if (!filename.is_valid()) {
        std::cerr << "ERROR: No any files are listed after file token!\n";
        return false;
    }

    if (multipage) {
        if (page_end < page_start) {
            int i = page_start;
            page_start = page_end;
            page_end = i;
        }

        if (!ref_only) {
            file_list_add_filename_with_filter(file_list, filename.getUTF82Native().getbuf(), page_start, page_end, opts.get());
        } else {
            file_list_add_ref_with_filter(file_list, &m_appOptions->file_list, filename.getUTF82Native().getbuf(), page_start, page_end);
        }
    } else {
        if (!ref_only) {
            file_list_add_filename(file_list, filename.getUTF82Native().getbuf(), opts.get());
        } else {
            file_list_add_ref(file_list, &m_appOptions->file_list, filename.getUTF82Native().getbuf());
        }
    }
    return true;
}

bool
SettingsReader::readInputFiles(struct FileList* file_list, bool ref_only)
{
    ParsingByteStream& pbs = **m_bs;

    GUTF8String token;
    int old_cnt = file_list->size;

    while (true) {
        token = pbs.get_token(true);
        int i = token.search(')');
        if (i > 0) {
            for (int j = token.length()-1; j >= i; j--) {
                pbs.unget(token[j]);
            }
            token = token.substr(0, i);
        }

        if (!token || token[0] == ')') {
            break;
        }

        if (token[0] == '(') {
            //nested list
            const int len = token.length();
            token = (len > 1) ? token.substr(1, len-1).downcase() // "(file"
                              : pbs.get_token(true).downcase();   // "( file"

            if (token == "file") {
                if (!readFile(file_list, ref_only)) {
                    return false;
                }
                continue;
            } else {
                std::cerr << "ERROR: only \"file\" list could be nested in \"input-files\" list, but " <<
                             token.getUTF82Native().getbuf() << " found!\n";
                return false;
            }
        }

        if (ref_only) {
            file_list_add_ref(file_list, &m_appOptions->file_list, token.getUTF82Native().getbuf());
        } else {
            file_list_add_filename(file_list, token.getUTF82Native().getbuf(), NULL);
        }
    }

    if (file_list->size - old_cnt == 0) {
        std::cerr << "No any files are listed after input-files token!\n";
        return false;
    }

    return true;
}


/*
 (djbz             # describes a block of page that should belong to
                   # the same djbz dictionary and its settings.
   id         0001 # Mandatory ID of the djbz. Should be unique.
                   # the djbz suffix will be added

   xtension   iff  # overrides default ("djbz") djbz extension, so
                   # the id will be "0001.iff"
   averaging  0    # overrides default averaging (0)
   aggression 100  # overrides default aggression (100)
   classifier 3    # overrides default classifier used to encode this block
   no-prototypes 0 # overides default no-prototypes
   erosion       0 # overrides erosion of glyphs in djbz dictionary
                   # (which is jb2 image by nature)
   (files          # a list of files that should be included in this djbz
                   # files MUST exists in (input-files ...)
                   # the structure is pretty same as in (input-files ...),
                   # but (file ...) lists in (files ...) must not include
                   # (image ...) options as they are provided in (input-files ...)

     path/file1
     (file
      path/file2
      ...
     )
   )
   ...
  )
 */


bool
SettingsReader::readDjbzOptions(struct DjbzOptions *djbz, bool is_defaults)
{
    ParsingByteStream& pbs = **m_bs;

    GUTF8String token; GUTF8String val; int i;

    while (true) {
        token = pbs.get_token(true).downcase();
        i = token.search(')');
        if (i > 0) {
            for (int j = token.length()-1; j >= i; j--) {
                pbs.unget(token[j]);
            }
            token = token.substr(0, i).downcase();
        }

        if (!token || token[0] == ')') {
            break;
        }

        if (token[0] == '(') {
            if (is_defaults) {
                std::cerr << "ERROR: no nested lists are allowed in \"default-djbz\"!\n";
                return false;
            }
            //nested list
            const int len = token.length();
            token = (len > 1) ? token.substr(1, len-1)          // "(files"
                              : pbs.get_token(true).downcase(); // "( files"

            if (token == "files") {
                if (!readInputFiles(&djbz->file_list_ref, true)) {
                    return false;
                }
                continue;
            } else {
                std::cerr << "ERROR: only \"files\" list could be nested in \"djbz\" list, but " <<
                             token.getUTF82Native().getbuf() << " found!\n";
                return false;
            }
        }


        if (token == "id") {
            if (is_defaults) {
                std::cerr << "ERROR: id option isn't allowed in \"default-djbz\"!\n";
                return false;
            }
            if (!readValStr("id", val)) return false;
            GNativeString str = val.getUTF82Native();
            copy_str_alloc(&djbz->id, str.getbuf());
        } else if (token == "averaging") {
            if (!readValInt("averaging", djbz->averaging)) return false;
        } else if (token == "aggression") {
            if (!readValInt("aggression", djbz->aggression, 0, 1000)) return false;
        } else if (token == "classifier") {
            if (!readValInt("classifier", djbz->classifier, 1, 3)) return false;
        } else if (token == "no-prototypes") {
            if (!readValInt("no-prototypes", djbz->no_prototypes)) return false;
        } else if (token == "xtension") {
            if (!readValStr("xtension", val)) return false;
            GNativeString str = val.getUTF82Native();
            copy_str_alloc(&djbz->dict_suffix, str.getbuf());
        } else if (token == "erosion") {
            if (!readValInt("erosion", djbz->erosion)) return false;
        } else {
            std::cerr << "Wrong djbz setting: " << token.getbuf() << "\n";
        }
    }


    if (!is_defaults) {
        // do some verifications
        if (!djbz->id) {
            std::cerr << "id is missing in djbz list!\n";
            return false;
        }

        if (djbz->file_list_ref.size == 0) {
            std::cerr << "djbz list with id \"" << djbz->id << "\" contains no files!\n";
            return false;
        }
    }

    for (int i = 0; i < djbz->file_list_ref.size; i ++) {
        djbz->file_list_ref.files[i]->djbz = djbz;
    }

    return true;
}

bool
SettingsReader::readValInt(const char* name, int& src, int min, int max)
{
    ParsingByteStream& pbs = **m_bs;
    GUTF8String val = pbs.get_token(true);

    int i = val.search(')');
    if (i > 0) {
        for (int j = val.length()-1; j >= i; j--) {
            pbs.unget(val[j]);
        }
        val = val.substr(0, i);
    }

    i = -99999;
    if (!!val && val.is_int()) {
        i = val.toInt();
    }
    if (i < min || i > max) {
        std::cerr << "Wrong \"" << name << "\" value: \"" << val.getbuf() <<
                     "\". Integer from " << min << " to " << max << " is expected.\n";
        return false;
    }
    src = i;
    return true;
}

bool
SettingsReader::readValStr(const char* name, GUTF8String &src)
{
    ParsingByteStream& pbs = **m_bs;
    GUTF8String val = pbs.get_token(true);

    int i = val.search(')');
    if (i > 0) {
        for (int j = val.length()-1; j >= i; j--) {
            pbs.unget(val[j]);
        }
        val = val.substr(0, i);
    }

    if (!val || val.length() == 0) {
        std::cerr << "\"" << name << "\" argument is missing!\n";
        return false;
    }
    src = val;
    return true;
}
