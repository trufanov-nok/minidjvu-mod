#ifndef SETTINGSREADER_H
#define SETTINGSREADER_H

#include "AppOptions.h"

class ParsingByteStream;
template <class T> class GP;
class GUTF8String;

class SettingsReader
{
public:
    SettingsReader(const char* settings_file, struct AppOptions * appOptions);
    bool readAllOptions();
    ~SettingsReader();
private:
//    bool readToken(GUTF8String& token);
    bool readAppOptions();
    bool readFile(struct FileList *file_list, bool ref_only);
    bool readInputFiles(struct FileList *file_list, bool ref_only);
    bool readDjbzOptions(struct DjbzOptions *djbz, bool is_defaults);
    bool readImageOptions(struct ImageOptions *opts);

    bool readValInt(const char* name, int &src, int min = 0, int max = 1);
    bool readValStr(const char* name, GUTF8String &src);
private:
    GP<ParsingByteStream>* m_bs;
    struct AppOptions* m_appOptions;
};

#endif // SETTINGSREADER_H
