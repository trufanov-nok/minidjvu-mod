#ifndef SETTINGSREADER_H
#define SETTINGSREADER_H

#include "DjbzSetting.h"

class ParsingByteStream;
template <class T> class GP;

class SettingsReader
{
public:
    SettingsReader(const char* settings_file);
    bool readSettings();
    const DjbzSetting* djbz(int idx) const;
    ~SettingsReader();
private:
    bool readDjbzSetting();
private:
    int m_settings_cnt;
    DjbzSetting** m_DjbzSettings;
    GP<ParsingByteStream>* m_bs;
};

#endif // SETTINGSREADER_H
