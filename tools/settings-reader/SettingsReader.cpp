#include "SettingsReader.h"
#include "ParsingBytestream.h"
#include "GSmartPointer.h"
#include "GURL.h"
#include "GString.h"
#include <iostream>

SettingsReader::SettingsReader(const char* settings_file):
    m_settings_cnt(0), m_DjbzSettings(nullptr), m_bs(nullptr)
{
    GNativeString str(settings_file);
    if (str.contains("file://") != -1) {
        m_bs = new GP<ParsingByteStream> (
                    ParsingByteStream::create(
                        ByteStream::create(GURL(str), "r")
                        )
                    );
    }
}

const DjbzSetting*
SettingsReader::djbz(int idx) const
{
    return (idx < m_settings_cnt) ? m_DjbzSettings[idx] : nullptr;
}

SettingsReader::~SettingsReader()
{
    if (m_bs) {
        delete m_bs;
    }

    for (int i = 0; i < m_settings_cnt; ++i) {
        delete m_DjbzSettings[0];
    }
    if (m_DjbzSettings) {
        free(m_DjbzSettings);
    }
}

bool
SettingsReader::readSettings()
{
    if (!m_bs) {
        return false;
    }

    GP<ParsingByteStream>& ptr = *m_bs;

    int c = ptr->get_spaces(true);
    if (c == EOF) {
        std::cerr << "Empty file or can't be read\n";
        return false;
    }

    while (c != EOF) {

        if (c!='(') {
            std::cerr << "Syntax error in settings file: expecting '('\n";
            return false;
        }

        if (ptr->get_token() == "djbz") {

        } else {
            std::cerr << "Syntax error in settings file: expecting '(djbz'\n";
            return false;
        }
        c = ptr->get_spaces(true);
    }

    return true;
}

bool
SettingsReader::readDjbzSetting()
{
    if (!m_bs) {
        return false;
    }
    GP<ParsingByteStream>& ptr = *m_bs;

    int c = ptr->get_spaces(true);
    if (c == EOF) {
        std::cerr << "Empty file or can't be read\n";
        return false;
    }

    GUTF8String id = ptr->get_token();
    if (!id.length()) {
        std::cerr << "Djbz id is empty\n";
        return false;
    }

//    DjbzSetting* djbz =
            create_DjbzSetting(id.UTF8ToNative().getbuf());
    return true;
}

