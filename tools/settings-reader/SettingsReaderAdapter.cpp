#include "SettingsReaderAdapter.h"
#include "SettingsReader.h"

int read_app_options_from_file(const char* fname, struct AppOptions* opts)
{
    int res = 0;
    if (opts) {
        SettingsReader reader(fname, opts);
        res = reader.readAllOptions();
    }

    return res;
}
