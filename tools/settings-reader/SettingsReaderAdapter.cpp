#include "SettingsReaderAdapter.h"
#include "SettingsReader.h"

int read_app_options_from_file(const char* fname, struct AppOptions* opts)
{
    if (opts) {
        SettingsReader reader(fname, opts);
        return reader.readAllOptions();
    }

    return 0;
}
