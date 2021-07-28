#include "SettingsReaderAdapter.h"
#include "SettingsReader.h"
#include "GException.h"

int read_app_options_from_file(const char* fname, struct AppOptions* opts)
{
    if (opts) {
        G_TRY {
            SettingsReader reader(fname, opts);
            return reader.readAllOptions();
        } G_CATCH(exc) {
            exc.perror();
            exit(-2);
        } G_ENDCATCH;
    }

    return 0;
}
