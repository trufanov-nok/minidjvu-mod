#ifndef SETTINGSREADERADAPTER_H
#define SETTINGSREADERADAPTER_H

#include "AppOptions.h"

#ifdef __cplusplus
extern "C" {
#endif

int read_app_options_from_file(const char* fname, struct AppOptions* opts);

#ifdef __cplusplus
}

#endif


#endif // SETTINGSREADERADAPTER_H
