#include "DjbzSetting.h"
#include <stdlib.h>
#include <string.h>

#include "ParsingBytestream.h"
#include "GSmartPointer.h"
#include "GURL.h"
#include "GString.h"

struct DjbzSetting* create_DjbzSetting(char * id)
{
    struct DjbzSetting* sett = (struct DjbzSetting*) calloc(sizeof(struct DjbzSetting), 1);
    int len = strlen(id)+1;
    sett->id = (char*) malloc(len);
    strncpy(sett->id, id, len);
    return sett;
}

void desitroy_DjbzSetting(struct DjbzSetting* sett)
{
    for (int i = 0; i < sett->data_size; ++i) {
        free(sett->file_list[i]);
    }
    if (sett->file_list) {
        free(sett->file_list);
    }
    if (sett->id) {
        free(sett->id);
    }
}

void reserve_DjbzSetting(struct DjbzSetting* sett, int n)
{
    if (n > sett->data_size) {
        char** new_data = (char**) calloc(sizeof(char*), n);
        memcpy(new_data, sett->file_list, sizeof(char*) * sett->data_size);
        free(sett->file_list);
        sett->file_list = new_data;
        sett->data_size = n;
    }
}

void add_file_to_DjbzSetting(struct DjbzSetting* sett, char* fname)
{
    if (sett->data_size <= sett->files) {
        reserve_DjbzSetting(sett, sett->data_size + 10);
    } else {
        sett->file_list[sett->data_size++] = fname;
    }
}


int read_settings(const char* fname, struct DjbzSetting** settings)
{
    int settings_cnt = 0;
    GP<ParsingByteStream>* m_bs;
}

char t[] = {'t', 'e', 's', 't', 0 };

char * test()
{
    return t;
}
