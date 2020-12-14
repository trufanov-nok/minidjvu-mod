#ifndef DJBZSETTING_H
#define DJBZSETTING_H

struct DjbzSetting
{
    char * id;
    int files;
    int data_size;
    char** file_list;
};

#ifdef __cplusplus
extern "C" {
#endif
    
struct DjbzSetting* create_DjbzSetting(char * id);
void desitroy_DjbzSetting(struct DjbzSetting* sett);
void reserve_DjbzSetting(struct DjbzSetting* sett, int n);
void add_file_to_DjbzSetting(struct DjbzSetting* sett, char* fname);

int read_settings(const char* fname, struct DjbzSetting** settings);
char* test();

#ifdef __cplusplus
}
#endif

#endif // DJBZSETTING_H
