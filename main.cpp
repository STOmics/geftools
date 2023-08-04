#include <iostream>

#include "main_bgef.h"
#include "main_cgef.h"
#include "main_view.h"

static int usage() {
    cerr << endl;
    cerr << "Program: geftools (Tools for manipulating GEFs)" << endl;
    cerr << "Version: " << GEFVERSION[0] << "." << GEFVERSION[1] << "." << GEFVERSION[2] << endl;
    cerr << "Usage:   geftools <command> [options]\n" << endl;
    cerr << "Command: bgef          Generate common bin GEF(.bgef) according to gem file or bin1 GEF" << endl;
    cerr << "         cgef          Generate cell bin GEF(.cgef) according to common bin GEF and mask file" << endl;
    cerr << "         view          View GEF, generate gem file" << endl;
    cerr << "\nNote: Please report issues at https://github.com/BGIResearch/geftools/issues" << endl;
    return 1;
}

int main(int argc, const char *argv[]) {
    time_t prev;
    time(&prev);

    int ret;
    if (argc < 2) return usage();

    for (int i = 0; i < argc; i++) {
        if (memcmp(argv[i], "-w", 2) == 0) {
            errorCode::isInSAWFlow = true;
        }
    }

    if (strcmp(argv[1], "bgef") == 0)
        ret = bgef(argc - 1, const_cast<char **>(argv + 1));
    else if (strcmp(argv[1], "cgef") == 0)
        ret = cgef(argc - 1, const_cast<char **>(argv + 1));
    else if (strcmp(argv[1], "view") == 0)
        ret = view(argc - 1, const_cast<char **>(argv + 1));
    else {
        log_error << errorCode::E_INVALIDPARAM << "[main] unrecognized command : " << argv[1];
        return 1;
    }
    return ret;
}
