/** @file main_view.h
    @brief Not yet implemented.

    Created by huangzhibo on 2021/12/14.
*/

#ifndef GEFTOOLS__MAIN_VIEW_H_
#define GEFTOOLS__MAIN_VIEW_H_

#include "utils.h"
#include "cgef_reader.h"
#include "bgef_reader.h"

using namespace std;

/**
 * @brief Parameters parsed from command lines.
 */
struct ViewOptions {
    string input_file;
    string output_gem;
    string output_mask;
    unsigned int region[4];
    vector<string> genes;
    int threads;
    int bin_size;
    bool restrict_region;
    bool force_genes;
    bool exclude; ///< Set the list of genes to exclude, not include.
    bool verbose;
};

/**
 * @brief Not yet implemented
 * @param argc
 * @param argv
 * @return
 */
int view(int argc, char *argv[]);

#endif //GEFTOOLS__MAIN_VIEW_H_
