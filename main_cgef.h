/*
 * @Author: zhaozijian
 * @Date: 2022-02-10 14:53:03
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-05-12 17:05:27
 * @Description: file content
 */
/** @file main_cgef.h
    @brief Main entrance of the geftools cgef command.

    Created by huangzhibo on 2021/12/14.
*/

#ifndef GEFTOOLS__MAIN_CGEF_H_
#define GEFTOOLS__MAIN_CGEF_H_

#include <string>
#include "cgef_writer.h"
#include "cxxopts.h"

using namespace std;

/**
 * @brief Main entrance of the geftools cgef command.
 *
 * Command line arguments for geftools cgef will be resolved here.
 * @param argc
 * @param argv
 * @return
 */
int cgef(int argc, char *argv[]);

/**
 * @brief Function to generate cgef (cell bin GEF) file.
 * @param cgef_file  The output file path of cgef (cell bin GEF).
 * @param bgef_file  The input file path of common bin GEF, only bin1 will be used.
 * @param mask_file  The input file path of mask file.
 * @param block_size The block size of cell dataset in cgef, including two elements, 0: x_block_size, 1: y_block_size
 * @param rand_celltype_num    Number of the randdom cell_type. default : 0.
 * @param verbose    Print the run time of this function. default : false.
 * @return
 */
int GEFTOOLS_API generateCgef(const string& cgef_file,
                 const string& bgef_file,
                 const string& mask_file,
                 const int* block_size,
                 int rand_celltype_num,
                 bool verbose = false);

int GEFTOOLS_API cgem2cgef(const string &strcgem, const string &strcgef, const int* block_size, int rand_celltype_num);

void GEFTOOLS_API AddClusterId4Cgef(const string &input_file,
                                    const string &output_file,
                                    const string &cluster_file);

#endif //GEFTOOLS__MAIN_CGEF_H_
