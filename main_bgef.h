/*
 * @Author: zhaozijian
 * @Date: 2022-02-10 14:53:03
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-05-16 14:17:38
 * @Description: file content
 */
/** @file main_cgef.h
    @brief Main entrance of the geftools cgef command.

    Created by huangzhibo on 2021/12/14.
*/

#ifndef GEFTOOLS__MAIN_BGEF_H_
#define GEFTOOLS__MAIN_BGEF_H_

#include <string>
#include <thread>
#include <zlib.h>
#include <queue>
#include "bgef_reader.h"
#include "bgef_writer.h"
#include "cxxopts.h"
#include "bgef_options.h"
#include "utils.h"

using namespace std;

//struct BgefOptions {
//  string input_file;
//  string output_file;
//  vector<int> bin_sizes;
//  int offset_x;
//  int offset_y;
//  int threads;
//  bool verbose;
//};

/**
 * @brief Main entrance of the geftools cgef command.
 *
 * Command line arguments for geftools cgef will be resolved here.
 * @param argc
 * @param argv
 * @return
 */
int bgef(int argc, char *argv[]);

/**
 * @brief Function to generate common bin GEF file(.bgef).
 * @param input_file  The input file path of gem file or bin1 bgef.
 * @param bgef_file   The output file path of common bin GEF (.bgef).
 * @param bin_sizes   Set bin size
 * @param n_thread    Number of thread
 * @param verbose
 * @return
 */
int GEFTOOLS_API generateBgef(const string &input_file,
                 const string &bgef_file,
                 const string &stromics,
                 int n_thread = 8,
                 vector<unsigned int> bin_sizes = vector<unsigned int>(),
                 vector<int> region = vector<int>(),
                 bool verbose = false,
                 bool bstat = true);

void GEFTOOLS_API gem2gef(BgefOptions *opts);

int GEFTOOLS_API mRead(BgefOptions *opts);


void GEFTOOLS_API writednb(BgefOptions *opts, BgefWriter &bgef_writer, int bin);

void GEFTOOLS_API StereoDataToGef(const string &output_file, int binsize, int sz, unsigned long *cellptr);

void GEFTOOLS_API MergeProteinAndRnaMatrices(const string &protein_raw_gef,
                                const string &rna_raw_gef,
                                const string &protein_output_gef,
                                const string &rna_output_gef);

void GEFTOOLS_API Gem2Image(const string &gem_path, const string &tif_path);

#endif //GEFTOOLS__MAIN_BGEF_H_
