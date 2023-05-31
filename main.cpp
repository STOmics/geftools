#include <iostream>
#include "main_bgef.h"
#include "main_cgef.h"
#include "main_view.h"

//using namespace cv;

static int usage()
{
    cerr << endl;
    cerr << "Program: geftools (Tools for manipulating GEFs)" << endl;
    cerr << "Version: " << GEFVERSION[0]<<"."<<GEFVERSION[1]<<"."<<GEFVERSION[2] << endl;
//    cerr << "Contact: Huang Zhibo <huangzhibo@genomics.cn>\n" << endl;
    cerr << "Usage:   geftools <command> [options]\n" << endl;
    cerr << "Command: bgef          Generate common bin GEF(.bgef) according to gem file or bin1 GEF" << endl;
    cerr << "         cgef          Generate cell bin GEF(.cgef) according to common bin GEF and mask file" << endl;
    cerr << "         view          View GEF, generate gem file" << endl;
//    fprintf(stderr, "         h5ls          scan h5 file\n");
//    fprintf(stderr, "         mask          manipulating mask file\n");
    cerr << "\nNote: Please report issues at https://github.com/BGIResearch/geftools/issues" << endl;
    return 1;
}

int main(int argc, const char* argv[]){
    time_t prev;
    time(&prev);

    // MergeProteinAndRnaMatrices("/jdfssz1/ST_BIGDATA/Stereomics_TestData/debug_tmp/stereotools_test/20221112_cellbin/02.count/SS200000135TL_D1.raw.gef",
    // "/jdfssz1/ST_BIGDATA/Stereomics_TestData/debug_tmp/stereotools_test/20221113_1g/02.count/SS200000135TL_D1.raw.gef",
    // "./protein.gef",
    // "./rna.gef");
    // return 0;
    
    int ret;
    if (argc < 2) return usage();

    // cxxopts::Options options("geftools main", "About: is in saw flow\n");
    // options
    // .set_width(120)
    // .add_options()
    // ("w,errorCode-file", "is in saw flow", cxxopts::value<bool>()->default_value("false"));
    // auto result = options.parse(argc, argv);
    // if (result.count("errorCode-file") == 1) {
    //     errorCode::isInSAWFlow = result["errorCode-file"].as<bool>();
    // }

    for(int i=0;i<argc;i++)
    {
        if(memcmp(argv[i], "-w", 2) == 0)
        {
            errorCode::isInSAWFlow = true;
        }
    }
    
    if (strcmp(argv[1], "bgef") == 0) ret = bgef(argc-1, const_cast<char **>(argv + 1));
    else if (strcmp(argv[1], "cgef") == 0) ret = cgef(argc-1, const_cast<char **>(argv + 1));
    else if (strcmp(argv[1], "view") == 0) ret = view(argc-1, const_cast<char **>(argv + 1));
    else {
        cerr << "[main] unrecognized command " << argv[1] << endl;
        char errMsg[32]={0};
        sprintf(errMsg, "[main] unrecognized command : %s", argv[1]);
        reportErrorCode2File(errorCode::E_INVALIDPARAM, errMsg);
        return 1;
    }
    return ret;
}
