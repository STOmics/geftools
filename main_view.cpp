//
// Created by huangzhibo on 2021/12/16.
//

#include "main_view.h"
#include "cxxopts.h"
#include "geftogem.h"

int view(int argc, char *argv[]) {
    cxxopts::Options options("geftools view",
                             "About:  Show the contents of cell bin GEF\n");
    //TODO support restrict gene_list and region for bGEF 
    options
        .set_width(120)
        .add_options()
            ("i,input-file", "Input bGEF/cGEF file [request]", cxxopts::value<std::string>(), "FILE")
            ("o,output-gem", "Output gem file ", cxxopts::value<std::string>()->default_value("stdout"), "FILE")
            ("d,exp_data", "Input bgef for cgem", cxxopts::value<std::string>()->default_value(""), "FILE")
            ("m,mask-file", "input mask file ", cxxopts::value<std::string>(), "FILE")
            // ("r,region", "Restrict to a rectangular region. The region is represented by the comma-separated list "
            //              "of two vertex coordinates (minX,maxX,minY,maxY). just support cGEF.",
            //              cxxopts::value<std::string>()->default_value(""), "STR")
            // ("g,genes", "Comma separated list of genes to include (or exclude with \"^\" prefix). just support cGEF.",
            //         cxxopts::value<std::string>(), "[^]STR")
            // ("G,genes-file", "File of genes to include (or exclude with \"^\" prefix)). just support cGEF.",
            //         cxxopts::value<std::string>(), "[^]FILE")
            // ("force-genes", "Only warn about unknown subset genes, just support cGEF.",
            //         cxxopts::value<bool>()->default_value("false"))
            ("b,bin-size", "Set bin size for bgef file, just support bGEF.", cxxopts::value<int>()->default_value("1"), "INT")
            ("s,serial-number", "Serial number [request]", cxxopts::value<std::string>(), "STR")
//            ("t,threads", "number of threads", cxxopts::value<int>()->default_value("1"), "INT")
            //("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
            ("e,exon", "whether or not output exon", cxxopts::value<int>()->default_value("1"), "INT")
            // ("w,errorCode-file", "is in saw flow", cxxopts::value<bool>()->default_value("false"))
            ("help", "Print help");

    auto result = options.parse(argc, argv);

    if (argc <= 1 || result.count("help"))
    {
        std::cerr << options.help() << std::endl;
        reportErrorCode2File(errorCode::E_INVALIDPARAM, "missing params");
        exit(1);
    }

    // if (result.count("errorCode-file") == 1) {
    //     errorCode::isInSAWFlow = result["errorCode-file"].as<bool>();
    // }

    if (result.count("input-file") != 1){
        std::cerr << "[ERROR] The -i,--input-file parameter must be given correctly.\n" << std::endl;
        std::cerr << options.help() << std::endl;
        reportErrorCode2File(errorCode::E_INVALIDPARAM, 
                            "[ERROR] The -i,--input-file parameter must be given correctly."); 
        exit(1);
    }
    if (result.count("serial-number") != 1){
        std::cerr << "[ERROR] The -s,--serial-number parameter must be given correctly.\n" << std::endl;
        std::cerr << options.help() << std::endl;
        reportErrorCode2File(errorCode::E_INVALIDPARAM, 
                            "[ERROR] The -s,--serial-number parameter must be given correctly."); 
        exit(1);
    }

    bool boutexon = result["exon"].as<int>();
    string strin = result["input-file"].as<string>();
    string snstr = result["serial-number"].as<string>();
    string strout = result["output-gem"].as<string>();

    geftogem gem(strout, snstr, boutexon);
    if(is_bgef(strin))
    {
        if(result.count("mask-file") != 1)
        {
            int bin_size = result["bin-size"].as<int>();
            gem.bgeftogem(strin, bin_size);
        }
        else
        {
            string strmask = result["mask-file"].as<string>();
            gem.bgeftocgem(strmask, strin);
        }
    }
    else
    {
        if (result.count("exp_data") == 1)
        {
            string strexp = result["exp_data"].as<string>();
            gem.cgeftogem(strin, strexp);
        }
        else
        {
            std::cerr << "[ERROR] The -d,--exp_data parameter must be given correctly.\n" << std::endl;
            std::cerr << options.help() << std::endl;
            reportErrorCode2File(errorCode::E_INVALIDPARAM, 
                            "[ERROR] The -d,--exp_data parameter must be given correctly."); 
            exit(1);
        }
    }
    
    return 0;
}
