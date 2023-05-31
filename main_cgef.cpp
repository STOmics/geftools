// Created by huangzhibo on 2021/12/14.
//

#include "main_cgef.h"
#include "opencv2/opencv.hpp"
#include "cgefParam.h"
#include "cgefCellgem.h"
#include "cellAdjust.h"
#include "cgef_reader.h"
#include "cgef3d.h"
int te1(const char *p1, const char *p2)
{
    cellAdjust cg;
    vector<sapBgefData> vecdata;
    vector<vector<int>> vecpos;

    vector<int> v1{12988,10307,12993,10307,12993,10313,12992,10313,12992,10315,12991,10315,12990,10316,12989,10316,12989,10317,12986,10317,12985,10316,12983,10316,12983,10315,12981,10315,12981,10314,12979,10314,12978,10313,12986,10306,12988,10306};
    vector<int> v2{12985,10309,12982,10314,12987,10314,12987,10313,12988,10313,12988,10310,12987,10310,12987,10309};
    vecpos.emplace_back(std::move(v1));
    vecpos.emplace_back(std::move(v2));


    cg.getSapRegion(p1, 1, 10, vecpos, vecdata);
    printf("cnt:%d\n", vecdata.size());
    return 0;
}

int cgef(int argc, char *argv[]) {
    cxxopts::Options options("geftools cgef",
                       "About:  Generate cell bin GEF (.cgef) according to"
                       " common bin GEF (.bgef) file and mask file\n");
    options
    .set_width(120)
    .add_options()
    ("i,input-file", "input GEF file [request]", cxxopts::value<std::string>(), "FILE")
    ("m,mask-file", "input mask file [request]", cxxopts::value<std::string>(), "FILE")
    ("o,output-file", "output cell bin GEF file (.cgef) [request]", cxxopts::value<std::string>(), "FILE")
    ("b,block", "Pre block size", cxxopts::value<std::string>()->default_value("256,256"), "FILE")
    ("r,rand-celltype", "number of random cell type", cxxopts::value<int>()->default_value("0"), "INT")
    ("t,threads", "number of threads", cxxopts::value<int>()->default_value("8"), "INT")
    ("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
    ("g,raw-gem", "raw gem file", cxxopts::value<std::string>(), "FILE")
    ("p,patch", "Create 3d group patch", cxxopts::value<int>()->default_value("0"))
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

    if (result.count("mask-file") != 1){
        cgefParam::GetInstance()->m_maskstr = "";
    }
    else
    {
        cgefParam::GetInstance()->m_maskstr = result["mask-file"].as<string>();
    }

    if (result.count("output-file") != 1){
        cgefParam::GetInstance()->m_outputstr = "";
    }
    else
    {
        cgefParam::GetInstance()->m_outputstr = result["output-file"].as<string>();
    }

    if(result.count("raw-gem") != 1)
    {
        cgefParam::GetInstance()->m_rawgemstr = "";
    }
    else
    {
        cgefParam::GetInstance()->m_rawgemstr = result["raw-gem"].as<string>();
    }

    int rand_celltype_num = result["rand-celltype"].as<int>();
    cgefParam::GetInstance()->m_inputstr = result["input-file"].as<string>();
    cgefParam::GetInstance()->m_threadcnt = result["threads"].as<int>();
    
    vector<string> block_size_tmp = split(result["block"].as<string>(), ',');
    if(block_size_tmp.size() != 2){
        std::cerr << "[ERROR] The -b,--block parameter must be given correctly.\n" << std::endl;
        std::cerr << options.help() << std::endl;
        reportErrorCode2File(errorCode::E_INVALIDPARAM, 
                            "[ERROR] The -b,--block parameter must be given correctly."); 
        exit(1);
    }
    cgefParam::GetInstance()->m_block_size[0] = static_cast<int>(strtol(block_size_tmp[0].c_str(), nullptr, 10));
    cgefParam::GetInstance()->m_block_size[1] = static_cast<int>(strtol(block_size_tmp[1].c_str(), nullptr, 10));

    int patch = result["patch"].as<int>();
    if(patch == 1)
    {
        cgef3dParam::GetInstance()->m_threadcnt = cgefParam::GetInstance()->m_threadcnt;
        cgef3d c3d;
        c3d.writeCgef(cgefParam::GetInstance()->m_inputstr, 
                    cgefParam::GetInstance()->m_rawgemstr,
                    cgefParam::GetInstance()->m_maskstr,
                    cgefParam::GetInstance()->m_outputstr);
    }
    else if(patch == 0)
    {
        generateCgef(cgefParam::GetInstance()->m_outputstr,
                    cgefParam::GetInstance()->m_inputstr, 
                    cgefParam::GetInstance()->m_maskstr,
                    cgefParam::GetInstance()->m_block_size,
                    rand_celltype_num);
    }
    else if(patch == 2)
    {
        cgem2cgef(cgefParam::GetInstance()->m_inputstr, cgefParam::GetInstance()->m_outputstr, 
                    cgefParam::GetInstance()->m_block_size, rand_celltype_num);
    }

    return 0;
}

int generateCgef(const string &cgef_file,
                 const string &bgef_file,
                 const string &mask_file,
                 const int* block_size,
                 int rand_celltype_num,
                 bool verbose) {
    unsigned long cprev=clock();
    CgefWriter cgef_writer(verbose);
    cgef_writer.setOutput(cgef_file);
    cgef_writer.setRandomCellTypeNum(rand_celltype_num);

    cgefCellgem cgem;
    cgem.writeFile(&cgef_writer, mask_file, bgef_file);

    if(verbose) printCpuTime(cprev, "generateCgef");
    return 0;
}

int cgem2cgef(const string &strcgem, const string &strcgef, const int* block_size, int rand_celltype_num)
{
    cgefParam::GetInstance()->m_block_size[0] = block_size[0];
    cgefParam::GetInstance()->m_block_size[1] = block_size[1];
    CgefWriter cgef_writer(false);
    cgef_writer.setOutput(strcgef);
    cgef_writer.setRandomCellTypeNum(rand_celltype_num);

    cgefCellgem cgef;
    cgef.cgem2cgef(&cgef_writer, strcgem);
    return 0;
}

std::map<unsigned int, unsigned short> ParseClusterFile(
    const string &cluster_file) {
    vector<string> cluster_info = readLines(cluster_file);
    if (cluster_info.empty()) {
        reportErrorCode2File(errorCode::E_INVALIDPARAM,
                             "[ERROR] Input file is not cellbin file.");
    }

    std::map<unsigned int, unsigned short> cell_cluster;
    for (uint64_t i = 1; i < cluster_info.size(); i++) {
        vector<string> cluster_id = split(cluster_info[i], '\t');

        if (!cluster_id.empty()) {
          unsigned int cellId = std::stoi(cluster_id[0]);
          cell_cluster[cellId] = std::stoi(cluster_id[1]);
        } else {
          reportErrorCode2File(errorCode::E_INVALIDPARAM,
                               "[ERROR] Input file is not cellbin file.");
        }
    }
    return cell_cluster;
}

void AddClusterId4Cgef(const string &input_file, const string &output_file,
                       const string &cluster_file) {
    std::cout << "parse cell_cluster file. " << std::endl;
    if (!copyFile(input_file, output_file)) {
        reportErrorCode2File(errorCode::E_INVALIDPARAM,
                             "[ERROR] Can not write cellbin file.");
        return;
    }
    auto cluster_info = ParseClusterFile(cluster_file);

    std::cout << "write cluster id to cgef. " << std::endl;
    CgefReader cgef_info = CgefReader(input_file);
    CellData *cell_array = cgef_info.getCell();
    for (uint32_t i = 0; i < cgef_info.getCellNum(); i++) {
        if (cluster_info.find(cell_array[i].id) != cluster_info.end()) {
          cell_array[i].cluster_id = cluster_info[cell_array[i].id];
        } else {
          cell_array[i].cluster_id = 32;
        }
    }

    hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_libver_bounds(fapl_id, H5F_LIBVER_V18, H5F_LIBVER_LATEST);
    H5Pset_fclose_degree(fapl_id, H5F_CLOSE_STRONG);
    hid_t file_id = H5Fopen(output_file.c_str(), H5F_ACC_RDWR, fapl_id);
    hid_t group_id = H5Gopen(file_id, "/cellBin", H5P_DEFAULT);
    hid_t cell_dataset_id = H5Dopen(group_id, "cell", H5P_DEFAULT);
    hid_t memtype = getMemtypeOfCellData();
    H5Dwrite(cell_dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
             cell_array);

    H5Tclose(memtype);
    H5Dclose(cell_dataset_id);
    H5Gclose(group_id);
    H5Fclose(file_id);
}
