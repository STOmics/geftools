/*
 * @Author: zhaozijian
 * @Date: 2022-02-10 14:53:03
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-05-16 14:15:54
 * @Description: file content
 */

#ifndef GENETOH5_COMMANDPARSE_H
#define GENETOH5_COMMANDPARSE_H

#include <unordered_map>
#include <zlib.h>
#include "gef.h"
#include "gene_info_queue.h"
#include "gene_queue.h"


class GEFTOOLS_API BgefOptions {
private:
    BgefOptions(){};
    ~BgefOptions(){};
public:
    static BgefOptions *GetInstance()
    {
        static BgefOptions instance;
        return &instance;
    }

    int thread_ = 8; //设置取bin线程数
//    int m_thread_dnb = 8; //设置dnbmerge线程数
    bool reverse_ = false; // true: gef to gem, false: gem to gef
    bool verbose_ = false;
    bool m_bexon = false;
    int m_stattype = 0; //0:不生成stat 1:只生成stat 2:生成stat和对应bin数据 

    string input_file_;
    string output_file_;
    vector<unsigned int> bin_sizes_;
    std::vector<int> region_;

    std::unordered_map<std::string, std::vector<Expression>> map_gene_exp_;
    //std::vector<GeneErank> vec_bin100_;
    std::vector<GeneStat> m_vec_bin100;
    unsigned long total_umicnt_ = 0;

    gzFile infile_;

    mutex dnbmtx_;
    DnbMatrix dnbmatrix_;
    std::vector<int> range_ = {INT_MAX, 0, INT_MAX, 0};
    GeneInfoQueue m_genes_queue;
    GefQueue<GeneInfo> m_geneinfo_queue;
    std::vector<Expression> expressions_;
    std::vector<Gene> genes_;
    int offset_x_ = 0; // offset of coordinate
    int offset_y_ = 0;
    string m_stromics; //物种组别
};

#endif //GENETOH5_COMMANDPARSE_H