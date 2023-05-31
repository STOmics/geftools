/*
 * @Author: zhaozijian
 * @Date: 2022-02-10 14:53:03
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-05-16 14:00:03
 * @Description: file content
 */

#ifndef GENETOH5_READTASK_H
#define GENETOH5_READTASK_H

#include <unordered_map>
#include "gef.h"
#include "utils.h"
#include "thread_pool.h"

using namespace std;

typedef struct
{
    int readlen; //想要读取的长度
    int reallen; //实际读取的长度
}RLen;


class ReadTask:public ITask
{
public:
    ReadTask(bool bexon, gzFile file, std::vector<int> &range, std::unordered_map<std::string, std::vector<Expression>> &map_gene_exp);
    ~ReadTask();
    void doTask();
private:
    void readbuf(RLen &rlen);
    int cuttail(char *pbuf);
    int getGeneInfo();
    int mergeGeneinfo();
    int getGeneInfo_exon();
private:
    bool m_bexon = false;
    int m_buflen = 0;
    int min_x = INT_MAX, min_y = INT_MAX, max_x = 0, max_y = 0;
    char *m_pbuf = nullptr;
    unordered_map<string, vector<Expression>> m_map_gege;

    gzFile m_file;
    std::vector<int> &m_range;
    std::unordered_map<std::string, std::vector<Expression>> &m_map_gene_exp;

    static string m_leftstr;
    static mutex m_readmtx; //读文件锁
    static mutex m_mergemtx; //合并锁
};


#endif //GENETOH5_READTASK_H