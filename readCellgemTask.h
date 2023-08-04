/*
 * @Author: zhaozijian
 * @Date: 2022-03-25 14:18:03
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-05-12 13:45:33
 * @Description: file content
 */
#ifndef GEFTOOLS_READCELLGEMTASK_H_
#define GEFTOOLS_READCELLGEMTASK_H_

#include <sstream>
#include <unordered_map>

#include "cgefUtil.h"
#include "gef.h"
#include "thread_pool.h"
#include "utils.h"
// #include <sys/time.h>

class readCellgemTask : public ITask {
  public:
    readCellgemTask();
    ~readCellgemTask();
    void doTask();

  protected:
    bool readbuf();
    int cuttail(char *pbuf);
    virtual int getInfo();
    virtual int mergeinfo();

  protected:
    int m_buflen = 0;
    char *m_pbuf = nullptr;
    unordered_map<int, cgef_cell *> m_map_cell;
    unordered_map<string, cgef_gene *> m_map_gene;  // 保存一个基因出现的细胞情况
    map<string, int> m_map_gene_id;

    static string m_leftstr;
    static mutex m_readmtx;   // 读文件锁
    static mutex m_mergemtx;  // 合并锁

    int m_min_x = INT_MAX, m_min_y = INT_MAX, m_max_x = 0, m_max_y = 0;
};

// geneID xPos yPos UMICount
class readCellgemTask_raw : public readCellgemTask {
  public:
    readCellgemTask_raw() {};
    ~readCellgemTask_raw() {};
    int getInfo();
    int mergeinfo();

  private:
    unordered_map<string, bgef_gene *> m_map_bgene;
};

class readCellgemTask_cell : public readCellgemTask {
  public:
    readCellgemTask_cell(bool exon) : m_bexon(exon) {};
    ~readCellgemTask_cell() {};
    int getInfo();

  private:
    int getdata();
    int getdata_exon();

  private:
    bool m_bexon = false;
};

#endif