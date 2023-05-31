#ifndef GEFTOOLS_READFLOATTASK_H_
#define GEFTOOLS_READFLOATTASK_H_

#include <map>
#include <unordered_map>
#include "utils.h"
#include "thread_pool.h"
#include "cgef3dParam.h"


class readFloatTask:public ITask
{
public:
    readFloatTask();
    ~readFloatTask();
    void doTask();
protected:
    bool readbuf();
    int cuttail(char *pbuf);
    virtual int getInfo();
    virtual int mergeinfo();

protected:
    int m_buflen = 0;
    char *m_pbuf = nullptr;

    unordered_map<uint32_t, cgef3d_cell*> m_map_cell;
    unordered_map<string, cgef3d_gene*> m_map_gene;
    unordered_map<string, cgef3d_ctype*> m_map_ctype;
    static string m_leftstr;
    static mutex m_readmtx; //读文件锁
    static mutex m_mergemtx; //合并锁
};


// class readFloatTask_5:public readFloatTask
// {
// public:
//     readFloatTask_5(){}
//     ~readFloatTask_5(){}
//     int getInfo();
// };

// class readFloatTask_6:public readFloatTask
// {
// public:
//     readFloatTask_6(){}
//     ~readFloatTask_6(){}
//     int getInfo();
// };
#endif