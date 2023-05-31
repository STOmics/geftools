
#ifndef GENETOH5_BINTASK_H
#define GENETOH5_BINTASK_H

#include <map>
#include <algorithm>
#include "thread_pool.h"
#include "bgef_options.h"

struct GEFTOOLS_API midcnt
{
    unsigned int total;
    unsigned int exon;
};


class GEFTOOLS_API BinTask:public ITask
{
public:
    BinTask(int bin, const char *ptr);
    ~BinTask();
    void doTask();
private:
    void bin1task();
    void bin100task();
    void othertask();
private:
    int m_bin;
    const char *m_geneid;
    BgefOptions *opts_;
    std::map<unsigned long long, midcnt> map_dnb;
    unsigned long long x, y, dnb;
    unsigned int m_maxexp = 0;
    unsigned int m_maxexon = 0;
};

#endif