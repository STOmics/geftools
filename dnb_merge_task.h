
#ifndef GENETOH5_DNBMERGETASK_H
#define GENETOH5_DNBMERGETASK_H

#include "gef.h"
#include "thread_pool.h"
#include "gene_info_queue.h"
#include "bgef_options.h"

class DnbMergeTask:public ITask
{
public:
    DnbMergeTask(int cnt, int tid, int binsize);
    ~DnbMergeTask();
    void doTask();
private:
    void doTask_nor();
    void doTask_Exon();

private:
    BgefOptions * opts_ = nullptr ;
    int m_genecnt = 0;
    int m_taskid = 0;
    int m_binsize = 0;
    int m_x_low = 0;
    int m_x_high = 0;
    int y_len = 0; //矩阵Y轴最大值和最小值之差
    static mutex m_mutex;
};

#endif