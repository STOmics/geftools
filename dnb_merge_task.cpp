/*
 * @Author: zhaozijian
 * @Date: 2022-02-10 14:53:03
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-04-08 15:22:58
 * @Description: file content
 */
#include "dnb_merge_task.h"
#include <assert.h>

mutex DnbMergeTask::m_mutex;

DnbMergeTask::DnbMergeTask(int cnt, int tid, int binsize):m_genecnt(cnt),m_taskid(tid),
    m_binsize(binsize)
{
    opts_ = BgefOptions::GetInstance();
    y_len = opts_->dnbmatrix_.dnb_attr.len_y;
    int len = opts_->dnbmatrix_.dnb_attr.len_x/opts_->thread_ + 1;
    m_x_low = tid * len;
    m_x_high = m_x_low+len;
}

DnbMergeTask::~DnbMergeTask() = default;

void DnbMergeTask::doTask_nor()
{
    unsigned int maxGene = 0;
    
    unsigned int idx = 0;
    while (idx < m_genecnt)
    {
        GeneS *pgeneinfo = opts_->m_genes_queue.getGeneInfo(idx);
        if(pgeneinfo == nullptr) 
        {
            printf("DnbMergeTask err\n");
            break;
        }
        idx++;

        std::vector<Expression> &exp_vec = *(pgeneinfo->vecptr);
        unsigned long x,y,col;
        if (m_binsize == 1)
        {
            BinStatUS *pmatrix = opts_->dnbmatrix_.pmatrix_us;
            for(auto exp : exp_vec)
            {
                x = exp.x;
                if(x >= m_x_low && x < m_x_high)
                {
                    y = exp.y;
                    col = x*y_len + y;

                    // printf("x %d y %d ylen %d col %d\n", x, y, y_len, col);
                    pmatrix[col].mid_count += exp.count;
                    pmatrix[col].gene_count += 1;
                    if (pmatrix[col].gene_count > maxGene)
                        maxGene = pmatrix[col].gene_count;
                }
            }
        }
        else
        {
            BinStat *pmatrix = opts_->dnbmatrix_.pmatrix;
            for(auto exp : exp_vec)
            {
                x = exp.x;
                if(x >= m_x_low && x < m_x_high)
                {
                    y = exp.y;
                    col = x*y_len + y;
                    pmatrix[col].mid_count += exp.count;
                    pmatrix[col].gene_count += 1;
                    if (pmatrix[col].gene_count > maxGene)
                        maxGene = pmatrix[col].gene_count;
                }
            }
        }
        
    }

    lock_guard<mutex> lock(m_mutex);
    opts_->dnbmatrix_.dnb_attr.max_gene = std::max(opts_->dnbmatrix_.dnb_attr.max_gene, maxGene);
}

void DnbMergeTask::doTask_Exon()
{
    unsigned int maxGene = 0, maxExon = 0;
    
    unsigned int idx = 0;
    while (idx < m_genecnt)
    {
        GeneS *pgeneinfo = opts_->m_genes_queue.getGeneInfo(idx);
        if(pgeneinfo == nullptr) 
        {
            printf("DnbMergeTask err\n");
            break;
        }
        idx++;

        std::vector<Expression> &exp_vec = *(pgeneinfo->vecptr);
        unsigned long x,y,col;
        if (m_binsize == 1)
        {
            BinStatUS *pmatrix = opts_->dnbmatrix_.pmatrix_us;
            unsigned short *exonptr = opts_->dnbmatrix_.pexon16;
            for(auto exp : exp_vec)
            {
                x = exp.x;
                if(x >= m_x_low && x < m_x_high)
                {
                    y = exp.y;
                    col = x*y_len + y;
                    // printf("x %d y %d ylen %d col %d\n", x, y, y_len, col);
                    pmatrix[col].mid_count += exp.count;
                    pmatrix[col].gene_count += 1;
                    exonptr[col] += exp.exon;
                    if (pmatrix[col].gene_count > maxGene)
                        maxGene = pmatrix[col].gene_count;
                    if(exonptr[col] > maxExon)
                    {
                        maxExon = exonptr[col];
                    }
                }
            }
        }
        else
        {
            BinStat *pmatrix = opts_->dnbmatrix_.pmatrix;
            unsigned int *exonptr = opts_->dnbmatrix_.pexon32;
            for(auto exp : exp_vec)
            {
                x = exp.x;
                if(x >= m_x_low && x < m_x_high)
                {
                    y = exp.y;
                    col = x*y_len + y;
                    pmatrix[col].mid_count += exp.count;
                    pmatrix[col].gene_count += 1;
                    exonptr[col] += exp.exon;
                    if (pmatrix[col].gene_count > maxGene)
                        maxGene = pmatrix[col].gene_count;
                    if(exonptr[col] > maxExon)
                    {
                        maxExon = exonptr[col];
                    }
                }
            }
        }
        
    }

    lock_guard<mutex> lock(m_mutex);
    opts_->dnbmatrix_.dnb_attr.max_exon = std::max(opts_->dnbmatrix_.dnb_attr.max_exon, maxExon);
    opts_->dnbmatrix_.dnb_attr.max_gene = std::max(opts_->dnbmatrix_.dnb_attr.max_gene, maxGene);
}

void DnbMergeTask::doTask()
{
    if(opts_->m_bexon)
    {
        doTask_Exon();
    }
    else
    {
        doTask_nor();
    }
}
