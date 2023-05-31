#include "bin_task.h"


BinTask::BinTask(int bin, const char *ptr):
    m_bin(bin),
    m_geneid(ptr),
    m_maxexp(0)
{
    opts_ = BgefOptions::GetInstance();
}

BinTask::~BinTask()
{
}

void BinTask::bin1task()
{
    std::vector<Expression> &exp_vec = opts_->map_gene_exp_[m_geneid];

    GeneS *pgenes = new GeneS(m_geneid);
    auto *pgeneinfo = new GeneInfo(m_geneid);
    pgeneinfo->vecptr = &exp_vec;

    for (auto& exp : exp_vec)
    {
        if (exp.count > m_maxexp)
        {
            m_maxexp = exp.count;
        }
        if(exp.exon > m_maxexon)
        {
            m_maxexon = exp.exon;
        }
    }
    pgeneinfo->maxexp = m_maxexp;
    pgeneinfo->maxexon = m_maxexon;
    pgenes->vecptr = pgeneinfo->vecptr;
    opts_->m_genes_queue.addqueue(pgenes);
    opts_->m_geneinfo_queue.addqueue(pgeneinfo);
}

void BinTask::bin100task()
{
    unsigned long umicnt = 0;
    std::vector<Expression> &exp_vec = opts_->map_gene_exp_[m_geneid];
    for(auto exp : exp_vec)
    {
        x = exp.x / m_bin;
        y = exp.y / m_bin;
        dnb = x<<32 | y;
        map_dnb[dnb].total+=exp.count;
        map_dnb[dnb].exon+=exp.exon;
        umicnt += exp.count;
    }

    GeneS *pgenes = new GeneS(m_geneid);
    auto *pgeneinfo = new GeneInfo(m_geneid);
    pgeneinfo->vecptr = new std::vector<Expression>;
    pgeneinfo->vecptr->reserve(map_dnb.size());

    pgeneinfo->umicnt = umicnt;
    auto itor_dnb = map_dnb.begin();
    Expression exp{0,0,0};
    for(;itor_dnb!=map_dnb.end();itor_dnb++)
    {
        exp.x = itor_dnb->first>>32;
        exp.y = itor_dnb->first & 0xFFFFFFFF;
        exp.count = itor_dnb->second.total;
        exp.exon = itor_dnb->second.exon;
        pgeneinfo->vecptr->emplace_back(exp);

        if (exp.count > m_maxexp)
            m_maxexp = exp.count;
        if(exp.exon > m_maxexon)
        {
            m_maxexon = exp.exon;
        }
    }
    pgeneinfo->maxexp = m_maxexp;
    pgeneinfo->maxexon = m_maxexon;

    std::sort(pgeneinfo->vecptr->begin(), pgeneinfo->vecptr->end(), 
            [](const Expression &a, const Expression &b){return a.count > b.count;});

    int j = 0;
    int sz = pgeneinfo->vecptr->size() * 0.1;
    unsigned long midcnt = 0;
    auto itor = pgeneinfo->vecptr->begin();
    for(;itor != pgeneinfo->vecptr->end() && j< sz;itor++,j++)
    {
        midcnt += itor->count;
    }
    pgeneinfo->e10 = (midcnt*1.0/umicnt)*100;

    // itor = pgeneinfo->vecptr->begin();
    // midcnt = 0;
    // j = 0;
    // for(;itor != pgeneinfo->vecptr->end();itor++,j++)
    // {
    //     midcnt += itor->count;
    //     if(midcnt*1.0/umicnt > 0.5) 
    //     {
    //         break;
    //     }
    // }

    // sz = pgeneinfo->vecptr->size();
    // pgenedata->c50 = (j*1.0/sz)*100;
    pgenes->vecptr = pgeneinfo->vecptr;
    opts_->m_genes_queue.addqueue(pgenes);
    opts_->m_geneinfo_queue.addqueue(pgeneinfo);
}

void BinTask::othertask()
{
    std::vector<Expression> &exp_vec = opts_->map_gene_exp_[m_geneid];
    for(auto exp : exp_vec)
    {
        x = exp.x / m_bin;
        y = exp.y / m_bin;
        dnb = x<<32 | y;
        map_dnb[dnb].total+=exp.count;
        map_dnb[dnb].exon+=exp.exon;
    }

    GeneS *pgenes = new GeneS(m_geneid);
    GeneInfo *pgeneinfo = new GeneInfo(m_geneid);
    pgeneinfo->vecptr = new std::vector<Expression>;
    pgeneinfo->vecptr->reserve(map_dnb.size());

    auto itor_dnb = map_dnb.begin();
    Expression exp{0,0,0};
    for(;itor_dnb!=map_dnb.end();itor_dnb++)
    {
        exp.x = itor_dnb->first>>32;
        exp.y = itor_dnb->first & 0xFFFFFFFF;
        exp.count = itor_dnb->second.total;
        exp.exon = itor_dnb->second.exon;
        pgeneinfo->vecptr->emplace_back(exp);

        if (exp.count > m_maxexp)
            m_maxexp = exp.count;
        if(exp.exon > m_maxexon)
        {
            m_maxexon = exp.exon;
        }
    }
    pgeneinfo->maxexp = m_maxexp;
    pgeneinfo->maxexon = m_maxexon;
    pgenes->vecptr = pgeneinfo->vecptr;
    opts_->m_genes_queue.addqueue(pgenes);
    opts_->m_geneinfo_queue.addqueue(pgeneinfo);
}

void BinTask::doTask()
{
    if(m_bin == 1)
    {
        bin1task();
    }
    else if(m_bin == 100)
    {
        bin100task();
    }
    else
    {
        othertask();
    }
}