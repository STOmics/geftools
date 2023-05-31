#ifndef GEFTOOLS_BGEFCREATER_H
#define GEFTOOLS_BGEFCREATER_H

#include "utils.h"
#include "gef.h"
#include "thread_pool.h"
#include "opencv2/opencv.hpp"
#include "gene_queue.h"
#include <unordered_map>
//using namespace cv;

struct GEFTOOLS_API gdata
{
    uint32_t geneid;
    vector<uint32_t> vecExp;
};

class GEFTOOLS_API bgefCreater
{
public:
    bgefCreater(int thcnt=8);
    ~bgefCreater();
    void createBgef(const string &strin, int maskbin, const string &strmask, const string &strout);
    void getStereoData(const string &strin, int maskbin, const string &strmask, vector<string> &vec_gene, vector<unsigned long long> &uniq_cells,
                                 vector<unsigned int> &cell_ind, vector<unsigned int> &gene_ind, vector<unsigned int> &count);
private:
    void readbgef(const string &strin);
    void readgem(const string &strin);
    void getmaskgenedata_bgef(vector<Gene> &vgene, vector<Expression> &vgExp, vector<uint8_t> &vexon);
    void getmaskgenedata_gem(vector<Gene> &vgene, vector<Expression> &vgExp, vector<uint8_t> &vexon);
    void writebgef(vector<Gene> &vgene, vector<Expression> &vgExp, vector<uint8_t> &vexon, const string &strout);
    
public:
    bool m_bexon = false;
    int m_bin = 1;
    int m_threadcnt = 8;
    uint32_t m_genencnt = 0;
    uint32_t m_geneexpcnt = 0;
    uint32_t m_maxExp = 0;
    uint32_t m_maxExon = 0;
    uint32_t m_resolution = 0;
    int m_min_x, m_min_y, m_max_x, m_max_y;
    Gene *m_genePtr = nullptr;
    Expression *m_expPtr = nullptr;
    char m_szomics[32]={0};
    cv::Mat m_outimg;
    GefQueue<gdata> m_queue;
    std::unordered_map<std::string, std::vector<Expression>> m_map_gene_exp;
    std::vector<int> m_range = {INT_MAX, 0, INT_MAX, 0};;
    gzFile m_file;
    std::vector<string> m_vecgenename;
    ThreadPool *m_tpoolPtr = nullptr;
};

class GEFTOOLS_API bgefmaskTask:public ITask
{
public:
    bgefmaskTask(uint32_t gid, bgefCreater * ptr):m_gid(gid),m_ptr(ptr)
    {}
    ~bgefmaskTask(){}
    void doTask()
    {
        uint32_t count = m_ptr->m_genePtr[m_gid].count;
        uint32_t offset = m_ptr->m_genePtr[m_gid].offset;
        Expression *ptr = m_ptr->m_expPtr+offset;
        gdata *pdata = new gdata();
        pdata->geneid = m_gid;
        int x = 0, y = 0, ret = 0;
        for(uint32_t i=0;i<count;i++)
        {
            x = (ptr[i].x/m_ptr->m_bin)*m_ptr->m_bin;
            y = (ptr[i].y/m_ptr->m_bin)*m_ptr->m_bin;
            ret = m_ptr->m_outimg.at<uchar>(y, x);
            if(ret)
            {
                pdata->vecExp.push_back(i+offset);//绝对位置
            }
        }
        m_ptr->m_queue.addqueue(pdata);
    }
private:
    uint32_t m_gid = 0;
    bgefCreater *m_ptr = nullptr;
};

class GEFTOOLS_API gemmaskTask:public ITask
{
public:
    gemmaskTask(uint32_t gid, bgefCreater * ptr):
    m_gid(gid),m_ptr(ptr)
    {}
    ~gemmaskTask(){}
    void doTask()
    {
        gdata *pdata = new gdata();
        pdata->geneid = m_gid;

        string &str = m_ptr->m_vecgenename[m_gid];
        auto &vecexp = m_ptr->m_map_gene_exp[str];
        uint32_t i = 0;
        int x = 0, y = 0, ret = 0;
        for(Expression &exp : vecexp)
        {
            x = (exp.x/m_ptr->m_bin)*m_ptr->m_bin;
            y = (exp.y/m_ptr->m_bin)*m_ptr->m_bin;
            ret = m_ptr->m_outimg.at<uchar>(y, x);
            if(ret)
            {
                pdata->vecExp.push_back(i);//相对位置
            }
            i++;
        }
        m_ptr->m_queue.addqueue(pdata);
    }
private:
    uint32_t m_gid = 0;
    bgefCreater *m_ptr = nullptr;
};


#endif