/*
 * @Author: zhaozijian
 * @Date: 2022-03-25 14:15:30
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-05-21 14:31:48
 * @Description: file content
 */

#include "cgefCellgem.h"
#include "cgefParam.h"
#include "readCellgemTask.h"
#include "getcellbinTask.h"
#include "timer.h"
#include <functional>
#include <set>


cgefCellgem::cgefCellgem()
{
    m_thpoolPtr = new ThreadPool(cgefParam::GetInstance()->m_threadcnt);
}

cgefCellgem::~cgefCellgem()
{
    delete m_thpoolPtr;
}

void cgefCellgem::gemPreAnalysis(const string &strmask, const string &strinput)
{
    if(H5Fis_hdf5(strinput.c_str()))
    {
        cgefParam::GetInstance()->m_intype = INPUTTYPE_BGEF_MASK;
        return;
    }

    cgefParam::GetInstance()->m_infile = gzopen(strinput.c_str(), "r");
    gzbuffer(cgefParam::GetInstance()->m_infile, READLEN);

    char buf[128] = {0};
    while (1)
    {
        gzgets(cgefParam::GetInstance()->m_infile, buf, 128);
        if(memcmp(buf, "geneID", 6) == 0)
        {
            break;
        }
    }
    
    int col = 1;
    int i = 0;
    while (buf[i] != 0)
    {
        if(buf[i] == '\t')
        {
            col++;
        }
        i++;
    }
    printf("%s %d\n", buf, col);
}

void cgefCellgem::cgem2cgef(CgefWriter *cwptr, const string &strin)
{
    m_cgefwPtr = cwptr;
    cgefParam::GetInstance()->m_infile = gzopen(strin.c_str(), "r");
    gzbuffer(cgefParam::GetInstance()->m_infile, READLEN);

    char buf[128] = {0};
    while (1)
    {
        gzgets(cgefParam::GetInstance()->m_infile, buf, 128);
        if(memcmp(buf, "geneID", 6) == 0)
        {
            break;
        }
    }
    
    int col = 1;
    int i = 0;
    while (buf[i] != 0)
    {
        if(buf[i] == '\t')
        {
            col++;
        }
        i++;
    }
    printf("%s %d\n", buf, col);

    if(col == 6)
    {
        m_bexon = true;
    }

    for(int i=0;i<cgefParam::GetInstance()->m_threadcnt;i++)
    {
        readCellgemTask *rtask = new readCellgemTask_cell(m_bexon);
        m_thpoolPtr->addTask(rtask);
    }
    m_thpoolPtr->waitTaskDone();
    gzclose(cgefParam::GetInstance()->m_infile);

    getCelldata_cgem();
    writeCell_cgem();
    writeGene_cgem();
    writeAttr();
}

void cgefCellgem::writeFile(CgefWriter *cwptr, const string &strmask, const string &strinput)
{
    m_cgefwPtr = cwptr;
    gemPreAnalysis(strmask, strinput);
    switch (cgefParam::GetInstance()->m_intype)
    {
    case INPUTTYPE_BGEF_MASK:
        readBgef_new(strinput);
        readmask_new(strmask);
        writeAttr();
        getCell();
        writeCell_new();
        writeGene_new();
        break;

    default:
        break;
    }
}

void cgefCellgem::writeAttr()
{
    CellBinAttr cell_bin_attr = {
            /*.version = */2,
            /*.resolution = */cgefParam::GetInstance()->m_resolution,
            /*.offsetX = */cgefParam::GetInstance()->m_min_x,
            /*.offsetY = */cgefParam::GetInstance()->m_min_y,
            /*.omics = */m_stromics
    };
    m_cgefwPtr->storeAttr(cell_bin_attr);
}

void cgefCellgem::addCellborder(int cx, int cy, vector<short> &vec_border, uint32_t idx)
{
    vector<cv::Point> &vborder = m_contours[idx];
    int i=0;
    int sz = vborder.size();
    if(sz > BORDERCNT)
    {
        vector<cv::Point> tmpborder;
        double epsilon = 0.01 * arcLength(vborder, true);
        approxPolyDP(vborder, tmpborder, epsilon, true);

        sz = tmpborder.size();
        for(;i<sz;i++)
        {
            vec_border.emplace_back(tmpborder[i].x - cx);
            vec_border.emplace_back(tmpborder[i].y - cy);
        }
    }
    else
    {
        for(;i<sz;i++)
        {
            vec_border.emplace_back(vborder[i].x - cx);
            vec_border.emplace_back(vborder[i].y - cy);
        }
    }
    
    for(;i<BORDERCNT;i++) //不足补0
    {
        vec_border.emplace_back(SHRT_MAX);
        vec_border.emplace_back(SHRT_MAX);
    }
}

// void cgefCellgem::writeCell()
// {
//     timer st(__FUNCTION__);
//     unsigned int cid = 0, gid = 0, offcnt = 0; 
//     unsigned short gene_count, exp_count, dnb_count, area, cell_type_id;
//     uint16_t maxExpmid = 0;
//     int cx, cy; //细胞质点
//     vector<short> vec_border;
//     vec_border.reserve(m_maskcellnum*2*BORDERCNT);
//     vector<unsigned int> vec_blkidx;
//     vec_blkidx.reserve(m_blocknum+1);
//     vector<unsigned int> vec_cellLabel;
//     vec_cellLabel.reserve(m_maskcellnum);

//     for(vector<celldata> &vec : m_vec_veccell)
//     {
//         int cnt = 0;
//         for(celldata & cdata : vec)
//         {
//             cgef_cell *cellptr = cgefParam::GetInstance()->m_map_cell[cdata.l_idx];
//             if(cellptr == nullptr) continue;
 
//             m_hash_clabel2cid.emplace(cdata.l_idx, cid);
//             vec_cellLabel.emplace_back(cdata.l_idx);
            
//             cx = m_centroids.at<double>(cdata.l_idx,0);
//             cy = m_centroids.at<double>(cdata.l_idx,1);
//             addCellborder(cx, cy, vec_border, cdata.c_idx);

//             gene_count = cellptr->m_map_cellexp.size();
//             exp_count = cellptr->expcnt;
//             dnb_count = cellptr->dnbcnt;
//             auto itor = cellptr->m_map_cellexp.begin();
//             for(;itor != cellptr->m_map_cellexp.end();itor++)
//             {
//                 gid = m_hash_gname2gid[itor->first];
//                 m_cgefwPtr->cell_exp_list_.emplace_back(gid, itor->second);
//                 maxExpmid = std::max(maxExpmid, itor->second);
//             }

//             area = m_stats.at<int>(cdata.l_idx, CC_STAT_AREA);
//             cell_type_id = m_cgefwPtr->random_cell_type_num_ == 0 ? 0 : rand()%(m_cgefwPtr->random_cell_type_num_ + 1);
//             CellData cell = {
//                     cid,
//                     cx,
//                     cy,
//                     m_cgefwPtr->expression_num_,
//                     gene_count,
//                     exp_count,
//                     dnb_count,
//                     area,
//                     cell_type_id
//             };
//             cnt++;
//             cid++;

//             m_cgefwPtr->cell_attr_.min_x = std::min(m_cgefwPtr->cell_attr_.min_x, cell.x);
//             m_cgefwPtr->cell_attr_.max_x = std::max(m_cgefwPtr->cell_attr_.max_x, cell.x);

//             m_cgefwPtr->cell_attr_.min_y = std::min(m_cgefwPtr->cell_attr_.min_y, cell.y);
//             m_cgefwPtr->cell_attr_.max_y = std::max(m_cgefwPtr->cell_attr_.max_y, cell.y);

//             m_cgefwPtr->cell_attr_.min_area = std::min(m_cgefwPtr->cell_attr_.min_area, area);
//             m_cgefwPtr->cell_attr_.max_area = std::max(m_cgefwPtr->cell_attr_.max_area, area);

//             m_cgefwPtr->cell_attr_.min_gene_count = std::min(m_cgefwPtr->cell_attr_.min_gene_count, gene_count);
//             m_cgefwPtr->cell_attr_.max_gene_count = std::max(m_cgefwPtr->cell_attr_.max_gene_count, gene_count);

//             m_cgefwPtr->cell_attr_.min_exp_count = std::min(m_cgefwPtr->cell_attr_.min_exp_count, exp_count);
//             m_cgefwPtr->cell_attr_.max_exp_count = std::max(m_cgefwPtr->cell_attr_.max_exp_count, exp_count);

//             m_cgefwPtr->cell_attr_.min_dnb_count = std::min(m_cgefwPtr->cell_attr_.min_dnb_count, dnb_count);
//             m_cgefwPtr->cell_attr_.max_dnb_count = std::max(m_cgefwPtr->cell_attr_.max_dnb_count, dnb_count);

//             m_cgefwPtr->expression_num_ += gene_count;
//             m_cgefwPtr->exp_count_sum_ += exp_count;
//             m_cgefwPtr->dnb_count_sum_ += dnb_count;
//             m_cgefwPtr->area_sum_ += area;

//             m_cgefwPtr->cell_list_.emplace_back(std::move(cell));
//         }
//         vec_blkidx.emplace_back(offcnt);
//         offcnt += cnt;
//     }
//     vec_blkidx.emplace_back(cid);
    
//     m_cgefwPtr->cell_num_ = cid;
//     m_cgefwPtr->max_mid_count_ = maxExpmid;
//     int effective_rect[4] ={m_min_x, m_min_y, m_max_x, m_max_y};
//     m_cgefwPtr->storeCellBorderWithAttr(vec_border.data(), cid, effective_rect);

//     m_cgefwPtr->storeCell(m_blocknum, vec_blkidx.data(), m_block_size);
//     m_cgefwPtr->storeCellExp();
//     m_cgefwPtr->storeCellTypeList();
//     m_cgefwPtr->storeCellLabel(vec_cellLabel);
// }

// void cgefCellgem::writeGene()
// {
//     timer st(__FUNCTION__);
//     m_cgefwPtr->gene_num_ = cgefParam::GetInstance()->m_map_gene.size();
//     GeneData *gene_data_list = static_cast<GeneData *>(calloc(m_cgefwPtr->gene_num_ , sizeof(GeneData)));
    
//     unsigned int exp_count, min_exp_count = UINT32_MAX, max_exp_count = 0, offset = 0;
//     unsigned int cell_count, min_cell_count = UINT32_MAX, max_cell_count = 0;
//     unsigned short max_MID_count = 0;
//     vector<GeneExpData> gene_exp_list;
//     gene_exp_list.reserve(m_cgefwPtr->expression_num_);

//     int cid = 0, i = 0;
//     auto itor = cgefParam::GetInstance()->m_map_gene.begin();
//     for(;itor != cgefParam::GetInstance()->m_map_gene.end();itor++,i++)
//     {
//         max_MID_count = 0;
//         cgef_gene *geneptr = itor->second;
//         auto itor_g = geneptr->m_map_geneexp.begin();
//         for(;itor_g != geneptr->m_map_geneexp.end();itor_g++)
//         {
//             cid = m_hash_clabel2cid[itor_g->first];
//             gene_exp_list.emplace_back(cid, itor_g->second);
//             max_MID_count = std::max(max_MID_count, itor_g->second);
//             m_cgefwPtr->max_mid_count_ = std::max(m_cgefwPtr->max_mid_count_, itor_g->second);
//         }

//         cell_count = geneptr->m_map_geneexp.size();
//         gene_data_list[i].cell_count = cell_count;
//         gene_data_list[i].exp_count = geneptr->expcnt;
//         memcpy(gene_data_list[i].gene_name, itor->first.c_str(), itor->first.length());
//         gene_data_list[i].max_mid_count = max_MID_count;
//         gene_data_list[i].offset = offset;
//         offset += cell_count;

//         min_exp_count = std::min(min_exp_count, geneptr->expcnt);
//         max_exp_count = std::max(max_exp_count, geneptr->expcnt);
//         min_cell_count = std::min(min_cell_count, cell_count);
//         max_cell_count = std::max(max_cell_count, cell_count);
//     }

//     m_cgefwPtr->storeGeneAndGeneExp(min_exp_count, max_exp_count, min_cell_count, max_cell_count,
//                                     gene_data_list, gene_exp_list);
//     free(gene_data_list);
// }

// void cgefCellgem::getCelldata()
// {
// //     uint32_t arry[3] = {36708, 36894, 36749};
// //     for(int i=0;i<3;i++)
// //     {
// //         char buf[10]={64};
// //         sprintf(buf, "%d_dnb.txt", i);
// //         FILE *f = fopen(buf, "w");
// //         if(f)
// //         {
// //             cgef_cell *pcell = cgefParam::GetInstance()->m_map_cell[arry[i]];

// //             // vector<Point> out;
// //             // double epsilon = 0.01 * arcLength(pcell->m_vecPoint, true);
// //             // approxPolyDP(pcell->m_vecPoint, out, epsilon, true);
// //             printf("%d %d\n", pcell->m_vecPoint.size(), arry[i]);
// //             for(Point &pt : pcell->m_vecPoint)
// //             {
// //                 fprintf(f, "%d %d\n", pt.x, pt.y);
// //             }
// //             fclose(f);
// //         }
// //     }

// // exit(0);
//     timer st(__FUNCTION__);
    
//     m_rows = cgefParam::GetInstance()->m_max_y - cgefParam::GetInstance()->m_min_y+1;
//     m_cols = cgefParam::GetInstance()->m_max_x - cgefParam::GetInstance()->m_min_x+1;

//     m_block_size[0] = cgefParam::GetInstance()->m_block_size[0];
//     m_block_size[1] = cgefParam::GetInstance()->m_block_size[1];
//     m_block_size[2] = ceil(m_cols * 1.0 / m_block_size[0]); //x_block_num
//     m_block_size[3] = ceil(m_rows * 1.0 / m_block_size[1]); //y_block_num
//     m_blocknum = m_block_size[2] * m_block_size[3];
//     m_vec_veccell.reserve(m_blocknum);
//     for(int i=0;i<m_blocknum;i++)
//     {
//         vector<celldata> vectmp;
//         m_vec_veccell.emplace_back(std::move(vectmp));
//     }

//     bool ret = false;
//     auto itor = cgefParam::GetInstance()->m_map_cell.begin();
//     for(;itor != cgefParam::GetInstance()->m_map_cell.end();itor++)
//     {
//         ret = itor->second->getCenter_border(m_block_size, cgefParam::GetInstance()->m_min_x, cgefParam::GetInstance()->m_min_y);
//         if(!ret) continue;
//         m_vec_veccell[itor->second->m_blkid].emplace_back(0, itor->first); //没有连通域,itor->first是文本分配的cid
//         assert(itor->first == itor->second->m_celllabel);
//         m_maskcellnum++;
//     }
//     printf("fn:%d cn:%d\n", cgefParam::GetInstance()->m_map_cell.size(), m_maskcellnum);
// }

// void cgefCellgem::getCelldata_celltype()
// {
//     timer st(__FUNCTION__);
    
//     m_rows = cgefParam::GetInstance()->m_max_y - cgefParam::GetInstance()->m_min_y +1;
//     m_cols = cgefParam::GetInstance()->m_max_x - cgefParam::GetInstance()->m_min_x +1;
//     m_cgefwPtr->m_x_len = m_cols;
//     m_cgefwPtr->m_y_len = m_rows;

//     m_block_size[0] = cgefParam::GetInstance()->m_block_size[0];
//     m_block_size[1] = cgefParam::GetInstance()->m_block_size[1];
//     m_block_size[2] = ceil(m_cols * 1.0 / m_block_size[0]); //x_block_num
//     m_block_size[3] = ceil(m_rows * 1.0 / m_block_size[1]); //y_block_num
//     m_blocknum = m_block_size[2] * m_block_size[3];
//     m_vec_veccell.reserve(m_blocknum);
//     for(int i=0;i<m_blocknum;i++)
//     {
//         vector<celldata> vectmp;
//         m_vec_veccell.emplace_back(std::move(vectmp));
//     }

//     bool ret = false;
//     int celltypeid = 0;
//     auto itor = cgefParam::GetInstance()->m_map_cell.begin();
//     for(;itor != cgefParam::GetInstance()->m_map_cell.end();itor++)
//     {
//         ret = itor->second->getCenter_median(m_block_size, cgefParam::GetInstance()->m_min_x, cgefParam::GetInstance()->m_min_y);
//         if(!ret) continue;
//         m_vec_veccell[itor->second->m_blkid].emplace_back(0, itor->first); //没有连通域
//         assert(itor->first == itor->second->m_celllabel);
//         if(m_hash_celltype.find(itor->second->m_type) == m_hash_celltype.end())
//         {
//             m_hash_celltype.emplace(itor->second->m_type, celltypeid++);
//             S32 cell_type(itor->second->m_type);
//             m_cgefwPtr->cell_type_list_.emplace_back(std::move(cell_type));
//         }
//         m_maskcellnum++;
//     }
//     printf("fn:%d cn:%d\n", cgefParam::GetInstance()->m_map_cell.size(), m_maskcellnum);
// }

// void cgefCellgem::writeCell_celltype()
// {
//     timer st(__FUNCTION__);
//     unsigned int cid = 0, gid = 0, offcnt = 0; 
//     unsigned short gene_count, exp_count, dnb_count, area, cell_type_id;
//     uint16_t maxExpmid = 0;
//     int cx, cy; //细胞质点
//     vector<unsigned int> vec_blkidx;
//     vec_blkidx.reserve(m_blocknum+1);
//     vector<unsigned int> vec_cellLabel;
//     vec_cellLabel.reserve(m_maskcellnum);

//     vector<short> vec_border;
//     vec_border.reserve(m_maskcellnum*2*BORDERCNT);

//     for(vector<celldata> &vec : m_vec_veccell)
//     {
//         int cnt = 0;
//         for(celldata & cdata : vec)
//         {
//             cgef_cell *cellptr = cgefParam::GetInstance()->m_map_cell[cdata.l_idx];
//             if(cellptr == nullptr) continue;
 
//             cx = cellptr->m_center.x;
//             cy = cellptr->m_center.y;
//             m_hash_clabel2cid.emplace(cdata.l_idx, cid);
//             vec_cellLabel.emplace_back(cdata.l_idx);

//             int bsz = cellptr->m_border.size();
//             int i=0;
//             for(;i<bsz;i++)
//             {
//                 vec_border.emplace_back(cellptr->m_border[i].x - cx);
//                 vec_border.emplace_back(cellptr->m_border[i].y - cy);
//             }

//             for(;i<BORDERCNT;i++) //不足补0
//             {
//                 vec_border.emplace_back(SHRT_MAX);
//                 vec_border.emplace_back(SHRT_MAX);
//             }

//             gene_count = cellptr->m_map_cellexp.size();
//             exp_count = cellptr->expcnt;
//             dnb_count = cellptr->dnbcnt;
//             auto itor = cellptr->m_map_cellexp.begin();
//             for(;itor != cellptr->m_map_cellexp.end();itor++)
//             {
//                 gid = m_hash_gname2gid[itor->first];
//                 m_cgefwPtr->cell_exp_list_.emplace_back(gid, itor->second);
//                 maxExpmid = std::max(maxExpmid, itor->second);
//             }

//             area = cellptr->m_area;
//             if(m_hash_celltype.empty())
//             {
//                 cell_type_id = m_cgefwPtr->random_cell_type_num_ == 0 ? 0 : rand()%(m_cgefwPtr->random_cell_type_num_ + 1);
//             }
//             else
//             {
//                 cell_type_id = m_hash_celltype[cellptr->m_type];
//             }
//             CellData cell = {
//                     cid,
//                     cx,
//                     cy,
//                     m_cgefwPtr->expression_num_,
//                     gene_count,
//                     exp_count,
//                     dnb_count,
//                     area,
//                     cell_type_id
//             };
//             cnt++;
//             cid++;

//             m_cgefwPtr->cell_attr_.min_x = std::min(m_cgefwPtr->cell_attr_.min_x, cell.x);
//             m_cgefwPtr->cell_attr_.max_x = std::max(m_cgefwPtr->cell_attr_.max_x, cell.x);

//             m_cgefwPtr->cell_attr_.min_y = std::min(m_cgefwPtr->cell_attr_.min_y, cell.y);
//             m_cgefwPtr->cell_attr_.max_y = std::max(m_cgefwPtr->cell_attr_.max_y, cell.y);

//             m_cgefwPtr->cell_attr_.min_area = std::min(m_cgefwPtr->cell_attr_.min_area, area);
//             m_cgefwPtr->cell_attr_.max_area = std::max(m_cgefwPtr->cell_attr_.max_area, area);

//             m_cgefwPtr->cell_attr_.min_gene_count = std::min(m_cgefwPtr->cell_attr_.min_gene_count, gene_count);
//             m_cgefwPtr->cell_attr_.max_gene_count = std::max(m_cgefwPtr->cell_attr_.max_gene_count, gene_count);

//             m_cgefwPtr->cell_attr_.min_exp_count = std::min(m_cgefwPtr->cell_attr_.min_exp_count, exp_count);
//             m_cgefwPtr->cell_attr_.max_exp_count = std::max(m_cgefwPtr->cell_attr_.max_exp_count, exp_count);

//             m_cgefwPtr->cell_attr_.min_dnb_count = std::min(m_cgefwPtr->cell_attr_.min_dnb_count, dnb_count);
//             m_cgefwPtr->cell_attr_.max_dnb_count = std::max(m_cgefwPtr->cell_attr_.max_dnb_count, dnb_count);

//             m_cgefwPtr->expression_num_ += gene_count;
//             m_cgefwPtr->exp_count_sum_ += exp_count;
//             m_cgefwPtr->dnb_count_sum_ += dnb_count;
//             m_cgefwPtr->area_sum_ += area;

//             m_cgefwPtr->cell_list_.emplace_back(std::move(cell));
//         }
//         vec_blkidx.emplace_back(offcnt);
//         offcnt += cnt;
//     }
//     vec_blkidx.emplace_back(cid);
    
//     m_cgefwPtr->cell_num_ = cid;
//     m_cgefwPtr->max_mid_count_ = maxExpmid;
//     int effective_rect[4] ={m_min_x, m_min_y, m_max_x, m_max_y};
//     m_cgefwPtr->storeCellBorderWithAttr(vec_border.data(), cid, effective_rect);

//     m_cgefwPtr->storeCell(m_blocknum, vec_blkidx.data(), m_block_size);
//     m_cgefwPtr->storeCellExp();
//     if(cgefParam::GetInstance()->m_intype == INPUTTYPE_GEM_6TYPE)
//     {
//         m_cgefwPtr->storeCellTypeList_N();
//     }
//     else
//     {
//         m_cgefwPtr->storeCellTypeList();
//     }
//     m_cgefwPtr->storeCellLabel(vec_cellLabel);

// }


/////////////////////////////////////////
void cgefCellgem::getCelldata_cgem()
{
    timer st(__FUNCTION__);
    m_rows = cgefParam::GetInstance()->m_max_y - cgefParam::GetInstance()->m_min_y+1;
    m_cols = cgefParam::GetInstance()->m_max_x - cgefParam::GetInstance()->m_min_x+1;

    m_block_size[0] = cgefParam::GetInstance()->m_block_size[0];
    m_block_size[1] = cgefParam::GetInstance()->m_block_size[1];
    m_block_size[2] = ceil(m_cols * 1.0 / m_block_size[0]); //x_block_num
    m_block_size[3] = ceil(m_rows * 1.0 / m_block_size[1]); //y_block_num
    m_blocknum = m_block_size[2] * m_block_size[3];
    m_vec_veccid.reserve(m_blocknum);
    for(int i=0;i<m_blocknum;i++)
    {
        vector<uint32_t> vectmp;
        m_vec_veccid.emplace_back(std::move(vectmp));
    }

    bool ret = false;
    auto itor = cgefParam::GetInstance()->m_map_cell.begin();
    for(;itor != cgefParam::GetInstance()->m_map_cell.end();itor++)
    {
        ret = itor->second->getCenter_border(m_block_size, cgefParam::GetInstance()->m_min_x, cgefParam::GetInstance()->m_min_y);
        if(!ret) continue;
        m_vec_veccid[itor->second->m_blkid].emplace_back(itor->first); //itor->first是文本分配的cid
        assert(itor->first == itor->second->m_celllabel);
        m_maskcellnum++;
    }

    int gid = 0;
    for (auto& p : cgefParam::GetInstance()->m_map_gene_id)//gene name to gid
    {
        p.second = gid++;
    }
    printf("fn:%d cn:%d gn:%d\n", cgefParam::GetInstance()->m_map_cell.size(), m_maskcellnum, cgefParam::GetInstance()->m_map_gene_id.size());
}

void cgefCellgem::writeCell_cgem()
{
    timer st(__FUNCTION__);
    unsigned int cid = 0, gid = 0, offcnt = 0;
    unsigned short gene_count, exp_count, dnb_count, area, cell_type_id;
    int cx, cy; //细胞质点
    vector<unsigned int> vec_blkidx;
    vec_blkidx.reserve(m_blocknum+1);
    vector<unsigned int> vec_cellLabel;
    vec_cellLabel.reserve(m_maskcellnum);

    vector<short> vec_border;
    vec_border.reserve(m_maskcellnum*2*BORDERCNT);

    vector<uint16_t> vec_cellexon;
    vector<uint16_t> vec_cellexon_exp;
    uint16_t minExon = USHRT_MAX, maxExon = 0, maxExpExon = 0, maxExpmid = 0;

    for(vector<uint32_t> &vec : m_vec_veccid)
    {
        int cnt = 0;
        for(uint32_t ocid : vec)
        {
            cgef_cell *cellptr = cgefParam::GetInstance()->m_map_cell[ocid];
            if(cellptr == nullptr) continue;
 
            minExon = std::min(minExon, cellptr->m_exon);
            maxExon = std::max(maxExon, cellptr->m_exon);

            cx = cellptr->m_center.x;
            cy = cellptr->m_center.y;
            m_hash_clabel2cid.emplace(ocid, cid);
            vec_cellLabel.emplace_back(ocid);

            int bsz = cellptr->m_border.size();
            int i=0;
            for(;i<bsz;i++)
            {
                vec_border.emplace_back(cellptr->m_border[i].x - cx);
                vec_border.emplace_back(cellptr->m_border[i].y - cy);
            }

            for(;i<BORDERCNT;i++) //不足补0
            {
                vec_border.emplace_back(SHRT_MAX);
                vec_border.emplace_back(SHRT_MAX);
            }

            gene_count = cellptr->m_map_cellexp.size();
            exp_count = cellptr->expcnt;
            dnb_count = cellptr->dnbcnt;
            vec_cellexon.push_back(cellptr->m_exon);

            auto itor = cellptr->m_map_cellexp.begin();
            for(;itor != cellptr->m_map_cellexp.end();itor++)
            {
                gid = cgefParam::GetInstance()->m_map_gene_id[itor->first];
                m_cgefwPtr->cell_exp_list_.emplace_back(gid, itor->second.umi);
                maxExpmid = std::max(maxExpmid, itor->second.umi);
                vec_cellexon_exp.emplace_back(itor->second.exon);
                maxExpExon = std::max(maxExpExon, itor->second.exon);     

                if(m_map_gene.find(gid) == m_map_gene.end())
                {
                    vector<cellt> tvec;
                    m_map_gene.emplace(gid, std::move(tvec));
                }
                m_map_gene[gid].emplace_back(itor->second.umi, itor->second.exon, cid);
            }

            area = cellptr->m_area;
            cell_type_id = m_cgefwPtr->random_cell_type_num_ == 0 ? 0 : rand()%(m_cgefwPtr->random_cell_type_num_ + 1);
            CellData cell{
                    cid,
                    cx,
                    cy,
                    m_cgefwPtr->expression_num_,
                    gene_count,
                    exp_count,
                    dnb_count,
                    area,
                    cell_type_id
            };
            cnt++;
            cid++;

            m_cgefwPtr->cell_attr_.min_x = std::min(m_cgefwPtr->cell_attr_.min_x, cell.x);
            m_cgefwPtr->cell_attr_.max_x = std::max(m_cgefwPtr->cell_attr_.max_x, cell.x);

            m_cgefwPtr->cell_attr_.min_y = std::min(m_cgefwPtr->cell_attr_.min_y, cell.y);
            m_cgefwPtr->cell_attr_.max_y = std::max(m_cgefwPtr->cell_attr_.max_y, cell.y);

            m_cgefwPtr->cell_attr_.min_area = std::min(m_cgefwPtr->cell_attr_.min_area, area);
            m_cgefwPtr->cell_attr_.max_area = std::max(m_cgefwPtr->cell_attr_.max_area, area);

            m_cgefwPtr->cell_attr_.min_gene_count = std::min(m_cgefwPtr->cell_attr_.min_gene_count, gene_count);
            m_cgefwPtr->cell_attr_.max_gene_count = std::max(m_cgefwPtr->cell_attr_.max_gene_count, gene_count);

            m_cgefwPtr->cell_attr_.min_exp_count = std::min(m_cgefwPtr->cell_attr_.min_exp_count, exp_count);
            m_cgefwPtr->cell_attr_.max_exp_count = std::max(m_cgefwPtr->cell_attr_.max_exp_count, exp_count);

            m_cgefwPtr->cell_attr_.min_dnb_count = std::min(m_cgefwPtr->cell_attr_.min_dnb_count, dnb_count);
            m_cgefwPtr->cell_attr_.max_dnb_count = std::max(m_cgefwPtr->cell_attr_.max_dnb_count, dnb_count);

            m_cgefwPtr->expression_num_ += gene_count;
            m_cgefwPtr->exp_count_sum_ += exp_count;
            m_cgefwPtr->dnb_count_sum_ += dnb_count;
            m_cgefwPtr->area_sum_ += area;

            m_cgefwPtr->cell_list_.emplace_back(std::move(cell));
        }
        vec_blkidx.emplace_back(offcnt);
        offcnt += cnt;
    }
    vec_blkidx.emplace_back(cid);
    
    m_cgefwPtr->cell_num_ = cid;
    m_cgefwPtr->max_mid_count_ = maxExpmid;
    // int effective_rect[4] ={m_min_x, m_min_y, m_max_x, m_max_y};
    int effective_rect[4] = {cgefParam::GetInstance()->m_min_x,cgefParam::GetInstance()->m_min_y,
                             cgefParam::GetInstance()->m_max_x,cgefParam::GetInstance()->m_max_y};
    m_cgefwPtr->storeCellBorderWithAttr(vec_border.data(), cid, effective_rect);

    m_cgefwPtr->storeCell(m_blocknum, vec_blkidx.data(), m_block_size);
    m_cgefwPtr->storeCellExp();
    m_cgefwPtr->storeCellTypeList();
    m_cgefwPtr->storeCellLabel(vec_cellLabel);
    if(m_bexon)
    {
        m_cgefwPtr->storeCellExon(minExon, maxExon, vec_cellexon, maxExpExon, vec_cellexon_exp);
    }
}

void cgefCellgem::writeGene_cgem()
{
    timer st(__FUNCTION__);
    m_cgefwPtr->gene_num_ = cgefParam::GetInstance()->m_map_gene_id.size();
    GeneData *gene_data_list = static_cast<GeneData *>(calloc(m_cgefwPtr->gene_num_ , sizeof(GeneData)));
    
    uint32_t *gene_exon_ptr = (uint32_t *)calloc(m_cgefwPtr->gene_num_, 4);
    vector<uint16_t> vec_exonExp;
    vec_exonExp.reserve(m_geneExpcnt);
    uint16_t maxExpExon = 0;
    uint32_t minExon = UINT_MAX, maxExon = 0;

    unsigned int exp_count, min_exp_count = UINT32_MAX, max_exp_count = 0, offset = 0;
    unsigned int cell_count, min_cell_count = UINT32_MAX, max_cell_count = 0;
    unsigned short max_MID_count = 0;
    vector<GeneExpData> gene_exp_list;
    gene_exp_list.reserve(m_cgefwPtr->expression_num_);

    int i = 0;
    auto itor = cgefParam::GetInstance()->m_map_gene_id.begin();
    for(;itor != cgefParam::GetInstance()->m_map_gene_id.end();itor++,i++)
    {
        max_MID_count = 0;
        vector<cellt> &tvec = m_map_gene[itor->second];
        uint32_t umisum = 0, exonsum = 0;
        for(cellt &ct : tvec)
        {
            gene_exp_list.emplace_back(ct.cid, ct.mid);
            max_MID_count = std::max(max_MID_count, ct.mid);
            m_cgefwPtr->max_mid_count_ = std::max(m_cgefwPtr->max_mid_count_, ct.mid);
            vec_exonExp.emplace_back(ct.exon);
            maxExpExon = std::max(maxExpExon, ct.exon);
            umisum += ct.mid;
            exonsum += ct.exon;
        }
        gene_exon_ptr[i] = exonsum;
        cell_count = tvec.size();
        gene_data_list[i].cell_count = cell_count;
        gene_data_list[i].exp_count = umisum;
        memcpy(gene_data_list[i].gene_name, itor->first.c_str(), itor->first.length());
        gene_data_list[i].max_mid_count = max_MID_count;
        gene_data_list[i].offset = offset;
        offset += cell_count;

        min_exp_count = std::min(min_exp_count, exonsum);
        max_exp_count = std::max(max_exp_count, exonsum);
        min_cell_count = std::min(min_cell_count, cell_count);
        max_cell_count = std::max(max_cell_count, cell_count);
    }

    m_cgefwPtr->storeGeneAndGeneExp(min_exp_count, max_exp_count, min_cell_count, max_cell_count,
                                    gene_data_list, gene_exp_list);

    if(m_bexon)
    {
        m_cgefwPtr->storeGeneExon(minExon, maxExon, gene_exon_ptr, maxExpExon, vec_exonExp);
    }

    free(gene_data_list);
    free(gene_exon_ptr);
}

///////////////////////////////////////////////////////
void cgefCellgem::readBgef_new(const string &strinput)
{
    timer st(__FUNCTION__);
    hid_t file_id = H5Fopen(strinput.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    hsize_t dims[1];
    hid_t gene_did = H5Dopen(file_id, "/geneExp/bin1/gene", H5P_DEFAULT);
    hid_t gene_sid = H5Dget_space(gene_did);
    H5Sget_simple_extent_dims(gene_sid, dims, nullptr);
    
    m_genecnt = dims[0];

    m_genePtr = (Gene*)malloc(dims[0] * sizeof(Gene));
    hid_t genememtype, strtype;
    strtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(strtype, 64);

    genememtype = H5Tcreate(H5T_COMPOUND, sizeof(Gene));
    H5Tinsert(genememtype, "gene", HOFFSET(Gene, gene), strtype);
    H5Tinsert(genememtype, "offset", HOFFSET(Gene, offset), H5T_NATIVE_UINT);
    H5Tinsert(genememtype, "count", HOFFSET(Gene, count), H5T_NATIVE_UINT);
    H5Dread(gene_did, genememtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_genePtr);
    H5Tclose(genememtype);
    H5Sclose(gene_sid);
    H5Dclose(gene_did);


    hid_t exp_did = H5Dopen(file_id, "/geneExp/bin1/expression", H5P_DEFAULT);
    hid_t exp_sid = H5Dget_space(exp_did);
    H5Sget_simple_extent_dims(exp_sid, dims, nullptr);

    m_geneExpcnt = dims[0];

    hid_t memtype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(Expression));
    H5Tinsert(memtype, "x", HOFFSET(Expression, x), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "y", HOFFSET(Expression, y), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "count", HOFFSET(Expression, count), H5T_NATIVE_UINT);

    m_expPtr = (Expression *) calloc(dims[0], sizeof(Expression));
    H5Dread(exp_did, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_expPtr);

    if(H5Lexists(file_id, "/geneExp/bin1/exon", H5P_DEFAULT)>0)
    {
        m_bexon = true;
        hsize_t edims[1];
        hid_t did = H5Dopen(file_id, "/geneExp/bin1/exon", H5P_DEFAULT);
        hid_t sid = H5Dget_space(did);
        H5Sget_simple_extent_dims(sid, edims, nullptr);
        assert(edims[0] == m_geneExpcnt);
        unsigned int *exonPtr = new unsigned int[edims[0]];
        H5Dread(did, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, exonPtr);
        H5Sclose(sid);
        H5Dclose(did);
        for(uint64_t i=0;i<m_geneExpcnt;i++)
        {
            m_expPtr[i].exon = exonPtr[i];
        }
        delete []exonPtr;
    }

    uint64_t l_id = 0 ;
    for(uint32_t i = 0; i < m_genecnt; i++)
    {
        Expression *ptr = m_expPtr + m_genePtr[i].offset;
        for(uint32_t j=0;j<m_genePtr[i].count;j++)
        {
            l_id = ptr[j].x;
            l_id = (l_id << 32) | ptr[j].y;
            if(m_hash_vecdnb.find(l_id) == m_hash_vecdnb.end())
            {
                vector<cellExp_Exon> tvec;
                m_hash_vecdnb.emplace(l_id, tvec);
            }
            m_hash_vecdnb[l_id].emplace_back(i, ptr[j].count, ptr[j].exon);
        }
    }
    free(m_expPtr);
    hid_t attr = H5Aopen(exp_did, "minX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &(cgefParam::GetInstance()->m_min_x));
    attr = H5Aopen(exp_did, "minY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &(cgefParam::GetInstance()->m_min_y));
    attr = H5Aopen(exp_did, "maxX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &(cgefParam::GetInstance()->m_max_x));
    attr = H5Aopen(exp_did, "maxY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &(cgefParam::GetInstance()->m_max_y));
    attr = H5Aopen(exp_did, "resolution", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &(cgefParam::GetInstance()->m_resolution));
    printf("minx:%d miny:%d maxx:%d maxy:%d\n", cgefParam::GetInstance()->m_min_x, 
    cgefParam::GetInstance()->m_min_y, cgefParam::GetInstance()->m_max_x,
    cgefParam::GetInstance()->m_max_y);
    H5Aclose(attr);
    H5Tclose(memtype);
    H5Sclose(exp_sid);
    H5Dclose(exp_did);

    if(H5Aexists(file_id, "omics"))
    {
        hid_t f_attr = H5Aopen(file_id, "omics", H5P_DEFAULT);
        char szbuf[128]={0};
        H5Aread(f_attr, strtype, szbuf);
        m_stromics.clear();
        m_stromics.append(szbuf);
        H5Aclose(f_attr);
    }
    H5Tclose(strtype);

    H5Fclose(file_id);
    printf("genecnt:%d geneExpcnt:%d hashcnt:%d\n", m_genecnt, m_geneExpcnt, m_hash_vecdnb.size());
}

void cgefCellgem::readmask_new(const string &strmask)
{
    timer ct(__FUNCTION__);
    cv::Mat img;
    tifread(img, strmask);
    if (img.empty()) {
        reportErrorCode2File(errorCode::E_INVALIDPARAM, "read mask file error ");
        std::cout << "read mask file error" << std::endl;
        exit(-1);
    }
    assert(!img.empty());
    m_rows = img.rows;
    m_cols = img.cols;
    assert(m_rows == cgefParam::GetInstance()->m_max_y - cgefParam::GetInstance()->m_min_y+1);
    assert(m_cols == cgefParam::GetInstance()->m_max_x - cgefParam::GetInstance()->m_min_x+1);
    if (m_rows != (cgefParam::GetInstance()->m_max_y -
                       cgefParam::GetInstance()->m_min_y + 1)) {
        reportErrorCode2File(errorCode::E_INVALIDPARAM,
                             "mask matrix dismatch gef ");
        std::cout << "mask matrix dismatch gef" << std::endl;
        exit(-1);
    }
    if (m_cols != (cgefParam::GetInstance()->m_max_x -
                      cgefParam::GetInstance()->m_min_x + 1)) {
        reportErrorCode2File(errorCode::E_INVALIDPARAM,
                             "mask matrix dismatch gef ");
        std::cout << "mask matrix dismatch gef" << std::endl;
        exit(-1);
    }

    m_block_size[0] = cgefParam::GetInstance()->m_block_size[0];
    m_block_size[1] = cgefParam::GetInstance()->m_block_size[1];
    m_block_size[2] = ceil(m_cols * 1.0 / m_block_size[0]); //x_block_num
    m_block_size[3] = ceil(m_rows * 1.0 / m_block_size[1]); //y_block_num
    m_blocknum = m_block_size[2] * m_block_size[3];

    vector<cv::Vec4i> hierarchy;
    cv::findContours(img, m_contours, hierarchy, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE);
    m_labelcnt = cv::connectedComponentsWithStats(img, m_outimg, m_stats, m_centroids);
}

void cgefCellgem::getCell()
{
    timer st(__FUNCTION__);
    m_vec_vec_cellunit.reserve(m_blocknum);
    for(int i=0;i<m_blocknum;i++)
    {
        vector<cellUnit*> vectmp;
        m_vec_vec_cellunit.emplace_back(std::move(vectmp));
    }

    int scnt = m_contours.size();
    std::unordered_map<cv::Rect, int, function<size_t (const cv::Rect &)>, function<bool (const cv::Rect &, const cv::Rect &)>> map_cidx(scnt, Rect_hash, Rectequal_to);
    for(int i=0;i<scnt;i++)
    {
        if(m_contours[i].size()>3)
        {
            const cv::Rect &rect = cv::boundingRect(m_contours[i]);
            map_cidx.emplace(rect, i);
        }
    }

    m_cellqueuePtr = new GefQueue<cellUnit>;
    int x, y, w, h, c_idx;
    int cnt = 0;
    for(int i=1;i<m_labelcnt;i++)
    {
        x = m_stats.at<int>(i, cv::CC_STAT_LEFT);
        y = m_stats.at<int>(i, cv::CC_STAT_TOP);
        w = m_stats.at<int>(i, cv::CC_STAT_WIDTH);
        h = m_stats.at<int>(i, cv::CC_STAT_HEIGHT);
        cv::Rect r1(x, y, w, h);

        if(map_cidx.find(r1) != map_cidx.end())
        {
            m_min_x = std::min(m_min_x, x);
            m_max_x = std::max(m_max_x, x+w);
            m_min_y = std::min(m_min_y, y);
            m_max_y = std::max(m_max_y, y+h);

            c_idx = map_cidx.at(r1);
            getcellbinTask *rtask = new getcellbinTask(this, i, r1, m_contours[c_idx]);
            m_thpoolPtr->addTask(rtask);
            cnt++;
        }
    }

    while (cnt--)
    {
        cellUnit *cptr = m_cellqueuePtr->getPtr();
        if(cptr->m_dnbcnt)
        {
            m_vec_vec_cellunit[cptr->m_blkid].emplace_back(cptr);
            m_maskcellnum++;
            m_borcnt += cptr->m_vecborder.size();
        }
        else
        {
            delete cptr;
        }
    }
    printf("borcnt:%d labcnt:%d maskcell %d\n",scnt, m_labelcnt, m_maskcellnum);
}

void cgefCellgem::writeCell_new()
{
    timer st(__FUNCTION__);
    unsigned int offcnt = 0, cid = 0; 
    unsigned short cell_type_id;
    uint32_t offset = 0;
    vector<uint16_t> vec_cellexon;
    vector<uint16_t> vec_cellexon_exp;
    uint16_t minExon = USHRT_MAX, maxExon = 0, maxExpExon = 0, maxExpmid = 0;

    vector<short> vec_border;
    vec_border.reserve(m_borcnt*2);
    m_vec_blkidx.reserve(m_blocknum+1);
    // vector<short> vec_borcnt;
    // vec_borcnt.reserve(m_maskcellnum);

    for(vector<cellUnit*> &vec : m_vec_vec_cellunit)
    {
        for(cellUnit* ptr : vec)
        {
            vec_cellexon.emplace_back(ptr->m_exoncnt);
            minExon = std::min(minExon, ptr->m_exoncnt);
            maxExon = std::max(maxExon, ptr->m_exoncnt);
            
            auto itor = ptr->m_map_gExp.begin();
            for(;itor != ptr->m_map_gExp.end();itor++)
            {
                m_cgefwPtr->cell_exp_list_.emplace_back(itor->first, itor->second.midcnt);
                vec_cellexon_exp.emplace_back(itor->second.exoncnt);
                maxExpmid = std::max(maxExpmid, itor->second.midcnt);
                maxExpExon = std::max(maxExpExon, itor->second.exoncnt);
                if(m_hash_geneunit.find(itor->first) == m_hash_geneunit.end())
                {
                    geneUnit * guptr = new geneUnit();
                    m_hash_geneunit.emplace(itor->first, guptr);
                }
                m_hash_geneunit[itor->first]->add(cid, itor->second.midcnt, itor->second.exoncnt);
            }
//vec_borcnt.push_back(ptr->m_vecborder.size());
            vec_border.insert(vec_border.end(), ptr->m_vecborder.begin(), ptr->m_vecborder.end());
            cell_type_id = m_cgefwPtr->random_cell_type_num_ == 0 ? 0 : rand()%(m_cgefwPtr->random_cell_type_num_ + 1);
            CellData cell {
                    cid++,
                    ptr->m_cx,
                    ptr->m_cy,
                    offset,
                    ptr->m_map_gExp.size(),
                    ptr->m_expcnt,
                    ptr->m_dnbcnt,
                    ptr->m_area,
                    cell_type_id
            };
            delete ptr;

            offset += ptr->m_map_gExp.size();
            m_cgefwPtr->cell_attr_.min_x = std::min(m_cgefwPtr->cell_attr_.min_x, cell.x);
            m_cgefwPtr->cell_attr_.max_x = std::max(m_cgefwPtr->cell_attr_.max_x, cell.x);

            m_cgefwPtr->cell_attr_.min_y = std::min(m_cgefwPtr->cell_attr_.min_y, cell.y);
            m_cgefwPtr->cell_attr_.max_y = std::max(m_cgefwPtr->cell_attr_.max_y, cell.y);

            m_cgefwPtr->cell_attr_.min_area = std::min(m_cgefwPtr->cell_attr_.min_area, cell.area);
            m_cgefwPtr->cell_attr_.max_area = std::max(m_cgefwPtr->cell_attr_.max_area, cell.area);

            m_cgefwPtr->cell_attr_.min_gene_count = std::min(m_cgefwPtr->cell_attr_.min_gene_count, cell.gene_count);
            m_cgefwPtr->cell_attr_.max_gene_count = std::max(m_cgefwPtr->cell_attr_.max_gene_count, cell.gene_count);

            m_cgefwPtr->cell_attr_.min_exp_count = std::min(m_cgefwPtr->cell_attr_.min_exp_count, cell.exp_count);
            m_cgefwPtr->cell_attr_.max_exp_count = std::max(m_cgefwPtr->cell_attr_.max_exp_count, cell.exp_count);

            m_cgefwPtr->cell_attr_.min_dnb_count = std::min(m_cgefwPtr->cell_attr_.min_dnb_count, cell.dnb_count);
            m_cgefwPtr->cell_attr_.max_dnb_count = std::max(m_cgefwPtr->cell_attr_.max_dnb_count, cell.dnb_count);

            m_cgefwPtr->exp_count_sum_ += cell.exp_count;
            m_cgefwPtr->dnb_count_sum_ += cell.dnb_count;
            m_cgefwPtr->area_sum_ += cell.area;

            m_cgefwPtr->cell_list_.emplace_back(std::move(cell));
            
        }
        m_vec_blkidx.emplace_back(offcnt);
        offcnt += vec.size();
    }
    
    m_vec_blkidx.emplace_back(cid);
    printf("cellcnt:%d maxmid:%d\n", cid, maxExpmid);
    m_cgefwPtr->cell_num_ = cid;
    m_cgefwPtr->expression_num_ = offset;
    m_cgefwPtr->max_mid_count_ = maxExpmid;
    int effective_rect[4] ={m_min_x, m_min_y, m_max_x, m_max_y};
    m_cgefwPtr->storeCellBorderWithAttr(vec_border.data(), cid, effective_rect);
//m_cgefwPtr->storeCellBorder_cnt(vec_borcnt);
    m_cgefwPtr->storeCell(m_blocknum, m_vec_blkidx.data(), m_block_size);
    m_cgefwPtr->storeCellExp();
    m_cgefwPtr->storeCellTypeList();

    if(m_bexon)
    {
        m_cgefwPtr->storeCellExon(minExon, maxExon, vec_cellexon, maxExpExon, vec_cellexon_exp);
    }
}

void cgefCellgem::writeGene_new()
{
    timer st(__FUNCTION__);
    GeneData *gene_data_list = static_cast<GeneData *>(calloc(m_genecnt, sizeof(GeneData)));
    vector<GeneExpData> gene_exp_list;
    gene_exp_list.reserve(m_geneExpcnt);

    uint32_t *gene_exon_ptr = (uint32_t *)calloc(m_genecnt, 4);
    vector<uint16_t> vec_exonExp;
    vec_exonExp.reserve(m_geneExpcnt);

    uint16_t maxExpExon = 0;
    uint32_t offset = 0, expsum = 0, exonsum = 0, cellcnt = 0;
    uint32_t minExon = UINT_MAX, maxExon = 0;
    printf("genecnt:%d hashcnt:%d\n", m_genecnt, m_hash_geneunit.size());
    for(uint32_t i=0;i<m_genecnt;i++)
    {
        auto itor = m_hash_geneunit.find(i);
        if(itor != m_hash_geneunit.end()) 
        {
            geneUnit *ptr = itor->second;
            cellcnt = ptr->m_vec_cexp.size();
            expsum = ptr->m_expcnt;
            exonsum = ptr->m_exoncnt;
            memcpy(gene_data_list[i].gene_name, m_genePtr[i].gene, 32);
            gene_data_list[i].cell_count = cellcnt;
            gene_data_list[i].exp_count = expsum;
            gene_data_list[i].max_mid_count = ptr->m_maxmid; 
            gene_data_list[i].offset = offset;
            offset += cellcnt;
            gene_exon_ptr[i] = exonsum;

            for(cexp &ce : ptr->m_vec_cexp)
            {
                gene_exp_list.emplace_back(ce.cellid, ce.midcnt);
                vec_exonExp.emplace_back(ce.exoncnt);
                maxExpExon = std::max(maxExpExon, ce.exoncnt);
            }
            m_cgefwPtr->max_mid_count_ = std::max(m_cgefwPtr->max_mid_count_, ptr->m_maxmid);
            delete ptr;
        }
        else
        {
            memcpy(gene_data_list[i].gene_name, m_genePtr[i].gene, 32);
            gene_data_list[i].cell_count = 0;
            gene_data_list[i].exp_count = 0;
            gene_data_list[i].max_mid_count = 0; 
            gene_data_list[i].offset = 0;
        }

        cgefParam::GetInstance()->m_minExp = std::min(cgefParam::GetInstance()->m_minExp, expsum);
        cgefParam::GetInstance()->m_maxExp = std::max(cgefParam::GetInstance()->m_maxExp, expsum);
        cgefParam::GetInstance()->m_minCell = std::min(cgefParam::GetInstance()->m_minCell, cellcnt);
        cgefParam::GetInstance()->m_maxCell = std::max(cgefParam::GetInstance()->m_maxCell, cellcnt);

        minExon = std::min(minExon, exonsum);
        maxExon = std::max(maxExon, exonsum);
    }

    m_cgefwPtr->gene_num_ = m_genecnt;
    m_cgefwPtr->expression_num_ = gene_exp_list.size();
    m_cgefwPtr->storeGeneAndGeneExp(cgefParam::GetInstance()->m_minExp, cgefParam::GetInstance()->m_maxExp,
                                    cgefParam::GetInstance()->m_minCell, cgefParam::GetInstance()->m_maxCell,
                                    gene_data_list, gene_exp_list);

    if(m_bexon)
    {
        m_cgefwPtr->storeGeneExon(minExon, maxExon, gene_exon_ptr, maxExpExon, vec_exonExp);
    }
    free(gene_data_list);
    free(gene_exon_ptr);
}
