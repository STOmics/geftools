#include "cgef3d.h"

cgef3d::cgef3d(/* args */)
{
    m_thpoolPtr = new ThreadPool(cgef3dParam::GetInstance()->m_threadcnt);
}

cgef3d::~cgef3d()
{
    delete m_thpoolPtr;
}

void cgef3d::writeCgef(const string &strin, const string &strtxt, const string &strmask, const string &strout)
{
    hid_t fileid = H5Fcreate(strout.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hid_t tgid = H5Gcreate(fileid, "/cellBin", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(tgid);
    m_gid_3d = H5Gcreate(fileid, "/3D", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(m_gid_3d > 0);
    int ret = gemAnalysis(strin);
    readgem_5();
    readtxt(strtxt);
    readmask(strmask);
    storeGene();
    storeCell();
    storeAttr(fileid);
    H5Gclose(m_gid_3d);
    H5Fclose(fileid);
}

int cgef3d::gemAnalysis(const string &strinput)
{
    cgef3dParam::GetInstance()->m_infile = gzopen(strinput.c_str(), "r");
    gzbuffer(cgef3dParam::GetInstance()->m_infile, READLEN);

    char buf[128] = {0};
    while (1)
    {
        gzgets(cgef3dParam::GetInstance()->m_infile, buf, 128);
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
    return col;
}

void cgef3d::readgem_5()
{
    for(int i=0;i<cgef3dParam::GetInstance()->m_threadcnt;i++)
    {
        readFloatTask *rtask = new readFloatTask();
        m_thpoolPtr->addTask(rtask);
    }
    m_thpoolPtr->waitTaskDone();
    gzclose(cgef3dParam::GetInstance()->m_infile);

    printf("genecnt:%ld cellcnt:%ld \n", cgef3dParam::GetInstance()->m_map_gene.size(),
                                    cgef3dParam::GetInstance()->m_map_cell.size());
}

void cgef3d::readtxt(const string &strtxt)
{
    if(strtxt.empty()) return;
    string strctype;
    uint16_t ctid = 0;
    map<string, uint16_t> tmap;
    gzFile fin = gzopen(strtxt.c_str(), "r");
    if(fin)
    {
        char buf[128]={0};
        gzgets(fin, buf, 128);
        uint32_t cid = 0;
        string strtmp;
        while (gzgets(fin, buf, 128)!=NULL)
        {
            int len = strlen(buf);
            buf[len-1] = '\0';
            char *p = strtok(buf, "\t");
            cid = atoi(p);
            p = strtok(NULL, "\t");
            strtmp.clear();
            strtmp.append(p);
            if(tmap.find(strtmp) == tmap.end())
            {
                tmap.emplace(strtmp, ctid++);
                strctype.append(strtmp).append("\t");
            }
            
            m_hash_cell2ctype.emplace(cid, tmap[strtmp]);
        }
        gzclose(fin);

        hsize_t dims[1] = {strctype.size()};
        hid_t did = h5DatasetWrite(m_gid_3d, H5T_STD_U8LE, H5T_NATIVE_CHAR, "cellTypeList", 1, dims, strctype.c_str());
        H5Dclose(did);
    }
}


// void cgef3d::readgem_6()
// {
//     for(int i=0;i<cgef3dParam::GetInstance()->m_threadcnt;i++)
//     {
//         readFloatTask_6 *rtask = new readFloatTask_6();
//         m_thpoolPtr->addTask(rtask);
//     }
//     m_thpoolPtr->waitTaskDone();
//     gzclose(cgef3dParam::GetInstance()->m_infile);

//     string strctype;
//     uint16_t ctid = 0;
//     auto itor = cgef3dParam::GetInstance()->m_map_ctype.begin();
//     for(;itor != cgef3dParam::GetInstance()->m_map_ctype.end();itor++)
//     {
//         for(uint32_t cid : itor->second->m_vec_cid)
//         {
//             m_hash_cell2ctype.emplace(cid, ctid);
//         }
//         strctype.append(itor->first).append("\t");
//         ctid++;
//     }

//     hsize_t dims[1] = {strctype.size()};
//     hid_t did = h5DatasetWrite(m_gid_3d, H5T_STD_U8LE, H5T_NATIVE_CHAR, "cellTypeList", 1, dims, strctype.c_str());
//     H5Dclose(did);
// }

void cgef3d::readmask(const string &strmask)
{
    if(strmask.empty()) return;
    cv::Mat img;
    tifread(img, strmask);
    if (img.empty())
    {
        reportErrorCode2File(errorCode::E_INVALIDPARAM, "read mask file error ");
        std::cout << "read mask file error" << std::endl;
        exit(-1);
    }
    assert(!img.empty());

    vector<vector<cv::Point>> contours;
    vector<cv::Vec4i> hierarchy;
    findContours(img, contours, hierarchy, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE);
    int scnt = contours.size();
    std::unordered_map<cv::Rect, int, function<size_t (const cv::Rect &)>, function<bool (const cv::Rect &, const cv::Rect &)>> map_cidx(scnt, Rect_hash, Rectequal_to);
    
    for(int i=0;i<scnt;i++)
    {
        if(contours[i].size()>3)
        {
            const cv::Rect &rect = cv::boundingRect(contours[i]);
            map_cidx.emplace(rect, i);
        }
    }
    cv::Mat outimg, stats, centroids;
    int count = connectedComponentsWithStats(img, outimg, stats, centroids);
    int  x, y, w, h, c_idx;
    uint32_t maskcellnum = 0;
    for(int i=1;i<count;i++)
    {
        x = stats.at<int>(i, cv::CC_STAT_LEFT);
        y = stats.at<int>(i, cv::CC_STAT_TOP);
        w = stats.at<int>(i, cv::CC_STAT_WIDTH);
        h = stats.at<int>(i, cv::CC_STAT_HEIGHT);
        cv::Rect r1(x, y, w, h);

        if(map_cidx.find(r1) != map_cidx.end())
        {
            c_idx = map_cidx.at(r1);
            cgef3d_cell *ptr = cgef3dParam::GetInstance()->m_map_cell[i];
            if(ptr != nullptr)
            {
                ptr->setCellInfo(centroids.at<double>(i,0), centroids.at<double>(i,1),
                stats.at<int>(i, cv::CC_STAT_AREA), contours[c_idx]);
            }
            maskcellnum++;
        }
    }
    printf("mask cellnum %d\n", maskcellnum);
}

void cgef3d::addCellborder(vector<float> &vec_border, vector<cv::Point2f> &vborder)
{
    int i=0;
    int sz = vborder.size();
    if(sz > BORDERCNT)
    {
        vector<cv::Point2f> tmpborder;
        double epsilon = 0.01 * cv::arcLength(vborder, true);
        cv::approxPolyDP(vborder, tmpborder, epsilon, true);

        sz = tmpborder.size();
        for(;i<sz;i++)
        {
            vec_border.emplace_back(tmpborder[i].x);
            vec_border.emplace_back(tmpborder[i].y);
        }
    }
    else
    {
        for(;i<sz;i++)
        {
            vec_border.emplace_back(vborder[i].x);
            vec_border.emplace_back(vborder[i].y);
        }
    }
    
    for(;i<BORDERCNT;i++)
    {
        vec_border.emplace_back(FLT_MAX);
        vec_border.emplace_back(FLT_MAX);
    }
}

void cgef3d::storeCell()
{
    uint32_t cid = 0;
    vector<float> vecborder;
    vector<cell_3d> veccell;
    map<uint32_t, std::vector<geneexp_3d>> map_gene;
    float minx = FLT_MAX, miny = FLT_MAX, maxx = FLT_MIN, maxy = FLT_MIN, maxumi = 0.0;
    auto itor = cgef3dParam::GetInstance()->m_map_cell.begin();
    for(;itor!=cgef3dParam::GetInstance()->m_map_cell.end();itor++)
    {
        uint16_t ctid = m_hash_cell2ctype.empty() ? 0:m_hash_cell2ctype[itor->first];
        uint16_t gcnt = m_hash_cell2gene[itor->first].size();
        cgef3d_cell *ptr = itor->second;
        bool ret = ptr->getCellInfo();
        if(!ret) continue;

        
        for(cellexp_3d &ce : m_hash_cell2gene[itor->first])
        {
            if(map_gene.find(ce.gid) == map_gene.end())
            {
                std::vector<geneexp_3d> tvec;
                map_gene.emplace(ce.gid, std::move(tvec));
            }
            map_gene[ce.gid].emplace_back(cid, ce.umi);
        }

        addCellborder(vecborder, ptr->m_border);
        minx = std::min(minx, ptr->m_x);
        miny = std::min(miny, ptr->m_y);
        maxx = std::max(maxx, ptr->m_x);
        maxy = std::max(maxy, ptr->m_y);
        maxumi = std::max(maxumi, ptr->m_sumumi);

        veccell.emplace_back(ctid, ptr->m_area, gcnt, ptr->m_dnbcnt, cid++, ptr->m_x, ptr->m_y, ptr->m_sumumi);
        delete itor->second;
    }

    hid_t filetype = H5Tcreate(H5T_COMPOUND, sizeof(cell_3d));
    H5Tinsert(filetype, "ctypeid", HOFFSET(cell_3d, ctypeid), H5T_STD_U16LE);
    H5Tinsert(filetype, "area", HOFFSET(cell_3d, area), H5T_STD_U16LE);
    H5Tinsert(filetype, "genecnt", HOFFSET(cell_3d, genecnt), H5T_STD_U16LE);
    H5Tinsert(filetype, "dnbcnt", HOFFSET(cell_3d, dnbcnt), H5T_STD_U16LE);
    H5Tinsert(filetype, "id", HOFFSET(cell_3d, id), H5T_STD_U32LE);
    H5Tinsert(filetype, "x", HOFFSET(cell_3d, x), H5T_IEEE_F32LE);
    H5Tinsert(filetype, "y", HOFFSET(cell_3d, y), H5T_IEEE_F32LE);
    H5Tinsert(filetype, "sumumi", HOFFSET(cell_3d, sumumi), H5T_IEEE_F32LE);

    hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(cell_3d));
    H5Tinsert(memtype, "ctypeid", HOFFSET(cell_3d, ctypeid), H5T_NATIVE_USHORT);
    H5Tinsert(memtype, "area", HOFFSET(cell_3d, area), H5T_NATIVE_USHORT);
    H5Tinsert(memtype, "genecnt", HOFFSET(cell_3d, genecnt), H5T_NATIVE_USHORT);
    H5Tinsert(memtype, "dnbcnt", HOFFSET(cell_3d, dnbcnt), H5T_NATIVE_USHORT);
    H5Tinsert(memtype, "id", HOFFSET(cell_3d, id), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "x", HOFFSET(cell_3d, x), H5T_NATIVE_FLOAT);
    H5Tinsert(memtype, "y", HOFFSET(cell_3d, y), H5T_NATIVE_FLOAT);
    H5Tinsert(memtype, "sumumi", HOFFSET(cell_3d, sumumi), H5T_NATIVE_FLOAT);

    hsize_t dims[1] = {veccell.size()};
    hid_t did = h5DatasetWrite(m_gid_3d, filetype, memtype, "cell", 1, dims, veccell.data());
    dims[0] = 1;
    h5AttrWrite(did, H5T_IEEE_F32LE, H5T_NATIVE_FLOAT, "minX", 1, dims, &minx);
    h5AttrWrite(did, H5T_IEEE_F32LE, H5T_NATIVE_FLOAT, "minY", 1, dims, &miny);
    h5AttrWrite(did, H5T_IEEE_F32LE, H5T_NATIVE_FLOAT, "maxX", 1, dims, &maxx);
    h5AttrWrite(did, H5T_IEEE_F32LE, H5T_NATIVE_FLOAT, "maxY", 1, dims, &maxy);
    h5AttrWrite(did, H5T_IEEE_F32LE, H5T_NATIVE_FLOAT, "maxumi", 1, dims, &maxumi);

    H5Tclose(filetype);
    H5Tclose(memtype);
    H5Dclose(did);

    dims[0] = vecborder.size();
    hid_t did_b = h5DatasetWrite(m_gid_3d, H5T_IEEE_F32LE, H5T_NATIVE_FLOAT, "cellBorder", 1, dims, vecborder.data());
    H5Dclose(did_b);

    vector<geneexp_3d> vecgeneexp;
    auto itor_gexp = map_gene.begin();
    for(;itor_gexp!=map_gene.end();itor_gexp++)
    {
        vecgeneexp.insert(vecgeneexp.end(), itor_gexp->second.begin(), itor_gexp->second.end());
    }

    hid_t f_gexp = H5Tcreate(H5T_COMPOUND, sizeof(geneexp_3d));
    H5Tinsert(f_gexp, "cid", HOFFSET(geneexp_3d, cid), H5T_STD_U32LE);
    H5Tinsert(f_gexp, "umi", HOFFSET(geneexp_3d, umi), H5T_IEEE_F32LE);
    hid_t m_gexp = H5Tcreate(H5T_COMPOUND, sizeof(geneexp_3d));
    H5Tinsert(m_gexp, "cid", HOFFSET(geneexp_3d, cid), H5T_NATIVE_UINT);
    H5Tinsert(m_gexp, "umi", HOFFSET(geneexp_3d, umi), H5T_NATIVE_FLOAT);

    dims[0] = vecgeneexp.size();
    hid_t didexp = h5DatasetWrite(m_gid_3d, f_gexp, m_gexp, "geneExp", 1, dims, vecgeneexp.data());
    H5Tclose(f_gexp);
    H5Tclose(m_gexp);
    H5Dclose(didexp);
}


void cgef3d::storeGene()
{
    uint32_t gid = 0, offset = 0;
    vector<gene_3d> vecgene;
    auto itor = cgef3dParam::GetInstance()->m_map_gene.begin();
    for(;itor!=cgef3dParam::GetInstance()->m_map_gene.end();itor++)
    {
        float maxumi = 0.0;
        auto &tmap = itor->second->m_map_cell;
        auto itor_t = tmap.begin();
        for(;itor_t != tmap.end();itor_t++)
        {
            maxumi = std::max(maxumi, itor_t->second);
            if(m_hash_cell2gene.find(itor_t->first) == m_hash_cell2gene.end())
            {
                std::vector<cellexp_3d> tmp;
                m_hash_cell2gene.emplace(itor_t->first, std::move(tmp));
            }
            m_hash_cell2gene[itor_t->first].emplace_back(gid, itor_t->second);
        }
        vecgene.emplace_back(offset, tmap.size(), itor->second->m_sumumi, maxumi, itor->first.c_str());
        offset += tmap.size();
        gid++;
        delete itor->second;
    }

    hid_t str32 = H5Tcopy(H5T_C_S1);
    H5Tset_size(str32, 32);

    hid_t f_gene = H5Tcreate(H5T_COMPOUND, sizeof(gene_3d));
    H5Tinsert(f_gene, "offset", HOFFSET(gene_3d, offset), H5T_STD_U32LE);
    H5Tinsert(f_gene, "cellcnt", HOFFSET(gene_3d, cellcnt), H5T_STD_U32LE);
    H5Tinsert(f_gene, "sumumi", HOFFSET(gene_3d, sumumi), H5T_IEEE_F32LE);
    H5Tinsert(f_gene, "maxumi", HOFFSET(gene_3d, maxumi), H5T_IEEE_F32LE);
    H5Tinsert(f_gene, "gene", HOFFSET(gene_3d, gene), str32);

    hid_t m_gene = H5Tcreate(H5T_COMPOUND, sizeof(gene_3d));
    H5Tinsert(m_gene, "offset", HOFFSET(gene_3d, offset), H5T_NATIVE_UINT);
    H5Tinsert(m_gene, "cellcnt", HOFFSET(gene_3d, cellcnt), H5T_NATIVE_UINT);
    H5Tinsert(m_gene, "sumumi", HOFFSET(gene_3d, sumumi), H5T_NATIVE_FLOAT);
    H5Tinsert(m_gene, "maxumi", HOFFSET(gene_3d, maxumi), H5T_NATIVE_FLOAT);
    H5Tinsert(m_gene, "gene", HOFFSET(gene_3d, gene), str32);

    hsize_t dims[1] = {vecgene.size()};
    hid_t did = h5DatasetWrite(m_gid_3d, f_gene, m_gene, "gene", 1, dims, vecgene.data());
    H5Tclose(str32);
    H5Tclose(f_gene);
    H5Tclose(m_gene);
    H5Dclose(did);
}

void cgef3d::storeAttr(hid_t fileid)
{
    hsize_t dimsAttr[1]={1};
    uint32_t version = 2;
    uint32_t resolution = 0;
    int offsetX = 0;
    int offsetY = 0;
    h5AttrWrite(fileid, H5T_STD_U32LE, H5T_NATIVE_UINT32, "version", 1, dimsAttr, &version);
    h5AttrWrite(fileid, H5T_STD_U32LE, H5T_NATIVE_UINT32, "resolution", 1, dimsAttr, &resolution);
    h5AttrWrite(fileid, H5T_STD_I32LE, H5T_NATIVE_INT32, "offsetX", 1, dimsAttr, &offsetX);
    h5AttrWrite(fileid, H5T_STD_I32LE, H5T_NATIVE_INT32, "offsetY", 1, dimsAttr, &offsetY);
    dimsAttr[0] = 3;
    h5AttrWrite(fileid, H5T_STD_U32LE, H5T_NATIVE_UINT32, "geftool_ver", 1, dimsAttr, GEFVERSION);
    string stromics(" ");
    dimsAttr[0] = stromics.size();
    hid_t str32 = H5Tcopy(H5T_C_S1);
    H5Tset_size(str32, 32);
    h5AttrWrite(fileid, str32, str32, "omics", 1, dimsAttr, stromics.c_str());
    H5Tclose(str32);
}