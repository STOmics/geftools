#include "geftogem.h"
#include <sstream>
#include <fstream>
#include "utils.h"

#define GEM_HEADER "#FileFormat=GEMv%d.%d\n" \
"#SortedBy=None\n" \
"#BinSize=%d\n" \
"#STOmicsChip=%s\n" \
"#OffsetX=%d\n" \
"#OffsetY=%d\n" 


geftogem::geftogem(const string &strout, const string &strsn, bool outexon):
    m_strout(strout), m_strsn(strsn), m_boutexon(outexon)
{
}

geftogem::~geftogem()
{
}

void geftogem::getBgefGene(hid_t file_id)
{
    char buf[128]={0};
    sprintf(buf, "/geneExp/bin%d/gene", m_bin);
    hsize_t dims[1];
    hid_t gene_did = H5Dopen(file_id, buf, H5P_DEFAULT);
    hid_t gene_sid = H5Dget_space(gene_did);
    H5Sget_simple_extent_dims(gene_sid, dims, nullptr);

    m_genencnt = dims[0];
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
    H5Tclose(strtype);
    H5Sclose(gene_sid);
    H5Dclose(gene_did);
}

void geftogem::getBgefExp(hid_t file_id)
{
    char buf[128]={0};
    sprintf(buf, "/geneExp/bin%d/expression", m_bin);
    hsize_t dims[1];
    hid_t exp_did = H5Dopen(file_id, buf, H5P_DEFAULT);
    hid_t exp_sid = H5Dget_space(exp_did);
    H5Sget_simple_extent_dims(exp_sid, dims, nullptr);

    m_geneexpcnt = dims[0];

    hid_t memtype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(Expression));
    H5Tinsert(memtype, "x", HOFFSET(Expression, x), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "y", HOFFSET(Expression, y), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "count", HOFFSET(Expression, count), H5T_NATIVE_UINT);

    m_expPtr = (Expression *) malloc(dims[0] * sizeof(Expression));
    H5Dread(exp_did, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_expPtr);

    int l = sprintf(buf, "/geneExp/bin%d/exon", m_bin);
    buf[l]='\0';
    if(H5Lexists(file_id, buf, H5P_DEFAULT)>0)
    {
        m_bexon = true;
        hsize_t edims[1];
        hid_t did = H5Dopen(file_id, buf, H5P_DEFAULT);
        hid_t sid = H5Dget_space(did);
        H5Sget_simple_extent_dims(sid, edims, nullptr);
        assert(edims[0] == m_geneexpcnt);
        unsigned int *exonPtr = new unsigned int[edims[0]];
        H5Dread(did, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, exonPtr);
        H5Sclose(sid);
        H5Dclose(did);
        for(uint64_t i=0;i<m_geneexpcnt;i++)
        {
            m_expPtr[i].exon = exonPtr[i];
        }
        delete []exonPtr;
    }

    hid_t attr = H5Aopen(exp_did, "minX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_min_x);
    attr = H5Aopen(exp_did, "minY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_min_y);
    attr = H5Aopen(exp_did, "maxX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_max_x);
    attr = H5Aopen(exp_did, "maxY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_max_y);
    attr = H5Aopen(exp_did, "resolution", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &m_resolution);
    //printf("minx:%d miny:%d maxx:%d maxy:%d\n", m_min_x, m_min_y, m_max_x, m_max_y);
    H5Aclose(attr);
    H5Tclose(memtype);
    H5Sclose(exp_sid);
    H5Dclose(exp_did);
}

void geftogem::bgef2gem()
{
    ostream *pout = nullptr;
    if(m_strout == "stdout")
    {
        pout = &std::cout;
    }
    else
    {
        pout = new fstream(m_strout.c_str(), ios_base::out);
    }

    stringstream sstrout;
    char buf[1024] = {0};
    sprintf(buf, GEM_HEADER, 0, 1, m_bin, m_strsn.c_str(), m_min_x, m_min_y);

    if(m_bexon && m_boutexon)
    {
        sstrout<<buf<<"geneID\tx\ty\tMIDCount\tExonCount\n";
        *pout<<sstrout.str();
        for (uint32_t i = 0; i < m_genencnt; ++i)
        {
            sstrout.clear();
            sstrout.str("");
            Expression *ptr = m_expPtr + m_genePtr[i].offset;
            for(uint32_t j=0;j<m_genePtr[i].count;j++)
            {
                sstrout<<m_genePtr[i].gene<<'\t'<<ptr[j].x<<'\t'<<ptr[j].y<<'\t'<<ptr[j].count<<'\t'<<ptr[j].exon<<'\n';
            }
            *pout<<sstrout.str();
        }
    }
    else
    {
        sstrout<<buf<<"geneID\tx\ty\tMIDCount\n";
        *pout<<sstrout.str();
        for (uint32_t i = 0; i < m_genencnt; ++i)
        {
            sstrout.clear();
            sstrout.str("");
            Expression *ptr = m_expPtr + m_genePtr[i].offset;
            for(uint32_t j=0;j<m_genePtr[i].count;j++)
            {
                sstrout<<m_genePtr[i].gene<<'\t'<<ptr[j].x<<'\t'<<ptr[j].y<<'\t'<<ptr[j].count<<'\n';
            }
            *pout<<sstrout.str();
        }
    }
    pout->flush();
    if(m_strout != "stdout")
    {
        delete pout;
    }
    free(m_genePtr);
    free(m_expPtr);
}

void geftogem::getdnb()
{
    uint64_t l_id = 0;
    if(m_bexon)
    {
        for(uint32_t i=0;i<m_genencnt;i++)
        {
            m_vecgenename.emplace_back(m_genePtr[i].gene);
            Expression *ptr = m_expPtr + m_genePtr[i].offset;
            for (uint32_t j = 0; j < m_genePtr[i].count; j++) {
                l_id = ptr[j].x;
                l_id = (l_id<<32) | ptr[j].y;
                
                if(m_hash_vecdnb_exon.find(l_id) == m_hash_vecdnb_exon.end())
                {
                    vector<Dnbs_exon> tvec;
                    m_hash_vecdnb_exon.emplace(l_id, tvec);
                }
                m_hash_vecdnb_exon[l_id].emplace_back(i, ptr[j].count, ptr[j].exon);
            }
        }
        printf("gene:%d geneexp:%d hashcnt:%d\n", m_genencnt, m_geneexpcnt, m_hash_vecdnb_exon.size());
    }
    else
    {
        for(uint32_t i=0;i<m_genencnt;i++)
        {
            m_vecgenename.emplace_back(m_genePtr[i].gene);
            Expression *ptr = m_expPtr + m_genePtr[i].offset;
            for(uint32_t j=0;j<m_genePtr[i].count;j++)
            {
                l_id = ptr[j].x;
                l_id = (l_id<<32) | ptr[j].y;
                
                if(m_hash_vecdnb.find(l_id) == m_hash_vecdnb.end())
                {
                    vector<Dnbs> tvec;
                    m_hash_vecdnb.emplace(l_id, tvec);
                }
                m_hash_vecdnb[l_id].emplace_back(i, ptr[j].count);
            }
        }
        printf("gene:%d geneexp:%d hashcnt:%d\n", m_genencnt, m_geneexpcnt, m_hash_vecdnb.size());
    }
    free(m_genePtr);
    free(m_expPtr);
}

void geftogem::readBgef(const string &strinput)
{
    hid_t file_id = H5Fopen(strinput.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    getBgefGene(file_id);
    getBgefExp(file_id);
    H5Fclose(file_id);
}

void geftogem::readCgef(const string &strinput)
{
    hid_t file_id = H5Fopen(strinput.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t cell_did = H5Dopen(file_id, "/cellBin/cell", H5P_DEFAULT);

    hsize_t dims[1];
    hid_t cell_sid = H5Dget_space(cell_did);
    H5Sget_simple_extent_dims(cell_sid, dims, nullptr);

    m_cellcnt = dims[0];
    hid_t memtype = getMemtypeOfCellData();
    CellData *cell_arrayptr = new CellData[dims[0]];
    H5Dread(cell_did, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell_arrayptr);

    H5Tclose(memtype);
    H5Sclose(cell_sid);
    H5Dclose(cell_did);

    hsize_t bdims[3];
    hid_t border_id = H5Dopen(file_id, "/cellBin/cellBorder", H5P_DEFAULT);
    hid_t border_sid = H5Dget_space(border_id);
    H5Sget_simple_extent_dims(border_sid, bdims, nullptr);

    short *borderdataPtr = (short*)calloc(bdims[0]*bdims[1]*bdims[2], 2);
    H5Dread(border_id, H5T_STD_I16LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, borderdataPtr);

    int x,y;
    vector<cv::Point> vecborder;
    vector<cv::Point> rvecborder;
    short *ptmp = borderdataPtr;
    for(uint32_t i=0;i<bdims[0];i++)
    {
        vecborder.clear();
        for(uint32_t j=0;j<bdims[1];j++)
        {
            x = ptmp[2*j];
            y = ptmp[2*j+1];
            if(x==SHRT_MAX && y==SHRT_MAX)
            {
                break;
            }
            x += cell_arrayptr[i].x;
            y += cell_arrayptr[i].y;
            vecborder.emplace_back(x,y);
        }

        if(!vecborder.empty())
        {
            rvecborder.clear();
            cv::Rect rtmp = boundingRect(vecborder);
            cv::Mat t = cv::Mat::zeros(rtmp.height, rtmp.width, CV_8UC1);
            for(cv::Point &pt : vecborder)
            {
                rvecborder.emplace_back(pt.x - rtmp.x, pt.y - rtmp.y);
            }
            fillPoly(t, rvecborder, 1);
            cellmat cm;
            cm.offsetx = rtmp.x;
            cm.offsety = rtmp.y;
            findNonZero(t, cm.vecP);
            m_hash_cellpoint.emplace(i, std::move(cm));
        }

        ptmp += BORDERCNT*2;
    }

    delete cell_arrayptr;
    free(borderdataPtr);
    //printf("cellcnt:%d hashrect %d\n", m_cellcnt, m_hash_cellpoint.size());


    int min_x, min_y, max_x, max_y;
    hid_t attr = H5Aopen(border_id, "minX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &min_x);
    attr = H5Aopen(border_id, "minY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &min_y);
    attr = H5Aopen(border_id, "maxX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &max_x);
    attr = H5Aopen(border_id, "maxY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &max_y);
    //printf("minx:%d miny:%d maxx:%d maxy:%d\n", min_x, min_y, max_x, max_y);

    attr = H5Aopen(file_id, "offsetX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT32, &m_offsetX);

    attr = H5Aopen(file_id, "offsetY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT32, &m_offsetY);
    //printf("offsetx:%d offsety:%d\n", m_offsetX, m_offsetY);
    H5Aclose(attr);
    H5Sclose(border_sid);
    H5Dclose(border_id);
    H5Fclose(file_id);
}

void geftogem::bgeftogem(const string &strbgef, int binsize)
{
    m_bin = binsize;
    readBgef(strbgef);
    bgef2gem();
}

void geftogem::cgeftogem(const string &strcgef, const string &strbgef)
{
    m_bin = 1;
    readBgef(strbgef);
    getdnb();
    readCgef(strcgef);

    if(m_bexon && m_boutexon)
    {
        cgef2gem_exon();
    }
    else
    {
        cgef2gem();
    }
}

void geftogem::cgef2gem()
{
    ostream *pout = nullptr;
    if(m_strout == "stdout")
    {
        pout = &std::cout;
    }
    else
    {
        pout = new fstream(m_strout.c_str(), ios_base::out);
    }
    stringstream sstrout;
    char buf[1024] = {0};
    sprintf(buf, GEM_HEADER, 0, 1, m_bin, m_strsn.c_str(), m_min_x, m_min_y);
    sstrout<<buf<<"geneID\tx\ty\tMIDCount\tCellID\n";
    *pout<<sstrout.str();
    uint64_t l_id = 0;
    int x,y;
    auto itor = m_hash_cellpoint.begin();
    for(;itor != m_hash_cellpoint.end();itor++)
    {
        vector<cv::Point> &vecpoint = itor->second.vecP;
        sstrout.clear();
        sstrout.str("");
        for(cv::Point &pt : vecpoint)
        {
            x = pt.x+itor->second.offsetx;
            y = pt.y+itor->second.offsety;
            l_id = x;
            l_id = (l_id << 32) | y;
            auto dnb_itor = m_hash_vecdnb.find(l_id);
            if(dnb_itor!= m_hash_vecdnb.end())
            {
                for(Dnbs &dnbs : dnb_itor->second)
                {
                    sstrout<<m_vecgenename[dnbs.geneid]<<'\t'<<x<<'\t'<<y<<'\t'<<dnbs.midcnt<<'\t'<<itor->first<<'\n';
                }
                m_hash_vecdnb.erase(l_id);
            }
        }
        *pout<<sstrout.str();
    }
    pout->flush();
    if(m_strout != "stdout")
    {
        delete pout;
    }
}

void geftogem::cgef2gem_exon()
{
    ostream *pout = nullptr;
    if(m_strout == "stdout")
    {
        pout = &std::cout;
    }
    else
    {
        pout = new fstream(m_strout.c_str(), ios_base::out);
    }

    stringstream sstrout;
    char buf[1024] = {0};
    sprintf(buf, GEM_HEADER, 0, 1, m_bin, m_strsn.c_str(), m_min_x, m_min_y);
    sstrout<<buf<<"geneID\tx\ty\tMIDCount\tExonCount\tCellID\n";
    *pout<<sstrout.str();

    uint64_t l_id = 0;
    int x,y;
    auto itor = m_hash_cellpoint.begin();
    for(;itor != m_hash_cellpoint.end();itor++)
    {
        vector<cv::Point> &vecpoint = itor->second.vecP;
        sstrout.clear();
        sstrout.str("");
        for(cv::Point &pt : vecpoint)
        {
            x = pt.x+itor->second.offsetx;
            y = pt.y+itor->second.offsety;
            l_id = x;
            l_id = (l_id << 32) | y;
            auto dnb_itor = m_hash_vecdnb_exon.find(l_id);
            if(dnb_itor!= m_hash_vecdnb_exon.end())
            {
                for(Dnbs_exon &dnbs : dnb_itor->second)
                {
                    sstrout<<m_vecgenename[dnbs.geneid]<<'\t'<<x<<'\t'<<y<<'\t'<<dnbs.midcnt<<'\t'<<dnbs.exon<<'\t'<<itor->first<<'\n';
                }
                m_hash_vecdnb_exon.erase(l_id);
            }
        }
        *pout<<sstrout.str();
    }
    pout->flush();
    if(m_strout != "stdout")
    {
        delete pout;
    }
}

void geftogem::bgeftocgem(const string &strmask, const string &strbgef)
{
    m_bin = 1;
    readBgef(strbgef);
    getdnb();
    readmask(strmask);

    if(m_bexon && m_boutexon)
    {
        cgef2gem_exon();
    }
    else
    {
        cgef2gem();
    }
}

void geftogem::readmask(const string &strmask)
{
    cv::Mat img;
    tifread(img, strmask);
    if (img.empty()) {
        reportErrorCode2File(errorCode::E_INVALIDPARAM, "read mask file error ");
        std::cout << "read mask file error" << std::endl;
        exit(-1);
    }
    assert(!img.empty());
    assert(img.rows == m_max_y - m_min_y+1);
    assert(img.cols == m_max_x - m_min_x+1);

    cv::Mat stats, fill_points, centroids;
    uint32_t labelcnt = connectedComponentsWithStats(img, fill_points, stats, centroids);

    int x, y, w, h;
    for(uint32_t i=1;i<labelcnt;i++)
    {
        x = stats.at<int>(i, cv::CC_STAT_LEFT);
        y = stats.at<int>(i, cv::CC_STAT_TOP);
        w = stats.at<int>(i, cv::CC_STAT_WIDTH);
        h = stats.at<int>(i, cv::CC_STAT_HEIGHT);

        cellmat cm;
        cm.offsetx = 0;
        cm.offsety = 0;

        for(uint32_t m=y;m<y+h;m++)
        {
            for(uint32_t n=x;n<x+w;n++)
            {
                if(fill_points.at<unsigned int>(m,n) == i)
                {
                    cm.vecP.emplace_back(n,m);
                }
            }
        }

        m_hash_cellpoint.emplace(i-1, std::move(cm));
    }
}