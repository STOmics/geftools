#include "bgefCreater.h"
#include "timer.h"
#include "read_task.h"


bgefCreater::bgefCreater(int thcnt):m_threadcnt(thcnt)
{
    m_tpoolPtr = new ThreadPool(m_threadcnt);
}

bgefCreater::~bgefCreater()
{
    delete m_tpoolPtr;
}

void bgefCreater::createBgef(const string &strin, int maskbin, const string &strmask, const string &strout)
{
    vector<Gene> vecgene;
    vector<Expression> vecgExp;
    vector<uint8_t> vecexon;

    m_bin = maskbin;
    tifread(m_outimg, strmask);

    if(H5Fis_hdf5(strin.c_str()))
    {
        readbgef(strin);
        vecgene.reserve(m_genencnt);
        vecgExp.reserve(m_geneexpcnt);
        if(m_bexon)
        {
            vecexon.reserve(m_geneexpcnt);
        }
        getmaskgenedata_bgef(vecgene, vecgExp, vecexon);
        free(m_genePtr);
        free(m_expPtr);
    }
    else
    {
        readgem(strin);
        vecgene.reserve(m_genencnt);
        vecgExp.reserve(m_geneexpcnt);
        if(m_bexon)
        {
            vecexon.reserve(m_geneexpcnt);
        }
        getmaskgenedata_gem(vecgene, vecgExp, vecexon);
    }
    
    writebgef(vecgene, vecgExp, vecexon, strout);
}

void bgefCreater::getStereoData(const string &strin, int maskbin, const string &strmask, vector<string> &vec_gene, vector<unsigned long long> &uniq_cells,
                                vector<unsigned int> &cell_ind, vector<unsigned int> &gene_ind, vector<unsigned int> &count)
{
    m_bin = maskbin;
    tifread(m_outimg, strmask);
    if(H5Fis_hdf5(strin.c_str()))
    {
        readbgef(strin);
        for(uint32_t i=0;i<m_genencnt;i++)
        {
            bgefmaskTask *ptask = new bgefmaskTask(i, this);
            m_tpoolPtr->addTask(ptask);
        }

        unsigned long long uniq_cell_id;
        uint32_t index = 0,gid = 0;
        std::unordered_map<unsigned long long, uint32_t> hash_map;
        uint32_t gcnt = m_genencnt;
        while (gcnt--)
        {
            gdata *ptr = m_queue.getPtr();
            if(ptr->vecExp.size())
            {
                vec_gene.emplace_back(m_genePtr[ptr->geneid].gene);
                for(uint32_t idx : ptr->vecExp)
                {
                    uniq_cell_id = m_expPtr[idx].x;
                    uniq_cell_id = (uniq_cell_id << 32) | m_expPtr[idx].y;

                    if(hash_map.find(uniq_cell_id) != hash_map.end())
                    {
                        cell_ind.push_back(hash_map[uniq_cell_id]);
                    }
                    else
                    {
                        cell_ind.push_back(index);
                        uniq_cells.emplace_back(uniq_cell_id);
                        hash_map.emplace(uniq_cell_id, index++);
                    }
                    count.push_back(m_expPtr[idx].count);
                    gene_ind.push_back(gid);
                }
                gid++;
            }
            
            delete ptr;
        }
        free(m_genePtr);
        free(m_expPtr);
    }
    else
    {
        readgem(strin);
        for(uint32_t i=0;i<m_genencnt;i++)
        {
            gemmaskTask *ptask = new gemmaskTask(i, this);
            m_tpoolPtr->addTask(ptask);
        }

        unsigned long long uniq_cell_id;
        uint32_t index = 0,gid = 0;
        std::unordered_map<unsigned long long, uint32_t> hash_map;
        uint32_t gcnt = m_genencnt;
        while (gcnt--)
        {
            gdata *ptr = m_queue.getPtr();
            if(ptr->vecExp.size())
            {
                string &str_gene = m_vecgenename[ptr->geneid];
                vec_gene.emplace_back(str_gene);
                auto &vecrawexp = m_map_gene_exp[str_gene];
                for(uint32_t idx : ptr->vecExp)
                {
                    uniq_cell_id = vecrawexp[idx].x;
                    uniq_cell_id = (uniq_cell_id << 32) | vecrawexp[idx].y;

                    if(hash_map.find(uniq_cell_id) != hash_map.end())
                    {
                        cell_ind.push_back(hash_map[uniq_cell_id]);
                    }
                    else
                    {
                        cell_ind.push_back(index);
                        uniq_cells.emplace_back(uniq_cell_id);
                        hash_map.emplace(uniq_cell_id, index++);
                    }
                    count.push_back(m_expPtr[idx].count);
                    gene_ind.push_back(gid);
                }
                gid++;
            }
            
            delete ptr;
        }
    }
}

void bgefCreater::readbgef(const string &strin)
{
    timer st(__FUNCTION__);
    hid_t file_id = H5Fopen(strin.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    hsize_t dims[1];
    hid_t gene_did = H5Dopen(file_id, "/geneExp/bin1/gene", H5P_DEFAULT);
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
    H5Sclose(gene_sid);
    H5Dclose(gene_did);


    hid_t exp_did = H5Dopen(file_id, "/geneExp/bin1/expression", H5P_DEFAULT);
    hid_t exp_sid = H5Dget_space(exp_did);
    H5Sget_simple_extent_dims(exp_sid, dims, nullptr);

    m_geneexpcnt = dims[0];

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
        //assert(edims[0] == m_geneexpcnt);
        unsigned int *exonPtr = new unsigned int[edims[0]];
        H5Dread(did, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, exonPtr);
        H5Sclose(sid);
        H5Dclose(did);
        for(int i=0;i<m_geneexpcnt;i++)
        {
            m_expPtr[i].exon = exonPtr[i];
        }
        delete []exonPtr;
    }
    H5Tclose(memtype);
    H5Sclose(exp_sid);
    

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
    printf("minx:%d miny:%d maxx:%d maxy:%d\n", m_min_x, m_min_y, m_max_x, m_max_y);
    H5Aclose(attr);
    H5Dclose(exp_did);


    if(H5Aexists(file_id, "omics"))
    {
        hid_t fattr = H5Aopen(file_id, "omics", H5P_DEFAULT);
        H5Aread(fattr, strtype, m_szomics);
    }
    H5Tclose(strtype);
    H5Fclose(file_id);
    printf("gene:%ld geneexp:%ld\n", m_genencnt, m_geneexpcnt);
    int x_len = (m_max_x - m_min_x)/m_bin*m_bin+1;
    int y_len = (m_max_y - m_min_y)/m_bin*m_bin+1;
    //assert(x_len == m_outimg.cols && y_len == m_outimg.rows);
}

void bgefCreater::readgem(const string &strin)
{
    m_resolution = parseResolutin(strin);
    m_file = gzopen(strin.c_str(), "r");
    gzbuffer(m_file, READLEN);

    int offset_x = 0, offset_y = 0;
    // Process the header lines
    std::string line;
    while (readline(m_file, line))
    {
        if (line[0] == '#')
        {
            // Skip the offset parameter
            if (line.substr(0, 9) == "#OffsetX=")
                offset_x = stoi(line.substr(9));
            else if (line.substr(0, 9) == "#OffsetY=")
                offset_y = stoi(line.substr(9));
            continue;
        }
        if (line.substr(0, 6) == "geneID") break;
    }

    int col = 1;
    for(char ch : line)
    {
        if (ch == '\t')
        {
            col++;
        }
    }
    printf("%s %d\n", line.c_str(), col);
    if(col == 5)
    {
        m_bexon = true;
    }

    for(int i=0;i<m_threadcnt;i++)
    {
        auto *rtask = new ReadTask(m_bexon, m_file, m_range, m_map_gene_exp);
        m_tpoolPtr->addTask(rtask);
    }

    m_tpoolPtr->waitTaskDone();
    gzclose(m_file);

    int minx = m_range[0];
    int miny = m_range[2];
    if (minx != 0 || miny != 0)
    {
        offset_x += minx;
        offset_y += miny;
        for (auto& p : m_map_gene_exp)
        {
            for (auto& g : p.second)
            {
                g.x -= minx;
                g.y -= miny;
            }
            m_geneexpcnt += p.second.size();
            m_vecgenename.emplace_back(p.first);
        }
    }
    else
    {
        for (auto& p : m_map_gene_exp)
        {
            m_geneexpcnt += p.second.size();
            m_vecgenename.emplace_back(p.first);
        }
    }
    m_min_x = offset_x;
    m_min_y = offset_y;
    m_max_x = offset_x + (m_range[1] - m_range[0]);
    m_max_y = offset_y + (m_range[3] - m_range[2]);
    printf("minx:%d miny:%d maxx:%d maxy:%d\n", m_min_x, m_min_y, m_max_x, m_max_y);
    m_genencnt = m_map_gene_exp.size();
    printf("gene:%ld geneexp:%ld\n", m_genencnt, m_geneexpcnt);
}

void bgefCreater::getmaskgenedata_bgef(vector<Gene> &vgene, vector<Expression> &vgExp, vector<uint8_t> &vexon)
{
    timer st(__FUNCTION__);
    for(uint32_t i=0;i<m_genencnt;i++)
    {
        bgefmaskTask *ptask = new bgefmaskTask(i, this);
        m_tpoolPtr->addTask(ptask);
    }

    st.showgap("thread time");
    uint32_t gcnt = m_genencnt;
    uint32_t cnt = 0, offset = 0, ngcnt = 0;
    if(m_bexon)
    {
        while (gcnt--)
        {
            gdata *ptr = m_queue.getPtr();
            if(ptr->vecExp.size())
            {
                for(uint32_t idx : ptr->vecExp)
                {
                    vgExp.emplace_back(m_expPtr[idx]);
                    m_maxExp = std::max(m_maxExp, m_expPtr[idx].count);
                    m_maxExon = std::max(m_maxExon, m_expPtr[idx].exon);
                    vexon.push_back(m_expPtr[idx].exon);
                }
                cnt = ptr->vecExp.size();
                vgene.emplace_back(m_genePtr[ptr->geneid].gene, offset, cnt);
                offset += cnt;
                ngcnt++;
            }
            delete ptr;
        }
    }
    else
    {
        while (gcnt--)
        {
            gdata *ptr = m_queue.getPtr();
            if(ptr->vecExp.size())
            {
                for(uint32_t idx : ptr->vecExp)
                {
                    vgExp.emplace_back(m_expPtr[idx]);
                    m_maxExp = std::max(m_maxExp, m_expPtr[idx].count);
                }
                cnt = ptr->vecExp.size();
                vgene.emplace_back(m_genePtr[ptr->geneid].gene, offset, cnt);
                offset += cnt;
                ngcnt++;
            }
            delete ptr;
        }
    }
    printf("new gcnt:%ld new gexp:%ld\n", ngcnt, offset);
}

void bgefCreater::getmaskgenedata_gem(vector<Gene> &vgene, vector<Expression> &vgExp, vector<uint8_t> &vexon)
{
    timer st(__FUNCTION__);
    for(uint32_t i=0;i<m_genencnt;i++)
    {
        gemmaskTask *ptask = new gemmaskTask(i, this);
        m_tpoolPtr->addTask(ptask);
    }

    st.showgap("thread time");
    uint32_t gcnt = m_genencnt, offset = 0, ngcnt = 0;
    if(m_bexon)
    {
        while (gcnt--)
        {
            gdata *ptr = m_queue.getPtr();
            if(ptr->vecExp.size())
            {
                string &str_gene = m_vecgenename[ptr->geneid];
                vgene.emplace_back(str_gene.c_str(), offset, ptr->vecExp.size());
                offset += ptr->vecExp.size();
                auto &vecrawexp = m_map_gene_exp[str_gene];
                for(uint32_t idx : ptr->vecExp)
                {
                    vgExp.emplace_back(vecrawexp[idx]);
                    m_maxExp = std::max(m_maxExp, vecrawexp[idx].count);
                    m_maxExon = std::max(m_maxExon, vecrawexp[idx].exon);
                    vexon.push_back(vecrawexp[idx].exon);
                }
                ngcnt++;
            }
            delete ptr;
        }
    }
    else
    {
        while (gcnt--)
        {
            gdata *ptr = m_queue.getPtr();
            if(ptr->vecExp.size())
            {
                string &str_gene = m_vecgenename[ptr->geneid];
                vgene.emplace_back(str_gene.c_str(), offset, ptr->vecExp.size());
                offset += ptr->vecExp.size();
                auto &vecrawexp = m_map_gene_exp[str_gene];
                for(uint32_t idx : ptr->vecExp)
                {
                    vgExp.emplace_back(vecrawexp[idx]);
                    m_maxExp = std::max(m_maxExp, vecrawexp[idx].count);
                }
                ngcnt++;
            }
            delete ptr;
        }
    }
    printf("new gcnt:%ld new gexp:%ld\n", ngcnt, offset);
}


void bgefCreater::writebgef(vector<Gene> &vgene, vector<Expression> &vgExp, vector<uint8_t> &vecexon, const string &strout)
{
    timer st(__FUNCTION__);
    //gene exp
    hid_t file_id = H5Fcreate(strout.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hid_t gene_exp = H5Gcreate(file_id, "geneExp", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t gene_exp_bin1 = H5Gcreate(gene_exp, "bin1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hsize_t dims[1] = {vgExp.size()};
    hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(Expression));
    H5Tinsert(memtype, "x", HOFFSET(Expression, x), H5T_NATIVE_INT);
    H5Tinsert(memtype, "y", HOFFSET(Expression, y), H5T_NATIVE_INT);
    H5Tinsert(memtype, "count", HOFFSET(Expression, count), H5T_NATIVE_UINT);
    hid_t filetype = H5Tcreate(H5T_COMPOUND, 8 + 1);
    H5Tinsert(filetype, "x", 0, H5T_STD_I32LE);
    H5Tinsert(filetype, "y", 4, H5T_STD_I32LE);
    H5Tinsert(filetype, "count", 8, H5T_STD_U8LE);
    hid_t exp_did = h5DatasetWrite(gene_exp_bin1, filetype, memtype, "expression", 1, dims, vgExp.data());
    dims[0] = 1;
    h5AttrWrite(exp_did, H5T_STD_I32LE, H5T_NATIVE_INT, "minX", 1, dims, &m_min_x);
    h5AttrWrite(exp_did, H5T_STD_I32LE, H5T_NATIVE_INT, "minY", 1, dims, &m_min_y);
    h5AttrWrite(exp_did, H5T_STD_I32LE, H5T_NATIVE_INT, "maxX", 1, dims, &m_max_x);
    h5AttrWrite(exp_did, H5T_STD_I32LE, H5T_NATIVE_INT, "maxY", 1, dims, &m_max_y);
    h5AttrWrite(exp_did, H5T_STD_I32LE, H5T_NATIVE_INT, "maxExp", 1, dims, &m_maxExp);
    h5AttrWrite(exp_did, H5T_STD_U32LE, H5T_NATIVE_UINT, "resolution", 1, dims, &m_resolution);
    H5Tclose(memtype);
    H5Tclose(filetype);
    H5Dclose(exp_did);

    //gene
    hid_t str32_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(str32_type, 32);
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(Gene));
    H5Tinsert(memtype, "gene", HOFFSET(Gene, gene), str32_type);
    H5Tinsert(memtype, "offset", HOFFSET(Gene, offset), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "count", HOFFSET(Gene, count), H5T_NATIVE_UINT);
    filetype = H5Tcreate(H5T_COMPOUND, 32 + 4 + 4);
    H5Tinsert(filetype, "gene", 0, str32_type);
    H5Tinsert(filetype, "offset", 32, H5T_STD_U32LE);
    H5Tinsert(filetype, "count", 32+4, H5T_STD_U32LE);
    dims[0] = vgene.size();
    hid_t gene_did = h5DatasetWrite(gene_exp_bin1, filetype, memtype, "gene", 1, dims, vgene.data());
    H5Tclose(memtype);
    H5Tclose(filetype);
    H5Dclose(gene_did);

    //exon
    if(m_bexon)
    {
        dims[0] = vecexon.size();
        hid_t exon_did = h5DatasetWrite(gene_exp_bin1, H5T_STD_U8LE, H5T_NATIVE_UCHAR, "exon", 1, dims, vecexon.data());
        dims[0] = 1;
        h5AttrWrite(exon_did, H5T_STD_I32LE, H5T_NATIVE_INT, "maxExon", 1, dims, &m_maxExon);
        H5Dclose(exon_did);
    }

    //attr
    dims[0] = 1;
    int bgefversion = 2;
    h5AttrWrite(file_id, H5T_STD_U32LE, H5T_NATIVE_UINT, "version", 1, dims, &bgefversion);
    h5AttrWrite(file_id, str32_type, str32_type, "omics", 1, dims, m_szomics);
    dims[0] = 3;
    h5AttrWrite(file_id, H5T_STD_I32LE, H5T_NATIVE_INT, "geftool_ver", 1, dims, &GEFVERSION);

    H5Tclose(str32_type);
    H5Gclose(gene_exp_bin1);
    H5Gclose(gene_exp);
    H5Fclose(file_id);
}
