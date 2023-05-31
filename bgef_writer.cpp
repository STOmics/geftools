
#include "bgef_writer.h"
#include <iostream>
#include <cstring>
#include <algorithm>


BgefWriter::BgefWriter(const string &output_filename, bool verbose, bool bexon, const string& stromics) {
    str32_type_ = H5Tcopy(H5T_C_S1);
    H5Tset_size(str32_type_, 32);

    str64_type_ = H5Tcopy(H5T_C_S1);
    H5Tset_size(str64_type_, 64);

    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);

    cerr << "create h5 file: " <<  output_filename << endl;
    file_id_ = H5Fcreate(output_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);

    verbose_ = verbose;
    m_bexon = bexon;
    raw_gef_ = false;

    hsize_t dimsAttr[1] = {1};
    hid_t attr;
    hid_t dataspace_id = H5Screate_simple(1, dimsAttr, nullptr);
    attr = H5Acreate(file_id_, "version", H5T_STD_U32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT, &version);
    H5Sclose(dataspace_id);
    H5Aclose(attr);

    hsize_t gef_dimsAttr[1] = {3};
    hid_t gef_dataspace_id = H5Screate_simple(1, gef_dimsAttr, nullptr);
    hid_t gef_attr = H5Acreate(file_id_, "geftool_ver", H5T_STD_U32LE, gef_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(gef_attr, H5T_NATIVE_UINT, GEFVERSION);
    H5Sclose(gef_dataspace_id);
    H5Aclose(gef_attr);

    hsize_t kind_dims[1] = {1};
    hid_t k_did = H5Screate_simple(1, kind_dims, nullptr);
    hid_t k_attr = H5Acreate(file_id_, "omics", str32_type_, k_did, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(k_attr, str32_type_, stromics.c_str());
    H5Sclose(k_did);
    H5Aclose(k_attr);

    gene_exp_group_id_ = H5Gcreate(file_id_, "geneExp", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    whole_exp_group_id_ = H5Gcreate(file_id_, "wholeExp", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(m_bexon)
    {
        m_wholeExpExon_id = H5Gcreate(file_id_, "wholeExpExon", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
}

void BgefWriter::SetGefArea(float &area)
{
    hsize_t dims[1] = {1};
    hid_t k_did = H5Screate_simple(1, dims, nullptr);
    hid_t k_attr = H5Acreate(file_id_, "gef_area", H5T_IEEE_F32LE, k_did, H5P_DEFAULT,
                             H5P_DEFAULT);
    H5Awrite(k_attr, H5T_NATIVE_FLOAT, &area);
    H5Sclose(k_did);
    H5Aclose(k_attr);
}

BgefWriter::BgefWriter(const string& output_filename, unsigned int raw_gef_version) {
    str32_type_ = H5Tcopy(H5T_C_S1);
    H5Tset_size(str32_type_, 32);

    str64_type_ = H5Tcopy(H5T_C_S1);
    H5Tset_size(str64_type_, 64);

    hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fclose_degree(fapl, H5F_CLOSE_STRONG);

    cerr << "create h5 file: " <<  output_filename << endl;
    file_id_ = H5Fcreate(output_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);

    hsize_t dimsAttr[1] = {1};
    hid_t attr;
    hid_t dataspace_id = H5Screate_simple(1, dimsAttr, nullptr);
    attr = H5Acreate(file_id_, "version", H5T_STD_U32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT, &version);
    H5Sclose(dataspace_id);
    H5Aclose(attr);

    gene_exp_group_id_ = H5Gcreate(file_id_, "geneExp", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    raw_gef_ = true;
}

BgefWriter::~BgefWriter(){
    if (!raw_gef_) {
        H5Gclose(whole_exp_group_id_);
        if(m_bexon) {
            H5Gclose(m_wholeExpExon_id);
        }
        H5Tclose(str32_type_);
        H5Tclose(str64_type_);
    }

    H5Gclose(gene_exp_group_id_);
    H5Fclose(file_id_);
}

bool BgefWriter::storeGene(vector<Expression>& exps, vector<Gene>& genes, DnbAttr& dnbAttr, unsigned int maxexp, int binsize)
{
    char buf[32]={0};
    sprintf(buf, "bin%d", binsize);
    hid_t gene_exp_bin_group_id = H5Gcreate(gene_exp_group_id_, buf, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Create expression compound
    int rank = 1;
    hsize_t dims[1] = {exps.size()};

    hid_t memtype, filetype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(Expression));
    H5Tinsert(memtype, "x", HOFFSET(Expression, x), H5T_NATIVE_INT);
    H5Tinsert(memtype, "y", HOFFSET(Expression, y), H5T_NATIVE_INT);
    H5Tinsert(memtype, "count", HOFFSET(Expression, count), H5T_NATIVE_UINT);

    if (maxexp > USHRT_MAX)
    {
        filetype = H5Tcreate(H5T_COMPOUND, 8 + 4);
        H5Tinsert(filetype, "x", 0, H5T_STD_I32LE);
        H5Tinsert(filetype, "y", 4, H5T_STD_I32LE);
        H5Tinsert(filetype, "count", 8, H5T_STD_U32LE);
    }
    else if (maxexp > UCHAR_MAX)
    {
        filetype = H5Tcreate(H5T_COMPOUND, 8 + 2);
        H5Tinsert(filetype, "x", 0, H5T_STD_I32LE);
        H5Tinsert(filetype, "y", 4, H5T_STD_I32LE);
        H5Tinsert(filetype, "count", 8, H5T_STD_U16LE);
    }
    else
    {
        filetype = H5Tcreate(H5T_COMPOUND, 8 + 1);
        H5Tinsert(filetype, "x", 0, H5T_STD_I32LE);
        H5Tinsert(filetype, "y", 4, H5T_STD_I32LE);
        H5Tinsert(filetype, "count", 8, H5T_STD_U8LE);       
    }

    hid_t dataspace_id = H5Screate_simple(rank, dims, nullptr);
    hid_t dataset_id = H5Dcreate(gene_exp_bin_group_id, "expression", filetype, dataspace_id, H5P_DEFAULT,
        H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &exps[0]);

    // Create expression attribute
    // ExpressionAttr expAttr{dnbAttr.min_x, dnbAttr.min_y, dnbAttr.min_x+(dnbAttr.len_x-1)*binsize, dnbAttr.min_y+(dnbAttr.len_y-1)*binsize, maxexp};
    ExpressionAttr expAttr{dnbAttr.min_x, dnbAttr.min_y, dnbAttr.max_x, dnbAttr.max_y, maxexp};
    hsize_t dimsAttr[1] = {1};
    hid_t attr;
    dataspace_id = H5Screate_simple(1, dimsAttr, nullptr);
    attr = H5Acreate(dataset_id, "minX", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT, &expAttr.min_x);
    attr = H5Acreate(dataset_id, "minY", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT, &expAttr.min_y);
    attr = H5Acreate(dataset_id, "maxX", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT, &expAttr.max_x);
    attr = H5Acreate(dataset_id, "maxY", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT, &expAttr.max_y);
    attr = H5Acreate(dataset_id, "maxExp", H5T_STD_U32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT, &expAttr.max_exp);
    attr = H5Acreate(dataset_id, "resolution", H5T_STD_U32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT, &resolution_);

    // Create gene compound
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(Gene));
    H5Tinsert(memtype, "gene", HOFFSET(Gene, gene), str64_type_);
    H5Tinsert(memtype, "offset", HOFFSET(Gene, offset), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "count", HOFFSET(Gene, count), H5T_NATIVE_UINT);

    filetype = H5Tcreate(H5T_COMPOUND, 64 + 4 + 4);
    H5Tinsert(filetype, "gene", 0, str64_type_);
    H5Tinsert(filetype, "offset", 64, H5T_STD_U32LE);
    H5Tinsert(filetype, "count", 64+4, H5T_STD_U32LE);

    dims[0] = genes.size();
    dataspace_id = H5Screate_simple(rank, dims, nullptr);
    dataset_id = H5Dcreate(gene_exp_bin_group_id, "gene", filetype, dataspace_id, H5P_DEFAULT,
        H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &genes[0]);

    // Create gene attribute
    // m_dataspace_id = H5Screate_simple(1, dimsAttr, NULL);
    // attr = H5Acreate(m_dataset_id, "number", H5T_STD_U32LE, m_dataspace_id, H5P_DEFAULT,
    //     H5P_DEFAULT);
    // H5Awrite(attr, H5T_NATIVE_UINT, &dims[0]);

    H5Aclose(attr);
    H5Tclose(memtype);
    H5Tclose(filetype);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Gclose(gene_exp_bin_group_id);
    return true;
}

bool BgefWriter::storeGeneExon(vector<Expression>& exps, unsigned int maxexon, int binsize)
{
    if(!m_bexon) return false;
    char buf[32]={0};
    sprintf(buf, "bin%d", binsize);
    hid_t gene_exp_bin_group_id = H5Gopen(gene_exp_group_id_, buf, H5P_DEFAULT);

    hsize_t dims[1] = {exps.size()};
    hid_t dataspace_id = H5Screate_simple(1, dims, nullptr);
    hid_t exon_did = 0;
    if(maxexon > USHRT_MAX)
    {
        exon_did = H5Dcreate(gene_exp_bin_group_id, "exon", H5T_STD_U32LE, dataspace_id, H5P_DEFAULT,
            H5P_DEFAULT, H5P_DEFAULT);
    }
    else if(maxexon > UCHAR_MAX)
    {
        exon_did = H5Dcreate(gene_exp_bin_group_id, "exon", H5T_STD_U16LE, dataspace_id, H5P_DEFAULT,
            H5P_DEFAULT, H5P_DEFAULT);
    }
    else
    {
        exon_did = H5Dcreate(gene_exp_bin_group_id, "exon", H5T_STD_U8LE, dataspace_id, H5P_DEFAULT,
            H5P_DEFAULT, H5P_DEFAULT);
    }
    vector<uint32_t> vecexon;
    for(Expression &exp : exps)
    {
        vecexon.push_back(exp.exon);
    }
    H5Dwrite(exon_did, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vecexon.data());

    hsize_t dimsAttr[1] = {1};
    hid_t a_sid = H5Screate_simple(1, dimsAttr, nullptr);
    hid_t attr = H5Acreate(exon_did, "maxExon", H5T_STD_I32LE, a_sid, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT, &maxexon);

    H5Aclose(attr);
    H5Sclose(a_sid);
    H5Sclose(dataspace_id);
    H5Dclose(exon_did);
    return true;
}

bool BgefWriter::storeDnb(DnbMatrix & dnb_matrix, int binsize){
    // Add compound dataset
    hid_t memtype, filetype;

    if (binsize == 1)
    {
        memtype = H5Tcreate(H5T_COMPOUND, sizeof(BinStatUS));
        H5Tinsert(memtype, "MIDcount", HOFFSET(BinStatUS, mid_count), H5T_NATIVE_USHORT);
        H5Tinsert(memtype, "genecount", HOFFSET(BinStatUS, gene_count), H5T_NATIVE_USHORT);

        filetype = H5Tcreate(H5T_COMPOUND, 4);
        H5Tinsert(filetype, "MIDcount", 0, H5T_STD_U16LE);
        H5Tinsert(filetype, "genecount", 2, H5T_STD_U16LE);
    }
    else
    {
        memtype = H5Tcreate(H5T_COMPOUND, sizeof(BinStat));
        H5Tinsert(memtype, "MIDcount", HOFFSET(BinStat, mid_count), H5T_NATIVE_UINT);
        H5Tinsert(memtype, "genecount", HOFFSET(BinStat, gene_count), H5T_NATIVE_USHORT);

        if (dnb_matrix.dnb_attr.max_mid > USHRT_MAX)
        {
            filetype = H5Tcreate(H5T_COMPOUND, 6);
            H5Tinsert(filetype, "MIDcount", 0, H5T_STD_U32LE);
            H5Tinsert(filetype, "genecount", 4, H5T_STD_U16LE);
        }
        else if (dnb_matrix.dnb_attr.max_mid > UCHAR_MAX)
        {
            filetype = H5Tcreate(H5T_COMPOUND, 4);
            H5Tinsert(filetype, "MIDcount", 0, H5T_STD_U16LE);
            H5Tinsert(filetype, "genecount", 2, H5T_STD_U16LE);
        }
        else
        {
            filetype = H5Tcreate(H5T_COMPOUND, 3);
            H5Tinsert(filetype, "MIDcount", 0, H5T_STD_U8LE);
            H5Tinsert(filetype, "genecount", 1, H5T_STD_U16LE);
        }
    }

    unsigned int y_len = dnb_matrix.dnb_attr.len_y;
    unsigned int x_len = dnb_matrix.dnb_attr.len_x;
    hsize_t dims[2];
    dims[0] = x_len;
    dims[1] = y_len;

    char dataName[32]={0};
    sprintf(dataName, "bin%d", binsize);
    // printf("%s %lld\n",dataName, (unsigned long)(y_len)*x_len);
    hid_t dataspace_id = H5Screate_simple(2, dims, nullptr);
    hid_t dataset_id = H5Dcreate(whole_exp_group_id_, dataName, filetype, dataspace_id,
                                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    if (binsize == 1)
        H5Dwrite(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, dnb_matrix.pmatrix_us);
    else
        H5Dwrite(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, dnb_matrix.pmatrix);

    // Create attribute
    hsize_t dimsAttr[1] = {1};
    hid_t attr;
    dataspace_id = H5Screate_simple(1, dimsAttr, nullptr);
    unsigned int len_x = dnb_matrix.dnb_attr.len_x*binsize;
    unsigned int len_y = dnb_matrix.dnb_attr.len_y*binsize;
    attr = H5Acreate(dataset_id, "minX", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT, &dnb_matrix.dnb_attr.min_x);
    attr = H5Acreate(dataset_id, "lenX", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT, &len_x);
    attr = H5Acreate(dataset_id, "minY", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT, &dnb_matrix.dnb_attr.min_y);
    attr = H5Acreate(dataset_id, "lenY", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT, &len_y);
    attr = H5Acreate(dataset_id, "maxMID", H5T_STD_U32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT, &dnb_matrix.dnb_attr.max_mid);
    attr = H5Acreate(dataset_id, "maxGene", H5T_STD_U32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT, &dnb_matrix.dnb_attr.max_gene);
    attr = H5Acreate(dataset_id, "number", H5T_STD_U64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_ULONG, &dnb_matrix.dnb_attr.number);
    attr = H5Acreate(dataset_id, "resolution", H5T_STD_U32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT, &resolution_);

    H5Aclose(attr);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Tclose(memtype);
    H5Tclose(filetype);

    return true;
}

bool BgefWriter::storeWholeExon(DnbMatrix & dnb_matrix, int binsize)
{
    if(!m_bexon) return false;
    char dataName[32]={0};
    sprintf(dataName, "bin%d", binsize);
    hsize_t dims[2] = {dnb_matrix.dnb_attr.len_x, dnb_matrix.dnb_attr.len_y};
    hid_t dataspace_id = H5Screate_simple(2, dims, nullptr);
    hid_t dataset_id = 0;
    if(dnb_matrix.dnb_attr.max_exon > USHRT_MAX)
    {
        dataset_id= H5Dcreate(m_wholeExpExon_id, dataName, H5T_STD_U32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    else if(dnb_matrix.dnb_attr.max_exon > UCHAR_MAX)
    {
        dataset_id= H5Dcreate(m_wholeExpExon_id, dataName, H5T_STD_U16LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    else
    {
        dataset_id= H5Dcreate(m_wholeExpExon_id, dataName, H5T_STD_U8LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    if (binsize == 1)
        H5Dwrite(dataset_id, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, dnb_matrix.pexon16);
    else
        H5Dwrite(dataset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, dnb_matrix.pexon32);

    // Create attribute
    hsize_t dimsAttr[1] = {1};
    hid_t attr_sid = H5Screate_simple(1, dimsAttr, nullptr);
    hid_t attr = H5Acreate(dataset_id, "maxExon", H5T_STD_U32LE, attr_sid, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT, &(dnb_matrix.dnb_attr.max_exon));
    
    H5Sclose(attr_sid);
    H5Aclose(attr);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    return true;
}

bool BgefWriter::storeStat(vector<GeneStat>& geneStat) const
{
    // printf("create group %s\n", buf);
    hid_t group_id = H5Gcreate(file_id_, "stat", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    if (geneStat.empty()) return false;

    // Create gene stat compound
    int rank = 1;
    hsize_t dims[1] = {geneStat.size()};

    hid_t memtype, filetype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(GeneStat));
    H5Tinsert(memtype, "gene", HOFFSET(GeneStat, gene), str64_type_);
    H5Tinsert(memtype, "MIDcount", HOFFSET(GeneStat, mid_count), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "E10", HOFFSET(GeneStat, E10), H5T_NATIVE_FLOAT);

    filetype = H5Tcreate(H5T_COMPOUND, 64 + 4 + 4);
    H5Tinsert(filetype, "gene", 0, str64_type_);
    H5Tinsert(filetype, "MIDcount", 64, H5T_STD_U32LE);
    H5Tinsert(filetype, "E10", 64 + 4, H5T_IEEE_F32LE);

    hid_t dataspace_id = H5Screate_simple(rank, dims, nullptr);
    hid_t dataset_id = H5Dcreate(group_id, "gene", filetype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, filetype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &geneStat[0]);

    // Create gene stat attribute
    float minE10 = geneStat[0].E10, maxE10 = geneStat[0].E10, cutoff = 0.1;
    for (auto& stat : geneStat)
    {
        minE10 = std::min(stat.E10, minE10);
        maxE10 = std::max(stat.E10, maxE10);
    }
    
    hsize_t dimsAttr[1] = {1};
    hid_t attr; 
    dataspace_id = H5Screate_simple(1, dimsAttr, nullptr);
    attr = H5Acreate(dataset_id, "minE10", H5T_IEEE_F32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_FLOAT, &minE10);
    attr = H5Acreate(dataset_id, "maxE10", H5T_IEEE_F32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_FLOAT, &maxE10);
    attr = H5Acreate(dataset_id, "cutoff", H5T_IEEE_F32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_FLOAT, &cutoff);

    H5Aclose(attr);
    H5Tclose(memtype);
    H5Tclose(filetype);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Gclose(group_id);

    return true;
}

unsigned int BgefWriter::getResolution() const {
    return resolution_;
}

void BgefWriter::setResolution(unsigned int resolution) {
    resolution_ = resolution;
}

void BgefWriter::StoreRawGef(Expression *exps, unsigned int exp_size, ExpressionAttr &exp_attr, Gene *genes, 
                             unsigned int gene_cnt, unsigned int *exons, unsigned int maxexon) {
    // Create expression compound
    hid_t gene_exp_bin_group_id = H5Gcreate(gene_exp_group_id_, "bin1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    int rank = 1;
    hsize_t dims[1] = {exp_size};

    hid_t memtype, filetype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(Expression));
    H5Tinsert(memtype, "x", HOFFSET(Expression, x), H5T_NATIVE_INT);
    H5Tinsert(memtype, "y", HOFFSET(Expression, y), H5T_NATIVE_INT);
    H5Tinsert(memtype, "count", HOFFSET(Expression, count), H5T_NATIVE_UINT);

    if (exp_attr.max_exp > USHRT_MAX)
    {
        filetype = H5Tcreate(H5T_COMPOUND, 8 + 4);
        H5Tinsert(filetype, "x", 0, H5T_STD_I32LE);
        H5Tinsert(filetype, "y", 4, H5T_STD_I32LE);
        H5Tinsert(filetype, "count", 8, H5T_STD_U32LE);
    }
    else if (exp_attr.max_exp > UCHAR_MAX)
    {
        filetype = H5Tcreate(H5T_COMPOUND, 8 + 2);
        H5Tinsert(filetype, "x", 0, H5T_STD_I32LE);
        H5Tinsert(filetype, "y", 4, H5T_STD_I32LE);
        H5Tinsert(filetype, "count", 8, H5T_STD_U16LE);
    }
    else
    {
        filetype = H5Tcreate(H5T_COMPOUND, 8 + 1);
        H5Tinsert(filetype, "x", 0, H5T_STD_I32LE);
        H5Tinsert(filetype, "y", 4, H5T_STD_I32LE);
        H5Tinsert(filetype, "count", 8, H5T_STD_U8LE);       
    }      

    hid_t dataspace_id = H5Screate_simple(rank, dims, nullptr);
    hid_t dataset_id = H5Dcreate(gene_exp_bin_group_id, "expression", filetype, dataspace_id, H5P_DEFAULT,
        H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, exps);

    // Create expression attribute
    hsize_t dimsAttr[1] = {1};
    hid_t attr;
    dataspace_id = H5Screate_simple(1, dimsAttr, nullptr);
    attr = H5Acreate(dataset_id, "minX", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT, &exp_attr.min_x);
    attr = H5Acreate(dataset_id, "minY", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT, &exp_attr.min_y);
    attr = H5Acreate(dataset_id, "maxX", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT, &exp_attr.max_x);
    attr = H5Acreate(dataset_id, "maxY", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT, &exp_attr.max_y);
    attr = H5Acreate(dataset_id, "maxExp", H5T_STD_U32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT, &exp_attr.max_exp);
    attr = H5Acreate(dataset_id, "resolution", H5T_STD_U32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT, &exp_attr.resolution);

    // Create gene compound
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(Gene));
    H5Tinsert(memtype, "gene", HOFFSET(Gene, gene), str64_type_);
    H5Tinsert(memtype, "offset", HOFFSET(Gene, offset), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "count", HOFFSET(Gene, count), H5T_NATIVE_UINT);

    filetype = H5Tcreate(H5T_COMPOUND, 64 + 4 + 4);
    H5Tinsert(filetype, "gene", 0, str64_type_);
    H5Tinsert(filetype, "offset", 64, H5T_STD_U32LE);
    H5Tinsert(filetype, "count", 64+4, H5T_STD_U32LE);

    dims[0] = gene_cnt;
    dataspace_id = H5Screate_simple(rank, dims, nullptr);
    dataset_id = H5Dcreate(gene_exp_bin_group_id, "gene", filetype, dataspace_id, H5P_DEFAULT,
        H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, genes);

    // Create gene attribute
    // m_dataspace_id = H5Screate_simple(1, dimsAttr, NULL);
    // attr = H5Acreate(m_dataset_id, "number", H5T_STD_U32LE, m_dataspace_id, H5P_DEFAULT,
    //     H5P_DEFAULT);
    // H5Awrite(attr, H5T_NATIVE_UINT, &dims[0]);

    H5Aclose(attr);
    H5Tclose(memtype);
    H5Tclose(filetype);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    if (exons) {
        hsize_t dims_exon[1] = {exp_size};
        dataspace_id = H5Screate_simple(1, dims_exon, nullptr);
        hid_t exon_did = 0;
        if (maxexon > USHRT_MAX) {
            exon_did =
                H5Dcreate(gene_exp_bin_group_id, "exon", H5T_STD_U32LE,
                          dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        } else if (maxexon > UCHAR_MAX) {
            exon_did =
                H5Dcreate(gene_exp_bin_group_id, "exon", H5T_STD_U16LE,
                          dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        } else {
            exon_did =
                H5Dcreate(gene_exp_bin_group_id, "exon", H5T_STD_U8LE,
                          dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        }
        H5Dwrite(exon_did, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                 exons);

        hsize_t dims_exon_attr[1] = {1};
        hid_t a_sid = H5Screate_simple(1, dims_exon_attr, nullptr);
        hid_t exon_attr = H5Acreate(exon_did, "maxExon", H5T_STD_I32LE, a_sid,
                                    H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(exon_attr, H5T_NATIVE_UINT, &maxexon);

        H5Aclose(exon_attr);
        H5Sclose(a_sid);
        H5Sclose(dataspace_id);
        H5Dclose(exon_did);
    }
    H5Gclose(gene_exp_bin_group_id);

}
