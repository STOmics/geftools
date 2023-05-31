//
// Created by huangzhibo on 2021/12/15.
//

#include "bgef_reader.h"
#include "khash.h"
#include "bin_task.h"
#include "dnb_merge_task.h"
#include "getdataTask.h"
#include "getBgefExpTask.h"

#define FILE_HEADER "#FileFormat=GEMv%d.%d\n" \
"#SortedBy=None\n" \
"#BinSize=%d\n" \
"#STOmicsChip=%s\n" \
"#OffsetX=%d\n" \
"#OffsetY=%d\n" \
"geneID\tx\ty\tMIDCount\n"

#define FILE_HEADER_EXON "#FileFormat=GEMv%d.%d\n" \
"#SortedBy=None\n" \
"#BinSize=%d\n" \
"#STOmicsChip=%s\n" \
"#OffsetX=%d\n" \
"#OffsetY=%d\n" \
"geneID\tx\ty\tMIDCount\texon\n"

KHASH_MAP_INIT_INT64(m64, unsigned int)

std::mutex getdataTask::m_mtx;

BgefReader::BgefReader(const string &filename, int bin_size, int n_thread, bool verbose) {
    file_id_ = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    printf("path:%s bin:%d\n", filename.c_str(), bin_size);
    if(file_id_ < 0)
    {
        printf("H5Fopen error\n");
        reportErrorCode2File(errorCode::E_FILEOPENERROR, "H5Fopen error ");
        exit(1);
    }
    bin_size_ = bin_size;
    verbose_ = verbose;
    n_thread_ = n_thread;

    char dname[128] = {0};
    sprintf(dname, "/geneExp/bin1/exon");
    if(H5Lexists(file_id_, dname, H5P_DEFAULT)>0)
    {
        m_bexon = true;
    }
    else
    {
        printf("%s is not exist\n", dname);
    }

    char binName[128]={0};
    sprintf(binName, "/geneExp/bin%d", bin_size_);
    if(H5Lexists(file_id_, binName, H5P_DEFAULT)>0){
        openExpressionSpace(bin_size_);
        openGeneSpace(bin_size_);
        if(m_bexon)
        {
            openExonSpace(bin_size_);
        }
    }else{
        openExpressionSpace(1);
        openGeneSpace(1);
        if(m_bexon)
        {
            openExonSpace(1);
        }
        generateGeneExp(bin_size_, n_thread);
    }

    hid_t attr;
    attr = H5Aopen(file_id_, "version", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &version_);
    H5Aclose(attr);

    if (H5Lexists(file_id_, "gef_area", H5P_DEFAULT) > 0) {
        hid_t area_attr = H5Aopen(file_id_, "gef_area", H5P_DEFAULT);
        H5Aread(area_attr, H5T_NATIVE_FLOAT, &gef_area_);
        H5Aclose(area_attr);
    } else {
        gef_area_ = 0.0f;
    }
}

BgefReader::~BgefReader() {
    if(genes_ != nullptr)
        free(genes_);
    if(cell_indices_ != nullptr)
        free(cell_indices_);
    if(expressions_ != nullptr)
        free(expressions_);
    if(reduce_expressions_ != nullptr)
        free(reduce_expressions_);
    
    H5Dclose(exp_dataset_id_);
    H5Sclose(exp_dataspace_id_);
    H5Dclose(gene_dataset_id_);
    H5Sclose(gene_dataspace_id_);
    if(whole_exp_dataset_id_ > 0)
        H5Dclose(whole_exp_dataset_id_);
    if(whole_exp_dataspace_id_ > 0)
        H5Sclose(whole_exp_dataspace_id_);
    if(m_exon_did)
        H5Dclose(m_exon_did);
    H5Fclose(file_id_);
//    if(opts_ != nullptr){
//        opts_->expressions_.clear();
//        opts_->genes_.clear();
//    }
}

void BgefReader::openExpressionSpace(int bin_size) {
    hsize_t dims[1];
    // Read raw data
    char expName[128]={0};
    sprintf(expName, "/geneExp/bin%d/expression", bin_size);
    exp_dataset_id_ = H5Dopen(file_id_, expName, H5P_DEFAULT);
    if (exp_dataset_id_ < 0)
    {
        cerr<<"failed open dataset: "<<expName<<endl;
        return;
    }
    exp_dataspace_id_ = H5Dget_space(exp_dataset_id_);
    H5Sget_simple_extent_dims(exp_dataspace_id_, dims, nullptr);
    expression_num_ = dims[0];

}

void BgefReader::openExonSpace(int bin_size)
{
    char expName[128]={0};
    sprintf(expName, "/geneExp/bin%d/exon", bin_size);
    m_exon_did = H5Dopen(file_id_, expName, H5P_DEFAULT);
    if (exp_dataset_id_ < 0)
    {
        cerr<<"failed open dataset: "<<expName<<endl;
        return;
    }
}

void BgefReader::openGeneSpace(int bin_size) {
    hsize_t dims[1];

    // Read index
    char idxName[128]={0};
    sprintf(idxName, "/geneExp/bin%d/gene", bin_size);
    gene_dataset_id_ = H5Dopen(file_id_, idxName, H5P_DEFAULT);
    if (gene_dataset_id_ < 0)
    {
        cerr<<"failed open dataset: "<<idxName<<endl;
        return;
    }
    gene_dataspace_id_ = H5Dget_space(gene_dataset_id_);
    H5Sget_simple_extent_dims(gene_dataspace_id_, dims, nullptr);
    gene_num_ = dims[0];
}

void BgefReader::openWholeExpSpace() {
    hsize_t dims[2];

    // Read index
    char idxName[128]={0};
    sprintf(idxName, "/wholeExp/bin%d", bin_size_);
//    if(!H5Lexists(file_id_, idxName, H5P_DEFAULT)){
//        generateWholeExp(bin_size_, n_thread_);
//    }
    whole_exp_dataset_id_ = H5Dopen(file_id_, idxName, H5P_DEFAULT);
    if (whole_exp_dataset_id_ < 0)
    {
        cerr<<"failed open wholeExp dataset: "<<idxName<<endl;
        return;
    }
    whole_exp_dataspace_id_ = H5Dget_space(whole_exp_dataset_id_);
    H5Sget_simple_extent_dims(whole_exp_dataspace_id_, dims, nullptr);

    whole_exp_matrix_shape_[0] = dims[0];
    whole_exp_matrix_shape_[1] = dims[1];
}

//hash较慢，应使用buildCellInfo2
void BgefReader::buildCellInfo() {
    unsigned long cprev=clock();
    if(cell_num_ != 0 && cell_indices_ != nullptr)
        return;

    hid_t memtype;

    memtype = H5Tcreate(H5T_COMPOUND, sizeof(Coordinate));
    H5Tinsert(memtype, "x", HOFFSET(Coordinate, pos[1]), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "y", HOFFSET(Coordinate, pos[0]), H5T_NATIVE_UINT);

    Coordinate * xy_id;
    xy_id = (Coordinate *) malloc(expression_num_ * sizeof(Coordinate));
    H5Dread(exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, xy_id);

    Coordinate uniq_cell_id{};
    unsigned int index = 0;

    cell_indices_ = (unsigned int *) malloc(expression_num_ * sizeof(unsigned int));

    int absent, is_missing;
    khint_t k;
    khash_t(m64) *h = kh_init(m64);  // allocate a hash table

    for (unsigned long long i = 0; i < expression_num_; ++i) {
        uniq_cell_id = xy_id[i];
        k = kh_get(m64, h, uniq_cell_id.pos_id);
        is_missing = (k == kh_end(h));
        if (!is_missing){
            cell_indices_[i] = kh_value(h, k);
        }else {
            cell_indices_[i] = index;
            cell_pos_.emplace_back(uniq_cell_id);
            k = kh_put(m64, h, uniq_cell_id.pos_id, &absent);  // insert a key to the hash table
            kh_value(h, k) = index;
            ++index;
        }
    }

    cell_num_ = index;
    kh_destroy(m64, h);
    H5Tclose(memtype);
    free(xy_id);
    if(verbose_) printCpuTime(cprev, "buildCellInfo");
}

void BgefReader::buildCellInfo2() {
    unsigned long cprev=clock();
    if(cell_num_ != 0 && cell_indices_ != nullptr)
        return;

    hid_t memtype;
    Coordinate * xy_id;
    xy_id = (Coordinate *) malloc(expression_num_ * sizeof(Coordinate));
    unsigned long cprev2=clock();

    if(expressions_ != nullptr){
        for(unsigned long long i = 0; i < expression_num_; i++){
            xy_id[i].pos[1] = expressions_[i].x;
            xy_id[i].pos[0] = expressions_[i].y;
        }
    }else{
        memtype = H5Tcreate(H5T_COMPOUND, sizeof(Coordinate));
        H5Tinsert(memtype, "x", HOFFSET(Coordinate, pos[1]), H5T_NATIVE_UINT);
        H5Tinsert(memtype, "y", HOFFSET(Coordinate, pos[0]), H5T_NATIVE_UINT);

        H5Dread(exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, xy_id);
        H5Tclose(memtype);
    }

    if(verbose_) printCpuTime(cprev2, "read");

    cell_indices_ = (unsigned int *) malloc(expression_num_ * sizeof(unsigned int));

    auto * exp_index = (unsigned int *) malloc(expression_num_ * sizeof(unsigned int));
    iota(exp_index, exp_index+expression_num_, 0);
    sort(exp_index, exp_index+expression_num_,
         [xy_id](int a,int b){return xy_id[a].pos_id < xy_id[b].pos_id; });

    Coordinate uniq_cell_id{}, pre_xy = xy_id[exp_index[0]];
    cell_pos_.emplace_back(pre_xy);
    unsigned int index = 0;
    cell_indices_[exp_index[0]] = 0;
    for (unsigned long long i = 1; i < expression_num_; ++i) {
        uniq_cell_id = xy_id[exp_index[i]];
        if (uniq_cell_id.pos_id != pre_xy.pos_id){
            cell_pos_.emplace_back(uniq_cell_id);
            ++index;
            pre_xy = uniq_cell_id;
        }
        cell_indices_[exp_index[i]] = index;
    }

    cell_num_ = cell_pos_.size();
    free(exp_index);
    free(xy_id);
    if(verbose_) printCpuTime(cprev, "buildCellInfo2");
}


int BgefReader::getBinSize() const {
    return bin_size_;
}

unsigned int BgefReader::getGeneNum() const {
    return gene_num_;
}

unsigned int BgefReader::getCellNum() {
    unsigned long cprev=clock();
    if(cell_num_ != 0 && cell_indices_ != nullptr)
        return cell_num_;

    buildCellInfo2();
    if(verbose_) printCpuTime(cprev, "getCellNum");
    return cell_num_;
}

unsigned int BgefReader::getExpressionNum() const {
    return expression_num_;
}

ExpressionAttr &BgefReader::getExpressionAttr() {
    if(expression_attr_init_)
        return expression_attr_;
    hid_t attr;
    attr = H5Aopen(exp_dataset_id_, "minX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &(expression_attr_.min_x));
    attr = H5Aopen(exp_dataset_id_, "minY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &(expression_attr_.min_y));
    attr = H5Aopen(exp_dataset_id_, "maxX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &(expression_attr_.max_x));
    attr = H5Aopen(exp_dataset_id_, "maxY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &(expression_attr_.max_y));
    attr = H5Aopen(exp_dataset_id_, "maxExp", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &(expression_attr_.max_exp));
    attr = H5Aopen(exp_dataset_id_, "resolution", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &(expression_attr_.resolution));

    expression_attr_init_ = true;

    H5Aclose(attr);
    return expression_attr_;
}

// void BgefReader::getSparseMatrixIndicesOfExp(vector<unsigned long long> &uniq_cells, unsigned int * cell_ind, unsigned int * count){
//     unsigned long cprev=clock();
//     unsigned long long uniq_cell_id;
//     Expression* expData = getExpression();

//     unsigned int index = 0;

//     int absent, is_missing;
//     khint_t k;
//     khash_t(m64) *h = kh_init(m64);  // allocate a hash table

//     for (int i = 0; i < expression_num_; ++i) {
//         uniq_cell_id = expData[i].x;
//         uniq_cell_id = uniq_cell_id << 32 | expData[i].y;

//         k = kh_get(m64, h, uniq_cell_id);
//         is_missing = (k == kh_end(h));
//         if (!is_missing){
//             cell_ind[i] = kh_value(h, k);
//         }else {
//             cell_ind[i] = index;
//             uniq_cells.push_back(uniq_cell_id);
//             k = kh_put(m64, h, uniq_cell_id, &absent);  // insert a key to the hash table
//             kh_value(h, k) = index;
//             ++index;
//         }

//         count[i] = expData[i].count;
//     }

//     cell_num_ = index;
//     kh_destroy(m64, h);
//     if(verbose_) printCpuTime(cprev, "getSparseMatrixIndicesOfExp");
//     //return uniq_cells;
// }

void BgefReader::getSparseMatrixIndicesOfExp(vector<unsigned long long> &uniq_cells, unsigned int * cell_ind, unsigned int * count)
{
    unsigned long long uniq_cell_id;
    Expression* expData = getExpression();

    //vector<unsigned long long> uniq_cells;
    uniq_cells.reserve(expression_num_);
    uint32_t index = 0;

    std::unordered_map<unsigned long long, uint32_t> hash_map;
    for(uint64_t i=0;i<expression_num_;i++)
    {
        uniq_cell_id = expData[i].x;
        uniq_cell_id = (uniq_cell_id << 32) | expData[i].y;

        if(hash_map.find(uniq_cell_id) != hash_map.end())
        {
            cell_ind[i] = hash_map[uniq_cell_id];
        }
        else
        {
            cell_ind[i] = index;
            uniq_cells.emplace_back(uniq_cell_id);
            hash_map.emplace(uniq_cell_id, index++);
        }
        count[i] = expData[i].count;
    }

    cell_num_ = index;
}

// void BgefReader::getSparseMatrixIndicesOfExp(vector<unsigned long long> &uniq_cells, unsigned int * cell_ind, unsigned int * pcount)
// {
//     unsigned long long uniq_cell_id;
//     Expression* expData = getExpression();

//     uint32_t index = 0;
//     unsigned long long *pcellid = new unsigned long long[expression_num_];
//     ThreadPool thpool(n_thread_);
//     uint32_t count = expression_num_/n_thread_;
//     for(int i=0;i<n_thread_;i++)
//     {
//         Expression* pdata_tmp = expData + i*count;
//         unsigned int *pcount_tmp = pcount + i*count;
//         unsigned long long *pcellid_tmp = pcellid + i*count;
//         if(i == n_thread_-1)
//         {
//             count = expression_num_ - (n_thread_-1)*count;
//         }
//         getBgefExpTask *ptask = new getBgefExpTask(count, pdata_tmp, pcount_tmp, pcellid_tmp);
//         thpool.addTask(ptask);
//     }
//     thpool.waitTaskDone();

//     std::unordered_map<unsigned long long, uint32_t> hash_map;
//     for(uint32_t i=0;i<expression_num_;i++)
//     {
//         if(hash_map.find(pcellid[i]) != hash_map.end())
//         {
//             cell_ind[i] = hash_map[pcellid[i]];
//         }
//         else
//         {
//             cell_ind[i] = index;
//             uniq_cells.emplace_back(pcellid[i]);
//             hash_map.emplace(pcellid[i], index++);
//         }
//     }
//     delete pcellid;
// }

vector<string> BgefReader::getSparseMatrixIndicesOfGene(unsigned int *gene_index) {
    Gene* gene_data = getGene();

    vector<string> uniq_genes;
    unsigned long long exp_len_index = 0;
    for (unsigned int i = 0; i < gene_num_; ++i)
    {
        const char* gene = gene_data[i].gene;
        uniq_genes.emplace_back(gene);
        unsigned int c = gene_data[i].count;
        for (int j = 0; j < c; ++j)
        {
            gene_index[exp_len_index++] = i;
        }
    }

    assert(exp_len_index == expression_num_);

    return uniq_genes;
}


void BgefReader::getSparseMatrixIndicesOfGene(unsigned int *gene_ind, char * gene_names) {
    Gene* gene_data = getGene();

    unsigned long long exp_len_index = 0;
    for (unsigned int i = 0; i < gene_num_; ++i)
    {
        memcpy(&gene_names[i*32], gene_data[i].gene, 32);
        unsigned int c = gene_data[i].count;
        for (int j = 0; j < c; ++j)
        {
            gene_ind[exp_len_index++] = i;
        }
    }
    assert(exp_len_index == expression_num_);
}

Expression *BgefReader::getExpression() {
    if(expressions_ != nullptr)
        return expressions_;

    hid_t memtype;

    memtype = H5Tcreate(H5T_COMPOUND, sizeof(Expression));
    H5Tinsert(memtype, "x", HOFFSET(Expression, x), H5T_NATIVE_INT);
    H5Tinsert(memtype, "y", HOFFSET(Expression, y), H5T_NATIVE_INT);
    H5Tinsert(memtype, "count", HOFFSET(Expression, count), H5T_NATIVE_UINT);

    expressions_ = (Expression *) malloc(expression_num_ * sizeof(Expression));
    H5Dread(exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, expressions_);

    H5Tclose(memtype);

    getGeneExon();
    if(m_exonPtr)
    {
        for(unsigned int i=0;i<expression_num_;i++)
        {
            expressions_[i].exon = m_exonPtr[i];
        }
    }
    return expressions_;
}

Gene *BgefReader::getGene() {
    if(genes_ != nullptr)
        return genes_;

    hid_t memtype, strtype;

    strtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(strtype, 64);

    memtype = H5Tcreate(H5T_COMPOUND, sizeof(Gene));
    H5Tinsert(memtype, "gene", HOFFSET(Gene, gene), strtype);
    H5Tinsert(memtype, "offset", HOFFSET(Gene, offset), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "count", HOFFSET(Gene, count), H5T_NATIVE_UINT);

    genes_ = (Gene*)malloc(gene_num_ * sizeof(Gene));
    H5Dread(gene_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, genes_);
    H5Tclose(strtype);
    H5Tclose(memtype);
    return genes_;
}

cv::Mat BgefReader::getWholeExpMatrix(cv::Rect roi){
    if(whole_exp_matrix_t_.empty())
        cacheWholeExpMatrix();
    return whole_exp_matrix_t_(roi);
}

void BgefReader::getBinGeneExpMap(
        map<unsigned long long int, pair<unsigned int, unsigned short>>& bin_exp_map,
        DnbExpression * dnb_exp_info) {
    unsigned long cprev=clock();
    hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(DnbExpression));
    H5Tinsert(memtype, "x", HOFFSET(DnbExpression, x), H5T_NATIVE_INT);
    H5Tinsert(memtype, "y", HOFFSET(DnbExpression, y), H5T_NATIVE_INT);
    H5Tinsert(memtype, "count", HOFFSET(DnbExpression, count), H5T_NATIVE_USHORT);
    H5Dread(exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, dnb_exp_info);

    Gene * gene_data = getGene();
    unsigned int n = 0;
    for(unsigned int i = 0; i < gene_num_; i++){
        for(unsigned int j = 0; j < gene_data[i].count; j++){
            dnb_exp_info[n++].gene_id = i;
        }
    }
    assert(n == expression_num_);

    sort(dnb_exp_info, dnb_exp_info+expression_num_, expressionComp);

    DnbExpression pre_dnb = dnb_exp_info[0];
    unsigned int offset = 0;
    n = 1;
    unsigned long long int bin_id;
    for (unsigned int i = 1; i < expression_num_; ++i){
        if (dnb_exp_info[i].x == pre_dnb.x && dnb_exp_info[i].y == pre_dnb.y){
            n++;
            continue;
        }

        bin_id = static_cast<unsigned long long int>(pre_dnb.x);
        bin_id = bin_id << 32 | static_cast<unsigned int>(pre_dnb.y);
        bin_exp_map.insert(map<unsigned long long int, pair<unsigned int, unsigned short>>::value_type (
                bin_id, make_pair(offset, n)));
        n=1;
        offset = i;
        pre_dnb = dnb_exp_info[i];
    }

    bin_id = static_cast<unsigned long long int>(pre_dnb.x);
    bin_id = bin_id << 32 | static_cast<unsigned int>(pre_dnb.y);
    bin_exp_map.insert(map<unsigned long long int, pair<unsigned int, unsigned short>>::value_type (
            bin_id, make_pair(offset, n)));

    cell_num_ = bin_exp_map.size();
    H5Tclose(memtype);
    if(verbose_) printCpuTime(cprev, "getBinGeneExpMap");
}

void BgefReader::clear() {
    if(genes_ != nullptr)
        free(genes_);
    if(expressions_ != nullptr)
        free(expressions_);
    whole_exp_matrix_t_.release();
}

void BgefReader::cacheWholeExpMatrix() {
    if(whole_exp_dataset_id_ == 0) openWholeExpSpace();

    hid_t memtype;
    memtype = H5Tcreate(H5T_COMPOUND, 1);
    // genecount的值大于255将读取为255
    whole_exp_matrix_t_ = cv::Mat::zeros(
            (int)whole_exp_matrix_shape_[0], (int)whole_exp_matrix_shape_[1], CV_8UC1);
    H5Tinsert(memtype, "genecount", 0, H5T_NATIVE_UCHAR);
    H5Dread(whole_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, whole_exp_matrix_t_.data);
    whole_exp_matrix_t_ = whole_exp_matrix_t_.t();
    H5Tclose(memtype);
}

void BgefReader::readWholeExpMatrix(string &key, unsigned char *matrix) {
    readWholeExpMatrix(0,
                       0,
                       whole_exp_matrix_shape_[0],
                       whole_exp_matrix_shape_[1],
                       key,
                       matrix);
}

void BgefReader::readWholeExpMatrix(unsigned int offset_x, unsigned int offset_y, unsigned int rows, unsigned int cols,
                                    string &key, unsigned char *matrix) {
    if(whole_exp_dataset_id_ == 0) openWholeExpSpace();
    hsize_t start[2] = {offset_x, offset_y},
            count[2] = {rows, cols},
            offset_out[2] = {0, 0};

    hid_t memtype;
    memtype = H5Tcreate(H5T_COMPOUND, 1);
    // genecount的值大于255将读取为255
    H5Tinsert(memtype, key.c_str(), 0, H5T_NATIVE_UCHAR);

    // Define memory dataspace.
    hid_t memspace = H5Screate_simple(2, count,nullptr);
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, nullptr, count, nullptr);

    H5Sselect_hyperslab(whole_exp_dataspace_id_, H5S_SELECT_SET, start, nullptr, count, nullptr);
    H5Dread(whole_exp_dataset_id_, memtype, memspace, whole_exp_dataspace_id_, H5P_DEFAULT, matrix);

    H5Tclose(memtype);
    H5Sclose(memspace);
}

int BgefReader::getVersion() const {
    return version_;
}

void BgefReader::getGeneNameList(vector<string> & gene_list) {
    Gene * genes = getGene();
    for(unsigned int i = 0; i < gene_num_; i++){
        string name = genes[i].gene;
        gene_list.emplace_back(name);
    }
}

//void BgefReader::getGeneNameList(char *gene_list) {
//    Gene * genes = getGene();
//    for(unsigned int i = 0; i < gene_num_; i++){
//        memcpy(&gene_list[i], genes[i].gene, 32);
//    }
//}


bool BgefReader::expressionComp(const DnbExpression& p1, const DnbExpression& p2) {
    return p1.x < p2.x || (p1.x == p2.x && p1.y < p2.y);
}

const unsigned int *BgefReader::getWholeExpMatrixShape() {
    if(whole_exp_dataset_id_ == 0) openWholeExpSpace();
    return whole_exp_matrix_shape_;
}


int BgefReader::getSparseMatrixIndices(unsigned int *indices, unsigned int *indptr, unsigned int *count) {
    unsigned long cprev=clock();

    if(cell_indices_ == nullptr) buildCellInfo2();
    memcpy(indices, cell_indices_, expression_num_ * sizeof(unsigned int));

    Gene * gene_data = getGene();

    indptr[0] = 0;
    for(unsigned int i = 1; i < gene_num_; i++){
        indptr[i] = gene_data[i].offset;
    }
    indptr[gene_num_] = gene_data[gene_num_-1].offset + gene_data[gene_num_-1].count;

    if(expressions_ != nullptr){
        for(unsigned long long i = 0; i < expression_num_; i++){
            count[i] = expressions_[i].count;
        }
    }else{
        hid_t memtype;
        memtype = H5Tcreate(H5T_COMPOUND, sizeof(unsigned int));
        H5Tinsert(memtype, "count", 0, H5T_NATIVE_UINT);
        H5Dread(exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, count);
        H5Tclose(memtype);
    }

    if(verbose_) printCpuTime(cprev, "getSparseMatrixIndices");
    return 0;
}

int BgefReader::getSparseMatrixIndices2(unsigned int * cell_ind, unsigned int * gene_ind, unsigned int * count){
    unsigned long cprev=clock();
    Gene * gene_data = getGene();

    if(cell_indices_ == nullptr) buildCellInfo2();
    memcpy(cell_ind, cell_indices_, expression_num_ * sizeof(unsigned int));

    if(expressions_ != nullptr){
        for(unsigned long long i = 0; i < expression_num_; i++){
            count[i] = expressions_[i].count;
        }
    }else{
        hid_t memtype;
        memtype = H5Tcreate(H5T_COMPOUND, sizeof(unsigned int));
        H5Tinsert(memtype, "count", 0, H5T_NATIVE_UINT);
        H5Dread(exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, count);
        H5Tclose(memtype);
    }

    unsigned int n = 0;
    for(unsigned int i = 0; i < gene_num_; i++){
        for(unsigned int j = 0; j < gene_data[i].count; j++){
            gene_ind[n++] = i;
        }
    }

    if(verbose_) printCpuTime(cprev, "getSparseMatrixIndices2");
    return 0;
}

void BgefReader::getGeneAndCount(unsigned short * gene_ind, unsigned short * count){
    unsigned long cprev=clock();
    Gene * gene_data = getGene();

    hid_t memtype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(unsigned short));
    H5Tinsert(memtype, "count", 0, H5T_NATIVE_USHORT);
    H5Dread(exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, count);

    unsigned int n = 0;
    for(unsigned int i = 0; i < gene_num_; i++){
        for(unsigned int j = 0; j < gene_data[i].count; j++){
            gene_ind[n++] = i;
        }
    }
    assert(n == expression_num_);

    H5Tclose(memtype);
    if(verbose_) printCpuTime(cprev, "getGeneAndCount");
}

void BgefReader::getCellNameList(unsigned long long int *cell_name_list) {
    memcpy(cell_name_list, cell_pos_.data(), cell_num_ * sizeof(unsigned long long int));
}

unsigned long long int * BgefReader::getCellPos() {
    return reinterpret_cast<unsigned long long int *>(cell_pos_.data());
}

//TODO support restrict gene_list and region
// unsigned int BgefReader::toGem(string &filename, string &sn) {
//     unsigned long cprev = clock();
//     Gene * gene_data = getGene();
//     Expression * expression = getExpression();
//     ExpressionAttr & expression_attr = getExpressionAttr();

//     // Create file header
//     FILE* outhandle;
//     if (filename == "stdout" || filename == "-"){
//         outhandle = stdout;
//     }
//     else{
//         outhandle = fopen(filename.c_str(), "w");
//         if (outhandle == nullptr)
//         {
//             cerr<<"failed create output file: "<< filename <<endl;
//             exit(4);
//         }
//     }

//     size_t pos = 0;
//     if(m_exonPtr)
//     {
//         fprintf(outhandle, FILE_HEADER_EXON, 0, 1, bin_size_, sn.c_str(), expression_attr.min_x, expression_attr.min_y);
//         // Write data line by line
        
//         for (int i = 0; i < gene_num_; ++i)
//         {
//             const char* gene = gene_data[i].gene;
//             size_t end = gene_data[i].offset + gene_data[i].count;
//             while (pos < end){
//                 Expression& exp = expression[pos];
//                 fprintf(outhandle, "%s\t%d\t%d\t%d\t%d\n", gene, exp.x, exp.y, exp.count, exp.exon);
//                 ++pos;
//             }
//         }
//     }
//     else
//     {
//         fprintf(outhandle, FILE_HEADER, 0, 1, bin_size_, sn.c_str(), expression_attr.min_x, expression_attr.min_y);
//         // Write data line by line

//         for (int i = 0; i < gene_num_; ++i)
//         {
//             const char* gene = gene_data[i].gene;
//             size_t end = gene_data[i].offset + gene_data[i].count;
//             while (pos < end){
//                 Expression& exp = expression[pos];
//                 fprintf(outhandle, "%s\t%d\t%d\t%d\n", gene, exp.x, exp.y, exp.count);
//                 ++pos;
//             }
//         }
//     }

//     fclose(outhandle);

//     if(verbose_) printCpuTime(cprev, "toGem");
//     return pos;
// }

void BgefReader::getGeneExpression(unordered_map<string, vector<Expression>> &gene_exp_map,
                                   const vector<int>& regions) {
    if(regions.empty()){
        getGeneExpression(gene_exp_map);
        return;
    }

    int min_x = regions[0];
    int max_x = regions[1];
    int min_y = regions[2];
    int max_y = regions[3];

    Gene * gene = getGene();
    Expression * expression = getExpression();


    for(unsigned int gene_id = 0; gene_id < gene_num_; gene_id++){
        vector<Expression> exps;
        exps.reserve(gene[gene_id].count);
        unsigned int end = gene[gene_id].offset + gene[gene_id].count;
        for(unsigned int i = gene[gene_id].offset; i < end; i++){
            Expression exp = expression[i];
            if(exp.x < min_x || exp.x > max_x || exp.y < min_y || exp.y > max_y){
                continue;
            }
            exp.x -= min_x;
            exp.y -= min_y;

            exps.emplace_back(exp);
        }

        if(exps.empty()) continue;
        gene_exp_map.insert(unordered_map<string, vector<Expression>>::value_type(gene[gene_id].gene, exps));
    }
}

void BgefReader::getGeneExpression(unordered_map<string, vector<Expression>> &gene_exp_map) {
    unsigned long cprev=clock();

    Gene * gene = getGene();
    Expression * expression = getExpression();


    for(unsigned int gene_id = 0; gene_id < gene_num_; gene_id++){
        vector<Expression> exps;
        exps.reserve(gene[gene_id].count);
        unsigned int end = gene[gene_id].offset + gene[gene_id].count;
        for(unsigned int i = gene[gene_id].offset; i < end; i++){
            exps.emplace_back(expression[i]);
        }

        gene_exp_map.insert(unordered_map<string, vector<Expression>>::value_type(gene[gene_id].gene, exps));
    }

    if(verbose_) printCpuTime(cprev, "getGeneExpression");
}

int BgefReader::generateGeneExp(int bin_size, int n_thread) {
    unsigned long cprev=clock();
    ExpressionAttr expression_attr{};

    hid_t attr;
    attr = H5Aopen(exp_dataset_id_, "minX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &(expression_attr.min_x));
    attr = H5Aopen(exp_dataset_id_, "minY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &(expression_attr.min_y));
    attr = H5Aopen(exp_dataset_id_, "maxX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &(expression_attr.max_x));
    attr = H5Aopen(exp_dataset_id_, "maxY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &(expression_attr.max_y));
    attr = H5Aopen(exp_dataset_id_, "maxExp", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &(expression_attr_.max_exp));
    attr = H5Aopen(exp_dataset_id_, "resolution", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT, &(expression_attr_.resolution));

    opts_ = BgefOptions::GetInstance();
    opts_->bin_sizes_.clear();
    opts_->range_.clear();
    opts_->map_gene_exp_.clear();
    opts_->bin_sizes_.emplace_back(bin_size);
    auto& dnb_attr = opts_->dnbmatrix_.dnb_attr;

    opts_->range_ = {expression_attr.min_x, expression_attr.max_x, expression_attr.min_y, expression_attr.max_y};
//    opts->offset_x_ = expression_attr.min_x;
//    opts->offset_y_ = expression_attr.min_y;

//    opts->region_ = std::move(region);
    opts_->verbose_ = verbose_;

//    if(opts->region_.empty()){
        getGeneExpression(opts_->map_gene_exp_);

    dnb_attr.len_x = int((float(expression_attr_.max_x) / bin_size) - (float(expression_attr_.min_x) / bin_size)) + 1;
    dnb_attr.len_y = int((float(expression_attr_.max_y) / bin_size) - (float(expression_attr_.min_y) / bin_size)) + 1;

    expression_attr_.min_x = (expression_attr.min_x / bin_size) * bin_size;
    expression_attr_.min_y = (expression_attr.min_y / bin_size) * bin_size;
    // expression_attr_.max_x = expression_attr_.min_x + (dnb_attr.len_x-1)*bin_size;
    // expression_attr_.max_y = expression_attr_.min_y + (dnb_attr.len_y-1)*bin_size;
    expression_attr_.max_x = (expression_attr.max_x / bin_size) * bin_size;
    expression_attr_.max_y = (expression_attr.max_y / bin_size) * bin_size;

    dnb_attr.min_x = expression_attr_.min_x;
    dnb_attr.min_x = expression_attr_.min_x;
    dnb_attr.max_x = expression_attr_.max_x;
    dnb_attr.max_y = expression_attr_.max_y;


//    }else{
//        bgef_reader.getGeneExpression(opts->map_gene_exp_, opts->region_);
//        unsigned int min_x = opts->region_[0];
//        unsigned int max_x = opts->region_[1];
//        unsigned int min_y = opts->region_[2];
//        unsigned int max_y = opts->region_[3];
//
//        opts->range_ = {expression_attr.min_x + min_x, min(expression_attr.max_x, max_x+expression_attr.min_x),
//                        expression_attr.min_y + min_y, min(expression_attr.max_y, max_y+expression_attr.min_y)};
//        opts->offset_x_ = expression_attr.min_x + min_x;
//        opts->offset_y_ = expression_attr.min_y + min_y;
//    }

    ThreadPool thpool(n_thread);

    auto itor = opts_->map_gene_exp_.begin();
    for(;itor != opts_->map_gene_exp_.end();itor++)
    {
        auto *task = new BinTask(bin_size, itor->first.c_str());
        thpool.addTask(task);
    }

    unsigned int offset = 0;
    unsigned int maxexp = 0;
    int genecnt = 0;
    while (true){
        GeneInfo *pgenedata = opts_->m_geneinfo_queue.getPtr();
        for (auto g : *pgenedata->vecptr){
            g.x *= bin_size;
            g.y *= bin_size;
            opts_->expressions_.push_back(std::move(g));
        }

        opts_->genes_.emplace_back(pgenedata->geneid, offset, static_cast<unsigned int>(pgenedata->vecptr->size()));
        offset += pgenedata->vecptr->size();
        maxexp = std::max(maxexp, pgenedata->maxexp);

        genecnt++;
        if(genecnt == opts_->map_gene_exp_.size()){
            break;
        }
    }

    thpool.waitTaskDone();

    expression_num_ = opts_->expressions_.size();
    gene_num_ = opts_->genes_.size();

    expressions_ = (Expression *) malloc(expression_num_ * sizeof(Expression));
    genes_ = (Gene*)malloc(gene_num_ * sizeof(Gene));

    memcpy(expressions_, &opts_->expressions_[0], expression_num_*sizeof(Expression));
    memcpy(genes_, &opts_->genes_[0], gene_num_*sizeof(Gene));

    //TODO 优化内存，取消opts_
    opts_->expressions_.clear();
    opts_->genes_.clear();

    cprev = printCpuTime(cprev, "generateBinInfo");
    return 0;
}


//TODO
void BgefReader::generateWholeExp(int bin_size, int thread) {

    unsigned long cprev=clock();
    ThreadPool thpool(n_thread_);
    auto& dnb_attr = opts_->dnbmatrix_.dnb_attr;
    auto& dnb_matrix = opts_->dnbmatrix_;

    unsigned long matrix_len = (unsigned long)(dnb_attr.len_x) * dnb_attr.len_y;
    if (bin_size == 1)
    {
        dnb_matrix.pmatrix_us = (BinStatUS*)calloc(matrix_len, sizeof(BinStatUS));
        assert(dnb_matrix.pmatrix_us);
    }
    else {
        dnb_matrix.pmatrix = (BinStat *) calloc(matrix_len, sizeof(BinStat));
    }

    for(int i=0; i < n_thread_; i++){
        auto *task = new DnbMergeTask(opts_->map_gene_exp_.size(), i, bin_size);
        thpool.addTask(task);
    }

    thpool.waitTaskDone();

//    special_bin sbin;
//    std::vector<int> vecdnb;
//    unsigned int x, y;
//    unsigned int y_len = dnbM.dnb_attr.len_y;
//    for(unsigned long i=0;i<matrix_len;i++)
//    {
//        if(dnbM.pmatrix[i].gene_count)
//        {
//            x = i/y_len;
//            y = i%y_len;
//
//            vecdnb.push_back(x*bin);
//            vecdnb.push_back(y*bin);
//            vecdnb.push_back(dnbM.pmatrix[i].mid_count);
//        }
//    }


    printCpuTime(cprev, "generateWholeExp");
}

Expression *BgefReader::getReduceExpression() {
    unsigned int cell_num = getCellNum();
    if(expressions_ == nullptr) getExpression();
    reduce_expressions_ = (Expression *)calloc(cell_num, sizeof(Expression));
    for(unsigned int i = 0; i < expression_num_; i++){
        reduce_expressions_[cell_indices_[i]].x = expressions_[i].x;
        reduce_expressions_[cell_indices_[i]].y = expressions_[i].y;
        reduce_expressions_[cell_indices_[i]].count += expressions_[i].count;
    }
    return reduce_expressions_;
    // if(expressions_ == nullptr) getExpression();
    // map<uint64_t, uint32_t> hash_cell;
    // uint64_t lid = 0;
    // for(unsigned int i = 0; i < expression_num_; i++)
    // {
    //     lid = expressions_[i].x;
    //     lid = lid <<32 | expressions_[i].y;
    //     hash_cell[lid] += expressions_[i].count;
    // }

    // Expression *treduce_expressions_ = (Expression *)calloc(hash_cell.size(), sizeof(Expression));
    // auto itor = hash_cell.begin();
    // uint32_t i = 0;
    // for(;itor!=hash_cell.end();itor++,i++)
    // {
    //     treduce_expressions_[i].x = itor->first>>32;
    //     treduce_expressions_[i].y = itor->first & 0xffffffff;
    //     treduce_expressions_[i].count = itor->second;
    // }
    // return treduce_expressions_;
}

void BgefReader::getfiltereddata(vector<int> &region, vector<string> &genelist,
                                 vector<string> &vec_gene, vector<unsigned long long> &uniq_cells,
                                 vector<unsigned int> &cell_ind, vector<unsigned int> &gene_ind, 
                                 vector<unsigned int> &count)
{
    int min_x = 0, max_x = 0, min_y = 0, max_y = 0;
    if(!region.empty())
    {
        min_x = region[0];
        max_x = region[1];
        min_y = region[2];
        max_y = region[3];
    }

    unsigned long long uniq_cell_id;
    uint32_t index = 0,gid = 0;
    std::unordered_map<unsigned long long, uint32_t> hash_map;

    Gene * gene = getGene();
    Expression * expression = getExpression();

    if(genelist.empty()&& (!region.empty()))
    {
        unordered_map<string, vector<Expression>> hash_exp;
        ThreadPool tpool(n_thread_);
        for(unsigned int gene_id = 0; gene_id < gene_num_; gene_id++)
        {
            getdataTask *ptask = new getdataTask(gene_id, gene, expression, hash_exp);
            ptask->setRange(min_x, min_y, max_x, max_y);
            tpool.addTask(ptask);
        }
        tpool.waitTaskDone();

        auto itor = hash_exp.begin();
        for(;itor != hash_exp.end();itor++)
        {
            vec_gene.emplace_back(itor->first);
            vector<Expression> &tmp = itor->second;
            for(Expression &expData : tmp)
            {
                uniq_cell_id = expData.x;
                uniq_cell_id = (uniq_cell_id << 32) | expData.y;

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
                count.push_back(expData.count);
                gene_ind.push_back(gid);
            }
            gid++;
        }
    }
    else if(region.empty() && (!genelist.empty()))
    {
        set<string> gset;
        for(string &str : genelist)
        {
            gset.insert(str);
        }
        for(unsigned int gene_id = 0; gene_id < gene_num_; gene_id++)
        {
            string str(gene[gene_id].gene);
            if(gset.find(str) != gset.end())
            {
                vec_gene.emplace_back(str);
                unsigned int end = gene[gene_id].offset + gene[gene_id].count;
                for(unsigned int i = gene[gene_id].offset; i < end; i++)
                {
                    uniq_cell_id = expression[i].x;
                    uniq_cell_id = (uniq_cell_id << 32) | expression[i].y;

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
                    count.push_back(expression[i].count);
                    gene_ind.push_back(gid); 
                }
                gid++;
            }
        }
    }
    else if((!region.empty()) && (!genelist.empty()))
    {
        set<string> gset;
        for(string &str : genelist)
        {
            gset.insert(str);
        }
        for(unsigned int gene_id = 0; gene_id < gene_num_; gene_id++)
        {
            string str(gene[gene_id].gene);
            if(gset.find(str) != gset.end())
            {
                vec_gene.emplace_back(str);
                unsigned int end = gene[gene_id].offset + gene[gene_id].count;
                for(unsigned int i = gene[gene_id].offset; i < end; i++)
                {
                    Expression &exp = expression[i];
                    if(exp.x < min_x || exp.x >= max_x || exp.y < min_y || exp.y >= max_y){
                        continue;
                    }

                    uniq_cell_id = expression[i].x;
                    uniq_cell_id = (uniq_cell_id << 32) | expression[i].y;

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
                    count.push_back(expression[i].count);
                    gene_ind.push_back(gid);                
                }
                gid++;
            }
        }
    }
    else 
    {
        for(unsigned int gene_id = 0; gene_id < gene_num_; gene_id++)
        {
            vec_gene.emplace_back(gene[gene_id].gene);
            unsigned int end = gene[gene_id].offset + gene[gene_id].count;
            for(unsigned int i = gene[gene_id].offset; i < end; i++)
            {
                uniq_cell_id = expression[i].x;
                uniq_cell_id = (uniq_cell_id << 32) | expression[i].y;

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
                count.push_back(expression[i].count);
                gene_ind.push_back(gene_id);
            }
        }
    }
}

void BgefReader::getfiltereddata_exon(vector<int> &region, vector<string> &genelist,
                                 vector<string> &vec_gene, vector<unsigned long long> &uniq_cells,
                                 vector<unsigned int> &cell_ind, vector<unsigned int> &gene_ind, 
                                 vector<unsigned int> &count, vector<unsigned int> &exon)
{
    int min_x = 0, max_x = 0, min_y = 0, max_y = 0;
    if(!region.empty())
    {
        min_x = region[0];
        max_x = region[1];
        min_y = region[2];
        max_y = region[3];
    }

    unsigned long long uniq_cell_id;
    uint32_t index = 0,gid = 0;
    std::unordered_map<unsigned long long, uint32_t> hash_map;

    Gene * gene = getGene();
    Expression * expression = getExpression();

    if(genelist.empty()&& (!region.empty()))
    {
        unordered_map<string, vector<Expression>> hash_exp;
        ThreadPool tpool(n_thread_);
        for(unsigned int gene_id = 0; gene_id < gene_num_; gene_id++)
        {
            getdataTask *ptask = new getdataTask(gene_id, gene, expression, hash_exp);
            ptask->setRange(min_x, min_y, max_x, max_y);
            tpool.addTask(ptask);
        }
        tpool.waitTaskDone();

        auto itor = hash_exp.begin();
        for(;itor != hash_exp.end();itor++)
        {
            vec_gene.emplace_back(itor->first);
            vector<Expression> &tmp = itor->second;
            for(Expression &expData : tmp)
            {
                uniq_cell_id = expData.x;
                uniq_cell_id = (uniq_cell_id << 32) | expData.y;

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
                exon.push_back(expData.exon);
                count.push_back(expData.count);
                gene_ind.push_back(gid);
            }
            gid++;
        }
    }
    else if(region.empty() && (!genelist.empty()))
    {
        set<string> gset;
        for(string &str : genelist)
        {
            gset.insert(str);
        }
        for(unsigned int gene_id = 0; gene_id < gene_num_; gene_id++)
        {
            string str(gene[gene_id].gene);
            if(gset.find(str) != gset.end())
            {
                vec_gene.emplace_back(str);
                unsigned int end = gene[gene_id].offset + gene[gene_id].count;
                for(unsigned int i = gene[gene_id].offset; i < end; i++)
                {
                    uniq_cell_id = expression[i].x;
                    uniq_cell_id = (uniq_cell_id << 32) | expression[i].y;

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
                    exon.push_back(expression[i].exon);
                    count.push_back(expression[i].count);
                    gene_ind.push_back(gid); 
                }
                gid++;
            }
        }
    }
    else if((!region.empty()) && (!genelist.empty()))
    {
        set<string> gset;
        for(string &str : genelist)
        {
            gset.insert(str);
        }
        for(unsigned int gene_id = 0; gene_id < gene_num_; gene_id++)
        {
            string str(gene[gene_id].gene);
            if(gset.find(str) != gset.end())
            {
                vec_gene.emplace_back(str);
                unsigned int end = gene[gene_id].offset + gene[gene_id].count;
                for(unsigned int i = gene[gene_id].offset; i < end; i++)
                {
                    Expression &exp = expression[i];
                    if(exp.x < min_x || exp.x >= max_x || exp.y < min_y || exp.y >= max_y){
                        continue;
                    }

                    uniq_cell_id = expression[i].x;
                    uniq_cell_id = (uniq_cell_id << 32) | expression[i].y;

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
                    exon.push_back(expression[i].exon);
                    count.push_back(expression[i].count);
                    gene_ind.push_back(gid);                
                }
                gid++;
            }
        }
    }
    else 
    {
        for(unsigned int gene_id = 0; gene_id < gene_num_; gene_id++)
        {
            vec_gene.emplace_back(gene[gene_id].gene);
            unsigned int end = gene[gene_id].offset + gene[gene_id].count;
            for(unsigned int i = gene[gene_id].offset; i < end; i++)
            {
                uniq_cell_id = expression[i].x;
                uniq_cell_id = (uniq_cell_id << 32) | expression[i].y;

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
                exon.push_back(expression[i].exon);
                count.push_back(expression[i].count);
                gene_ind.push_back(gene_id);
            }
        }
    }
}


void BgefReader::getOffset(int *data)
{
    if(data)
    {
        ExpressionAttr & expression_attr = getExpressionAttr();
        data[0] = expression_attr.min_x;
        data[1] = expression_attr.min_y;
    }
}

void BgefReader::getExpAttr(int *data)
{
    if(data)
    {
        ExpressionAttr & expression_attr = getExpressionAttr();
        data[0] = expression_attr.min_x;
        data[1] = expression_attr.min_y;
        data[2] = expression_attr.max_x;
        data[3] = expression_attr.max_y;
        data[4] = expression_attr.max_exp;
        data[5] = expression_attr.resolution;
    }
}

////////////////////////////
unsigned int *BgefReader::getGeneExon()
{
    if(m_bexon)
    {
        if(m_exonPtr) return m_exonPtr;

        hsize_t dims[1];
        hid_t sid = H5Dget_space(m_exon_did);
        H5Sget_simple_extent_dims(sid, dims, nullptr);
        assert(dims[0] == expression_num_);
        m_exonPtr = new unsigned int[dims[0]];
        H5Dread(m_exon_did, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_exonPtr);
        H5Sclose(sid);
        return m_exonPtr;
    }
    return nullptr;
}

unsigned int BgefReader::getGeneExonAttr() {
    if (!isExonExist()) return 0;
    hid_t attr;
    attr = H5Aopen(m_exon_did, "maxExon", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &max_exon_);

    H5Aclose(attr);
    return max_exon_;
}

Expression *BgefReader::getExpression_abs() {
    if(expressions_ != nullptr)
        return expressions_;

    ExpressionAttr & expression_attr = getExpressionAttr();

    hid_t memtype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(Expression));
    H5Tinsert(memtype, "x", HOFFSET(Expression, x), H5T_NATIVE_INT);
    H5Tinsert(memtype, "y", HOFFSET(Expression, y), H5T_NATIVE_INT);
    H5Tinsert(memtype, "count", HOFFSET(Expression, count), H5T_NATIVE_UINT);

    expressions_ = (Expression *) malloc(expression_num_ * sizeof(Expression));
    H5Dread(exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, expressions_);

    H5Tclose(memtype);

    getGeneExon();
    if(m_exonPtr)
    {
        for(unsigned long long i=0;i<expression_num_;i++)
        {
            expressions_[i].x += expression_attr.min_x;
            expressions_[i].y += expression_attr.max_y;
            expressions_[i].exon = m_exonPtr[i];
        }
    }
    else
    {
        for(unsigned long long i=0;i<expression_num_;i++)
        {
            expressions_[i].x += expression_attr.min_x;
            expressions_[i].y += expression_attr.max_y;
        }
    }
    return expressions_;
}