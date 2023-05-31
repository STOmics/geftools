
#include "cgef_reader.h"

CgefReader::CgefReader(const string &filename, bool verbose) {
    str32_type_ = H5Tcopy(H5T_C_S1);
    H5Tset_size(str32_type_, 32);
    verbose_ = verbose;

    file_id_ = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    group_id_ = H5Gopen(file_id_, "/cellBin", H5P_DEFAULT);
    getAttr();
    cell_dataset_id_ = openCellDataset(group_id_);
    cell_exp_dataset_id_ = openCellExpDataset(group_id_);
    isOldCellExpVersion = isOlderCellExpDataVersion(file_id_);
    
    gene_dataset_id_ = openGeneDataset(group_id_);
    gene_exp_dataset_id_ = openGeneExpDataset(group_id_);
    gene_exp_dataspace_id_ = H5Dget_space(gene_exp_dataset_id_);

    hsize_t dims[1];
    cell_exp_dataspace_id_ = H5Dget_space(cell_exp_dataset_id_);
    H5Sget_simple_extent_dims(cell_exp_dataspace_id_, dims, nullptr);
    expression_num_ = dims[0];
    expression_num_current_ = dims[0];

    cell_dataspace_id_ = H5Dget_space(cell_dataset_id_);
    H5Sget_simple_extent_dims(cell_dataspace_id_, dims, nullptr);
    cell_num_ = dims[0];
    cell_num_current_ = dims[0];

    gene_array_ = loadGene();

    char dname[128] = {0};
    sprintf(dname, "/cellBin/cellExon");
    if(H5Lexists(file_id_, dname, H5P_DEFAULT)>0)
    {
        m_bexon = true;
    }
}

CgefReader::~CgefReader() {
    closeH5();
}

void CgefReader::closeH5()
{
    if(gene_array_ == nullptr) return;
    H5Tclose(str32_type_);
    H5Dclose(cell_dataset_id_);
    H5Dclose(gene_dataset_id_);
    H5Dclose(cell_exp_dataset_id_);
    H5Dclose(gene_exp_dataset_id_);
    H5Sclose(cell_dataspace_id_);
    H5Sclose(cell_exp_dataspace_id_);
    H5Sclose(gene_exp_dataspace_id_);
    H5Gclose(group_id_);
    H5Fclose(file_id_);
    free(gene_array_);
    gene_array_ = nullptr;
    if (cell_array_ != nullptr)
        free(cell_array_);
    if (cell_array_current_ != nullptr) free(cell_array_current_);
    if (cell_id_array_current_ != nullptr) free(cell_id_array_current_);
    if (cell_id_to_index_ != nullptr) free(cell_id_to_index_);

    // if(m_borderdata_currentPtr!=nullptr) free(m_borderdata_currentPtr);
    // if(m_borderdataPtr!=nullptr) free(m_borderdataPtr);
}

void CgefReader::getAttr()
{
    if(m_ver == 0)
    {
        hid_t attr = H5Aopen(file_id_, "version", H5P_DEFAULT);
        H5Aread(attr, H5T_NATIVE_UINT32, &m_ver);

        attr = H5Aopen(file_id_, "resolution", H5P_DEFAULT);
        H5Aread(attr, H5T_NATIVE_UINT32, &m_resolution);

        attr = H5Aopen(file_id_, "offsetX", H5P_DEFAULT);
        H5Aread(attr, H5T_NATIVE_INT32, &offsetX);

        attr = H5Aopen(file_id_, "offsetY", H5P_DEFAULT);
        H5Aread(attr, H5T_NATIVE_INT32, &offsetY);

        attr = H5Aopen(file_id_, "geftool_ver", H5P_DEFAULT);
        H5Aread(attr, H5T_NATIVE_UINT32, m_ver_tool);

        H5Aclose(attr);
    }
}

hid_t CgefReader::openCellDataset(hid_t group_id) {
    cell_dataset_id_ = H5Dopen(group_id, "cell", H5P_DEFAULT);
    if (cell_dataset_id_ < 0) {
        cerr << "failed open dataset: cell" << endl;
        reportErrorCode2File(errorCode::E_MISSINGFILEINFO, "failed open dataset: cell");
        exit(3);
    }

    hid_t s1_tid = H5Dget_type(cell_dataset_id_);
    int nmemb = H5Tget_nmembers(s1_tid);

    if(nmemb < 9){
        cerr << "Please use geftools(>=0.6) to regenerate this cgef file." << endl;
        reportErrorCode2File(errorCode::E_LOWVERSION, 
                            "Please use geftools(>=0.6) to regenerate this cgef file.");
        exit(2);
    }

    if(H5Aexists(cell_dataset_id_, "blockIndex"))
    {
        hsize_t dims_attr[1];
        hid_t attr, attr_dataspace;
        attr = H5Aopen(cell_dataset_id_, "blockIndex", H5P_DEFAULT);
        attr_dataspace = H5Aget_space(attr);
        H5Sget_simple_extent_dims(attr_dataspace, dims_attr, nullptr);

        block_index_ = static_cast<unsigned int *>(malloc(dims_attr[0] * sizeof(unsigned int)));

        H5Aread(attr, H5T_NATIVE_UINT32, block_index_);

        attr = H5Aopen(cell_dataset_id_, "blockSize", H5P_DEFAULT);
        H5Aread(attr, H5T_NATIVE_UINT32, block_size_);

        H5Aclose(attr);
        H5Sclose(attr_dataspace);
    }
    else
    {
        hid_t d_id = 0;
        if(H5Lexists(group_id, "blockIndex", H5P_DEFAULT)>0)
        {
            d_id = H5Dopen(group_id, "blockIndex", H5P_DEFAULT);
        }
        else if(H5Lexists(group_id, "blkidx", H5P_DEFAULT)>0)
        {
            d_id = H5Dopen(group_id, "blkidx", H5P_DEFAULT);
        }
        
        hsize_t dims[1];
        hid_t s_id = H5Dget_space(d_id);
        H5Sget_simple_extent_dims(s_id, dims, nullptr);
        block_index_ = static_cast<unsigned int *>(calloc(dims[0], sizeof(unsigned int)));
        H5Dread(d_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, block_index_);

        H5Sclose(s_id);
        H5Dclose(d_id);

        d_id = H5Dopen(group_id, "blockSize", H5P_DEFAULT);
        H5Dread(d_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, block_size_);
        H5Dclose(d_id);
    }

    return cell_dataset_id_;
}

hid_t CgefReader::openGeneDataset(hid_t group_id) {
    hsize_t dims[1];
    gene_dataset_id_ = H5Dopen(group_id, "gene", H5P_DEFAULT);
    if (gene_dataset_id_ < 0) {
        cerr << "failed open dataset: gene" << endl;
        return gene_dataset_id_;
    }
    hid_t gene_dataspace_id = H5Dget_space(gene_dataset_id_);
    H5Sget_simple_extent_dims(gene_dataspace_id, dims, nullptr);
    gene_num_ = dims[0];
    gene_num_current_ = dims[0];
    H5Sclose(gene_dataspace_id);
    return gene_dataset_id_;
}

hid_t CgefReader::openCellExpDataset(hid_t group_id) {
    cell_exp_dataset_id_ = H5Dopen(group_id, "cellExp", H5P_DEFAULT);
    if (cell_exp_dataset_id_ < 0) {
        cerr << "failed open dataset: cellExp" << endl;
        reportErrorCode2File(errorCode::E_MISSINGFILEINFO, "failed open dataset: cellExp");
        exit(3);
    }
    return cell_exp_dataset_id_;
}

hid_t CgefReader::openGeneExpDataset(hid_t group_id) {
    gene_exp_dataset_id_ = H5Dopen(group_id, "geneExp", H5P_DEFAULT);
    if (gene_exp_dataset_id_ < 0) {
        cerr << "failed open dataset: geneExp" << endl;
        return gene_exp_dataset_id_;
    }
    return gene_exp_dataset_id_;
}

unsigned int CgefReader::getCellNum() const {
    return cell_num_current_;
}

unsigned int CgefReader::getGeneNum() const {
    return gene_num_current_;
}

unsigned long long int CgefReader::getExpressionNum() const {
    return expression_num_current_;
}

GeneData *CgefReader::loadGene(bool reload) {
    unsigned long cprev = clock();
    if (gene_array_ != nullptr) {
        if (reload) free(gene_array_);
        else return gene_array_;
    }

    hid_t memtype = getMemtypeOfGeneData();
    gene_array_ = (GeneData *) malloc(gene_num_ * sizeof(GeneData));
    H5Dread(gene_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, gene_array_);

    for (unsigned int i = 0; i < gene_num_; i++) {
        genename_to_id_[gene_array_[i].gene_name] = i;
    }

    gene_id_to_index_ = (int *) malloc(gene_num_ * sizeof(int));
    iota(gene_id_to_index_, gene_id_to_index_+gene_num_, 0);

    if (verbose_) printCpuTime(cprev, "loadGene");
    return gene_array_;
}

CellData *CgefReader::loadCell(bool reload) {
    unsigned long cprev = clock();
    if (cell_array_ != nullptr) {
        if (reload) free(cell_array_);
        else return cell_array_;
    }

    hid_t memtype = getMemtypeOfCellData();
    cell_array_ = (CellData *) malloc(cell_num_ * sizeof(CellData));
    H5Dread(cell_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell_array_);
    if (verbose_) printCpuTime(cprev, "getCell");
    return cell_array_;
}

void CgefReader::getGeneNameList(vector<string> &gene_list) {
    for(unsigned int gene_id = 0; gene_id < gene_num_; gene_id++) {
        if(gene_id_to_index_[gene_id]<0) continue;
        gene_list.emplace_back(gene_array_[gene_id].gene_name);
    }
}

void CgefReader::getGeneNames(char *gene_list) {
    int i = 0;
    for(unsigned int gene_id = 0; gene_id < gene_num_; gene_id++) {
        if(gene_id_to_index_[gene_id]<0) continue;
        memcpy(&gene_list[i * 32], gene_array_[gene_id].gene_name, 32);
        i++;
    }
}

void CgefReader::getCellNameList(unsigned long long int *cell_pos_list) {
    if (restrict_region_) {
        for (unsigned int i = 0; i < cell_num_current_; i++) {
            cell_pos_list[i] = cell_array_current_[i].x;
            cell_pos_list[i] = cell_pos_list[i] << 32 | cell_array_current_[i].y;
        }
    } else {
        CellData *cells = loadCell();
        for (unsigned int i = 0; i < cell_num_; i++) {
            cell_pos_list[i] = cells[i].x;
            cell_pos_list[i] = cell_pos_list[i] << 32 | cells[i].y;
        }
    }
}

//indptr length = gene_num_ + 1
int CgefReader::getSparseMatrixIndices(unsigned int *indices, unsigned int *indptr, unsigned int *count,
                                       const char *order) {
    if (order[0] == 'g') {
        if (restrict_region_ || gene_num_current_ < gene_num_) {
            unsigned int n = 0, rows = 0;
            indptr[0] = 0;

            auto *gene_exp_data = static_cast<GeneExpData *>(malloc(expression_num_current_ * sizeof(GeneExpData)));

            for(unsigned int gene_id = 0; gene_id < gene_num_; gene_id++) {
                if(gene_id_to_index_[gene_id]<0) continue;
                GeneData gene_data = gene_array_[gene_id];
                if(gene_data.cell_count ==0 ) {
                    indptr[rows+1] = indptr[rows];
                    rows ++;
                    continue;
                }

                selectGeneExp(gene_data.offset, gene_data.cell_count, gene_exp_data);

                unsigned int c_count = 0;
                for (unsigned int j = 0; j < gene_data.cell_count; j++) {
                    unsigned int cid = gene_exp_data[j].cell_id;
                    if(restrict_region_){
                        if(!isInRegion(cid)) continue;
                        indices[n] = cell_id_to_index_[cid-start_cell_id];
                    }else{
                        indices[n] = cid;
                    }

                    count[n] = gene_exp_data[j].count;
                    n++;
                    c_count++;
                }
                indptr[rows+1] = indptr[rows] + c_count;
                rows ++;
            }
            assert(rows == gene_num_current_);
            assert(n == expression_num_current_);
            free(gene_exp_data);
        } else {
            hid_t memtype;
            memtype = H5Tcreate(H5T_COMPOUND, sizeof(unsigned int));
            H5Tinsert(memtype, "count", 0, H5T_NATIVE_USHORT);
            H5Dread(gene_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, count);
            memtype = H5Tcreate(H5T_COMPOUND, sizeof(unsigned int));
            H5Tinsert(memtype, "cellID", 0, H5T_NATIVE_UINT);
            H5Dread(gene_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, indices);
            for (unsigned int i = 0; i < gene_num_; i++) {
                indptr[i] = gene_array_[i].offset;
            }
            indptr[gene_num_] = gene_array_[gene_num_ - 1].offset + gene_array_[gene_num_ - 1].cell_count;
            H5Tclose(memtype);
        }

    } else if (order[0] == 'c') {
        if (restrict_region_ || gene_num_current_ < gene_num_) {
            unsigned int n = 0, rows = 0;
            indptr[0] = 0;

            if (isOldCellExpVersion) {
                auto *cell_exp_data = static_cast<olderCellExpData *>(malloc(expression_num_current_ * sizeof(olderCellExpData)));
                CellData * cell_datas = getCell();
                for (unsigned int i = 0; i < cell_num_current_; i++) {
                    CellData cell = cell_datas[i];
                    selectOlderCellExp(cell.offset, cell.gene_count, cell_exp_data);

                    unsigned short gene_id, g_count = 0;
                    for (unsigned int j = 0; j < cell.gene_count; j++) {
                        gene_id = cell_exp_data[j].gene_id;
                        if(gene_id_to_index_[gene_id]<0) continue;
                        indices[n] = gene_id_to_index_[gene_id];
                        count[n] = cell_exp_data[j].count;
                        n++;
                        g_count++;
                    }
                    indptr[rows+1] = indptr[rows] + g_count;
                    rows++;
                }
                assert(n == expression_num_current_);
                free(cell_exp_data);
            } else {
                auto *cell_exp_data = static_cast<CellExpData *>(malloc(expression_num_current_ * sizeof(CellExpData)));
                CellData * cell_datas = getCell();
                for (unsigned int i = 0; i < cell_num_current_; i++) {
                    CellData cell = cell_datas[i];
                    selectCellExp(cell.offset, cell.gene_count, cell_exp_data);

                    unsigned short gene_id, g_count = 0;
                    for (unsigned int j = 0; j < cell.gene_count; j++) {
                        gene_id = cell_exp_data[j].gene_id;
                        if(gene_id_to_index_[gene_id]<0) continue;
                        indices[n] = gene_id_to_index_[gene_id];
                        count[n] = cell_exp_data[j].count;
                        n++;
                        g_count++;
                    }
                    indptr[rows+1] = indptr[rows] + g_count;
                    rows++;
                }
                assert(n == expression_num_current_);
                free(cell_exp_data);
            }
        } else {
            // todo wanruiwen
            hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(unsigned int));
            H5Tinsert(memtype, "count", 0, H5T_NATIVE_USHORT);
            H5Dread(cell_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, count);
            memtype = H5Tcreate(H5T_COMPOUND, sizeof(unsigned int));
            H5Tinsert(memtype, "geneID", 0, H5T_NATIVE_USHORT);
            H5Dread(cell_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, indices);

            CellData *cell_data = loadCell();
            indptr[0] = 0;
            for (unsigned int i = 1; i < cell_num_; i++) {
                indptr[i] = cell_data[i].offset;
            }
            indptr[cell_num_] = cell_data[cell_num_ - 1].offset + cell_data[cell_num_ - 1].gene_count;
            H5Tclose(memtype);
        }
    } else {
        return -1;
    }

    return 0;
}

int CgefReader::getSparseMatrixIndices2(unsigned int *cell_ind, unsigned int *gene_ind, unsigned int *count) {
    hid_t memtype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(unsigned int));
    H5Tinsert(memtype, "count", 0, H5T_NATIVE_USHORT);
    H5Dread(gene_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, count);
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(unsigned int));
    H5Tinsert(memtype, "cellID", 0, H5T_NATIVE_UINT);
    H5Dread(gene_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell_ind);

    unsigned int n = 0;
    for (unsigned int i = 0; i < gene_num_; i++) {
        unsigned int cell_count = gene_array_[i].cell_count;
        for (unsigned int j = 0; j < cell_count; j++) {
            gene_ind[n++] = i;
        }
    }

    H5Tclose(memtype);
    return 0;
}

void CgefReader::getCellIdAndCount(unsigned int *cell_id, unsigned short *count) const{
    hid_t memtype = getMemtypeOfGeneExpData();
    GeneExpData* gene_exp_data;
    gene_exp_data = (GeneExpData*)malloc(expression_num_ * sizeof(GeneExpData));
    H5Dread(gene_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, gene_exp_data);
    for(unsigned int i = 0; i < expression_num_; i++){
        cell_id[i] = gene_exp_data->cell_id;
        count[i] = gene_exp_data->count;
    }
    free(gene_exp_data);
}

void CgefReader::getGeneIdAndCount(unsigned int *gene_id, unsigned short *count) const{
    if(isOldCellExpVersion) {
        hid_t memtype = getMemtypeOfOlderCellExpData();
        olderCellExpData* cell_exp_data;
        cell_exp_data = (olderCellExpData*)malloc(expression_num_ * sizeof(olderCellExpData));
        H5Dread(cell_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell_exp_data);
        for (unsigned long long i = 0; i < expression_num_; i++) {
            gene_id[i] = cell_exp_data->gene_id;
            count[i] = cell_exp_data->count;
        }
        free(cell_exp_data);
    } else {
        hid_t memtype = getMemtypeOfCellExpData();
        CellExpData* cell_exp_data;
        cell_exp_data = (CellExpData*)malloc(expression_num_ * sizeof(CellExpData));
        H5Dread(cell_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell_exp_data);
        for (unsigned long long i = 0; i < expression_num_; i++) {
            gene_id[i] = cell_exp_data->gene_id;
            count[i] = cell_exp_data->count;
        }
        free(cell_exp_data);
    }
}

unsigned int CgefReader::getExpressionCountByGene(string &gene_name, GeneExpData *expressions) {
    int gene_id = getGeneId(gene_name);
    if (gene_id < 0) {
        cerr << "Gene ID < 0 : " << gene_id << endl;
        reportErrorCode2File(errorCode::E_FILEDATAERROR, "Gene ID < 0 : ");
        exit(2);
    }
    return getExpressionCountByGeneId(gene_id, expressions);
}

unsigned int CgefReader::getExpressionCountByGeneId(unsigned int gene_id, GeneExpData *expressions) {
    unsigned int cell_count = gene_array_[gene_id].cell_count;
    selectGeneExp(gene_array_[gene_id].offset, cell_count, expressions);
    if(restrict_region_){
        unsigned int offset_cell_id, new_cell_count = 0;
        for(unsigned int i = 0; i < cell_count; i++){
            if(isInRegion(expressions[i].cell_id)){
                memmove(&expressions[new_cell_count], &expressions[i], sizeof(GeneExpData));
                new_cell_count ++;
            }
        }
        memset(&expressions[new_cell_count], 0, sizeof(GeneExpData));  // TODO
        return new_cell_count;
    }
    return cell_count;
}

GeneData CgefReader::getGeneDataByGeneId(unsigned int gene_id) {
    return gene_array_[gene_id];
}

CellData CgefReader::getCell(unsigned int cell_id) const {
    CellData cell_data{};
    selectCells(cell_id, 1, &cell_data);
    return cell_data;
}

CellData *CgefReader::getCell() {
    if (cell_array_current_ != nullptr) {
        return cell_array_current_;
    }
    return loadCell();
}

string CgefReader::getGeneName(unsigned int gene_id) {
    return gene_array_[gene_id].gene_name;
}

GeneData CgefReader::getGene(unsigned int gene_id) const {
    return gene_array_[gene_id];
}

GeneData *CgefReader::getGene() {
    if (gene_array_current_ != nullptr) {
        return gene_array_current_;
    } else if (gene_num_current_ < gene_num_) {
        gene_array_current_ = static_cast<GeneData *>(malloc(gene_num_current_ * sizeof(GeneData)));
        int i = 0;
        for (unsigned int gene_id = 0; gene_id < gene_num_; gene_id++) {
            if(gene_id_to_index_[gene_id]<0) continue;
            memcpy(&gene_array_current_[i], &gene_array_[gene_id], sizeof(GeneData));
            i++;
        }
        assert(i == gene_num_current_);
        return gene_array_current_;
    }
    return gene_array_;
}

int CgefReader::getGeneId(string &gene_name) {
    auto iter = genename_to_id_.find(gene_name);
    if (iter != genename_to_id_.end()) {
        return iter->second;
    }
    return -1;
}

// unsigned int CgefReader::toGem(string &filename, const vector<string> &gene_name_list, bool force_genes, bool exclude) {
//     unsigned long cprev = clock();
//     unsigned short gene_ids[gene_num_];
//     unsigned short n = 0;
//     unsigned int cell_num = 0;
//     unordered_map<string, bool> genename_map;
//     for (const auto &gene_name: gene_name_list) {
//         genename_map[gene_name] = true;
//     }
//     if (exclude) {
//         for (unsigned int i = 0; i < gene_num_; i++) {
//             if (genename_map.find(gene_array_[i].gene_name) == genename_map.end()) {
//                 gene_ids[n++] = i;
//             }
//         }
//     } else {
//         for (const auto &gene_name: gene_name_list) {
//             if (genename_to_id_.find(gene_name) == genename_to_id_.end()) {
//                 cerr << "Gene ( " << gene_name << " ) does not exist." << endl;
//                 if (!force_genes) exit(2);
//                 continue;
//             }
//             gene_ids[n++] = genename_to_id_[gene_name];
//         }
//     }

//     ofstream fout;

//     bool to_file = true;
//     if (filename != "stdout") {
//         fout.open(filename);
//         if (!fout.is_open()) cerr << "Fail to open file : " << filename << endl;
//         fout << "#geneName\tx\ty\tcount\tcellID" << endl;
//     } else {
//         to_file = false;
//         cout << "#geneName\tx\ty\tcount\tcellID" << endl;
//     }

//     if (restrict_region_) {
//         if (verbose_) cerr << "toGem restrict_region_ true" << endl;
//         auto *cell_exp_data = static_cast<CellExpData *>(malloc(gene_num_ * sizeof(CellExpData)));
//         for (unsigned int i = 0; i < cell_num_current_; i++) {
//             unsigned int cell_id = cell_id_array_current_[i];
//             CellData cell = cell_array_current_[i];
//             selectCellExp(cell.offset, cell.gene_count, cell_exp_data);

//             for (unsigned int j = 0; j < cell.gene_count; j++) {
//                 string gene_name = gene_array_[cell_exp_data[j].gene_id].gene_name;

//                 if (!gene_name_list.empty() && (
//                         (exclude && genename_map.find(gene_name) != genename_map.end())
//                         || (!exclude && genename_map.find(gene_name) == genename_map.end())
//                 ))
//                     continue;

//                 cell_num++;
//                 if (to_file) {
//                     fout << gene_name << "\t" << cell.x << "\t" << cell.y
//                          << "\t" << cell_exp_data[j].count << "\t" << cell_id << endl;
//                 } else {
//                     cout << gene_name << "\t" << cell.x << "\t" << cell.y
//                          << "\t" << cell_exp_data[j].count << "\t" << cell_id << endl;
//                 }
//             }
//         }

//         free(cell_exp_data);
//     } else {
//         unsigned int max_malloc = 1000;
//         auto *expression = static_cast<GeneExpData *>(malloc(max_malloc * sizeof(GeneExpData)));

//         for (unsigned short i = 0; i < n; i++) {
//             unsigned short gene_id = gene_ids[i];
//             GeneData gene_data = getGeneDataByGeneId(gene_id);
//             unsigned int cell_count = gene_data.cell_count;

//             cell_num += cell_count;

//             if (cell_count > max_malloc) {
//                 expression = static_cast<GeneExpData *>(realloc(expression, cell_count));
//                 max_malloc = cell_count + 1000;
//             }

//             selectGeneExp(gene_data.offset, cell_count, expression);

//             for (unsigned int j = 0; j < cell_count; j++) {
//                 CellData cell_data = getCell(expression[j].cell_id);
//                 if (to_file) {
//                     fout << gene_data.gene_name << "\t" << cell_data.x << "\t" << cell_data.y
//                          << "\t" << expression[j].count << "\t" << expression[j].cell_id << endl;
//                 } else {
//                     cout << gene_data.gene_name << "\t" << cell_data.x << "\t" << cell_data.y
//                          << "\t" << expression[j].count << "\t" << expression[j].cell_id << endl;
//                 }
//             }
//         }
//         free(expression);
//     }

//     if (verbose_) printCpuTime(cprev, "toGem");

//     if (to_file) fout.close();
//     return cell_num;
// }

bool CgefReader::isVerbose() const {
    return verbose_;
}

void CgefReader::setVerbose(bool verbose) {
    verbose_ = verbose;
}

void CgefReader::selectCells(unsigned int offset,
                             unsigned int cell_count,
                             CellData *cell) const {
    hsize_t start[1] = {offset},
            count[1] = {cell_count},
            offset_out[1] = {0};

    hid_t memtype;
    memtype = getMemtypeOfCellData();

    // Define memory dataspace.
    hid_t memspace = H5Screate_simple(1, count, nullptr);
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, nullptr, count, nullptr);

    H5Sselect_hyperslab(cell_dataspace_id_, H5S_SELECT_SET, start, nullptr, count, nullptr);
    H5Dread(cell_dataset_id_, memtype, memspace, cell_dataspace_id_, H5P_DEFAULT, cell);
}

void CgefReader::selectCellExp(unsigned int offset,
                               unsigned int gene_count,
                               CellExpData *cell_exp_data) const {
    hsize_t start[1] = {offset},
            count[1] = {gene_count},
            offset_out[1] = {0};

    hid_t memtype;
    memtype = getMemtypeOfCellExpData();

    // Define memory dataspace.
    hid_t memspace = H5Screate_simple(1, count, nullptr);
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, nullptr, count, nullptr);

    H5Sselect_hyperslab(cell_exp_dataspace_id_, H5S_SELECT_SET, start, nullptr, count, nullptr);
    H5Dread(cell_exp_dataset_id_, memtype, memspace, cell_exp_dataspace_id_, H5P_DEFAULT, cell_exp_data);
}

void CgefReader::selectOlderCellExp(unsigned int offset,
                               unsigned int gene_count,
                               olderCellExpData *cell_exp_data) const {
    hsize_t start[1] = {offset},
            count[1] = {gene_count},
            offset_out[1] = {0};

    hid_t memtype = getMemtypeOfOlderCellExpData();

    // Define memory dataspace.
    hid_t memspace = H5Screate_simple(1, count, nullptr);
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, nullptr, count, nullptr);

    H5Sselect_hyperslab(cell_exp_dataspace_id_, H5S_SELECT_SET, start, nullptr, count, nullptr);
    H5Dread(cell_exp_dataset_id_, memtype, memspace, cell_exp_dataspace_id_, H5P_DEFAULT, cell_exp_data);
}

void CgefReader::selectGeneExp(unsigned int offset,
                               unsigned int cell_count,
                               GeneExpData *gene_exp_data) const {
    hsize_t start[1] = {offset},
            count[1] = {cell_count},
            offset_out[1] = {0};

    hid_t memtype = getMemtypeOfGeneExpData();

    hid_t memspace = H5Screate_simple(1, count, nullptr);
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, nullptr, count, nullptr);

    H5Sselect_hyperslab(gene_exp_dataspace_id_, H5S_SELECT_SET, start, nullptr, count, nullptr);
    H5Dread(gene_exp_dataset_id_, memtype, memspace, gene_exp_dataspace_id_, H5P_DEFAULT, gene_exp_data);
}

bool CgefReader::isRestrictRegion() const {
    return restrict_region_;
}

bool CgefReader::isRestrictGene() const {
    return restrict_gene_;
}

void CgefReader::restrictRegion(unsigned int min_x, unsigned int max_x, unsigned int min_y, unsigned int max_y) {
    unsigned long cprev = clock();
    if (restrict_gene_ || restrict_region_) {
        cerr << "Please call freeRestriction first, or call restrictRegion function before restrictGene." << endl;
        reportErrorCode2File(errorCode::E_STEPERROR, 
        "Please call freeRestriction first, or call restrictRegion function before restrictGene.");
        exit(2);
    }

    restrict_region_ = true;
    unsigned int x_block_num = block_size_[2];
    unsigned int y_block_num = block_size_[3];
    unsigned int min_block_x = min_x / block_size_[0];
    unsigned int max_block_x = max_x / block_size_[0];
    unsigned int min_block_y = min_y / block_size_[1];
    unsigned int max_block_y = max_y / block_size_[1];

    max_block_x = max_block_x > x_block_num ? x_block_num : max_block_x;
    max_block_y = max_block_y > y_block_num ? y_block_num : max_block_y;

    unsigned int block_id_y, offset, cell_num_tmp = 0, cell_id_tmp;
    for (unsigned int y = min_block_y; y <= max_block_y; y++) {
        block_id_y = y * x_block_num;
        cell_num_tmp += block_index_[max_block_x + block_id_y + 1] - block_index_[min_block_x + block_id_y];
    }

    start_cell_id = block_index_[min_block_x+min_block_y*x_block_num];
    end_cell_id = block_index_[max_block_x+max_block_y*x_block_num + 1];

    // support rerun this method
    cell_num_current_ = 0;
    expression_num_current_ = 0;
    cell_array_current_ = static_cast<CellData *>(malloc(cell_num_tmp * sizeof(CellData)));
    cell_id_array_current_ = static_cast<unsigned int *>(malloc(cell_num_tmp * sizeof(unsigned int)));
    cell_id_to_index_ = static_cast<int *>(malloc((end_cell_id - start_cell_id) * sizeof(int)));
    memset(cell_id_to_index_, -1, (end_cell_id - start_cell_id) * sizeof(int));

    for (unsigned int y = min_block_y; y <= max_block_y; y++) {
        block_id_y = y * x_block_num;
        offset = block_index_[min_block_x + block_id_y];
        // new cell_num_tmp
        cell_num_tmp = block_index_[max_block_x + block_id_y + 1] - offset;
        selectCells(offset, cell_num_tmp, &cell_array_current_[cell_num_current_]);

        unsigned int start = cell_num_current_;
        for (unsigned int j = 0; j < cell_num_tmp; j++) {
            CellData cell = cell_array_current_[start + j];
            if (cell.x < min_x || cell.x > max_x || cell.y < min_y || cell.y > max_y)
                continue;
            memmove(&(cell_array_current_[cell_num_current_]), &cell, sizeof(CellData));
            cell_id_tmp = offset + j;
            cell_id_array_current_[cell_num_current_] = cell_id_tmp;
            cell_id_to_index_[cell_id_tmp-start_cell_id] = static_cast<int>(cell_num_current_);
            cell_num_current_++;
            expression_num_current_ += cell.gene_count;
        }
    }

    if (verbose_) printCpuTime(cprev, "restrictRegion");
}

void CgefReader::restrictGene(vector<string> &gene_list, bool exclude) {
    restrict_gene_ = true;
    auto *gene_id_bool_tmp = static_cast<bool*>(malloc(gene_num_ * sizeof(bool)));
    memset(gene_id_bool_tmp, exclude, gene_num_ * sizeof(bool));
    for (const auto &gene_name: gene_list) {
        if (genename_to_id_.find(gene_name) != genename_to_id_.end()) {
            gene_id_bool_tmp[genename_to_id_[gene_name]]  = !exclude;
        }
    }

    unsigned int new_gene_num_current = 0;
    for(unsigned int i = 0; i < gene_num_; i++){
        if(!gene_id_bool_tmp[i]) gene_id_to_index_[i] = -1;
        if(gene_id_to_index_[i]>=0) {
            gene_id_to_index_[i] = new_gene_num_current;
            new_gene_num_current++;
        }
    }
    gene_num_current_ = new_gene_num_current;
    free(gene_id_bool_tmp);
}

void CgefReader::updateGeneInfo() {
    if (isOldCellExpVersion)
    {
        auto *cell_exp_data = static_cast<olderCellExpData *>(malloc(gene_num_ * sizeof(olderCellExpData)));
        auto *gene_id_bool_tmp = static_cast<bool*>(calloc(gene_num_, sizeof(bool)));
        for (unsigned int i = 0; i < cell_num_current_; i++) {
            CellData cell = cell_array_current_[i];
            selectOlderCellExp(cell.offset, cell.gene_count, cell_exp_data);
            for (unsigned int j = 0; j < cell.gene_count; j++) {
                gene_id_bool_tmp[cell_exp_data[j].gene_id] = true;
            }
        }

        unsigned int new_gene_num_current = 0;
        for(unsigned int i = 0; i < gene_num_; i++){
            if(!gene_id_bool_tmp[i]) gene_id_to_index_[i] = -1;
            if(gene_id_to_index_[i]>=0) {
                gene_id_to_index_[i] = new_gene_num_current;
                new_gene_num_current++;
            }
        }

        gene_num_current_ = new_gene_num_current;
        free(cell_exp_data);
        free(gene_id_bool_tmp);
    } else {
        auto *cell_exp_data = static_cast<CellExpData *>(malloc(gene_num_ * sizeof(CellExpData)));
        auto *gene_id_bool_tmp = static_cast<bool*>(calloc(gene_num_, sizeof(bool)));
        for (unsigned int i = 0; i < cell_num_current_; i++) {
            CellData cell = cell_array_current_[i];
            selectCellExp(cell.offset, cell.gene_count, cell_exp_data);
            for (unsigned int j = 0; j < cell.gene_count; j++) {
                gene_id_bool_tmp[cell_exp_data[j].gene_id] = true;
            }
        }

        unsigned int new_gene_num_current = 0;
        for(unsigned int i = 0; i < gene_num_; i++){
            if(!gene_id_bool_tmp[i]) gene_id_to_index_[i] = -1;
            if(gene_id_to_index_[i]>=0) {
                gene_id_to_index_[i] = new_gene_num_current;
                new_gene_num_current++;
            }
        }

        gene_num_current_ = new_gene_num_current;
        free(cell_exp_data);
        free(gene_id_bool_tmp);
    }
}

unsigned int CgefReader::getCellCount(string &gene_name) {
    auto iter = genename_to_id_.find(gene_name);
    if (iter == genename_to_id_.end())
        return 0;
    return gene_array_[iter->second].cell_count;
}

unsigned int CgefReader::getCellCount(unsigned int gene_id) {
    if (gene_id >= gene_num_)
        return 0;
    return gene_array_[gene_id].cell_count;
}

unsigned short CgefReader::getGeneCount(unsigned int cell_id) const {
    if (cell_id >= cell_num_)
        return 0;
    return getCell(cell_id).gene_count;
}

void CgefReader::freeRestriction() {
    restrict_region_ = false;
    restrict_gene_ = false;
    if (cell_array_current_ != nullptr) {
        free(cell_array_current_);
        cell_array_current_ = nullptr;
    }
    if (cell_id_array_current_ != nullptr) {
        free(cell_id_array_current_);
        cell_id_array_current_ = nullptr;
    }
    if (cell_id_to_index_ != nullptr) {
        free(cell_id_to_index_);
        cell_id_to_index_ = nullptr;
    }
    iota(&gene_id_to_index_[0], gene_id_to_index_+gene_num_, 0);
    cell_num_current_ = cell_num_;
    gene_num_current_ = gene_num_;
    expression_num_current_ = expression_num_;
}

bool CgefReader::isInRegion(unsigned int cell_id) {
    if(cell_id < start_cell_id || cell_id >= end_cell_id || cell_id_to_index_[cell_id-start_cell_id] < 0)
        return false;
    return true;
}


void CgefReader::getCellBorders(vector<unsigned int> &cell_ind, vector<short> &border, vector<short> &borcnt)
{
    unsigned long cprev = clock();
    if(m_borderdataPtr_s == nullptr)
    {
        hsize_t dims[1];
        hid_t dataset_id = H5Dopen(group_id_, "cellBorder", H5P_DEFAULT);
        hid_t dataspace_id = H5Dget_space(dataset_id);
        H5Sget_simple_extent_dims(dataspace_id, dims, nullptr);
        m_bordercnt = dims[0];
        m_borderdataPtr_s = (short*)calloc(dims[0], 2);
        H5Dread(dataset_id, H5T_NATIVE_SHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_borderdataPtr_s);

        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);

        hid_t d_id = H5Dopen(group_id_, "cellBordercnt", H5P_DEFAULT);
        m_pborcnt = (short*)calloc(cell_num_, 2);
        H5Dread(d_id, H5T_NATIVE_SHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_pborcnt);
        H5Dclose(d_id);
    }

    vector<short> tmp(m_borderdataPtr_s, m_borderdataPtr_s+m_bordercnt);
    border.swap(tmp);
    vector<short> tmpcnt(m_pborcnt, m_pborcnt+cell_num_);
    borcnt.swap(tmpcnt);

}

int CgefReader::getCellBorders(vector<unsigned int> &cell_ind, vector<short> &border)
{
    unsigned long cprev = clock();
    if(m_borderdataPtr_s == nullptr)
    {
        hsize_t dims[3];
        hid_t dataset_id = H5Dopen(group_id_, "cellBorder", H5P_DEFAULT);
        hid_t dataspace_id = H5Dget_space(dataset_id);
        H5Sget_simple_extent_dims(dataspace_id, dims, nullptr);

        m_borderdataPtr_s = (short*)calloc(dims[0]*dims[1]*dims[2], 2);
        H5Dread(dataset_id, H5T_NATIVE_SHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_borderdataPtr_s);

        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        m_bordercnt = dims[1];
    }

    int cnt = 2*m_bordercnt;
    if(cell_ind.empty())
    {
        uint32_t total = cell_num_*m_bordercnt*2;
        vector<short> tmp(m_borderdataPtr_s, m_borderdataPtr_s+total);
        border.swap(tmp);
    }
    else
    {
        for(uint32_t cid : cell_ind)
        {
            short *psrc = m_borderdataPtr_s + cid*cnt;
            for(int i=0;i<cnt;i++)
            {
                border.push_back(psrc[i]);
            }
        }
    }
    return cnt;
}

void CgefReader::getfiltereddata(vector<int> &region, vector<string> &genelist,
                vector<string> &vec_gene, vector<unsigned long long> &uniq_cells,
                vector<unsigned int> &cell_ind, vector<unsigned int> &gene_ind, vector<unsigned int> &count)
{
    int min_x = 0, max_x = 0, min_y = 0, max_y = 0;
    if(!region.empty())
    {
        min_x = region[0];
        max_x = region[1];
        min_y = region[2];
        max_y = region[3];
    }

    CellData *cdata = loadCell();
    GeneData *gdata = loadGene();

    if(genelist.empty()&& (!region.empty()))
    {
        if (isOldCellExpVersion) {
            hid_t memtype = getMemtypeOfOlderCellExpData();
            olderCellExpData* cell_exp_data = (olderCellExpData*)malloc(expression_num_ * sizeof(olderCellExpData));
            H5Dread(cell_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell_exp_data);

            unsigned long long uniq_cell_id;
            uint32_t gid = 0, cid = 0;
            map<uint32_t, uint32_t> gidmap;//新旧gid对照
            for(uint32_t i=0;i<cell_num_;i++)
            {
                if(cdata[i].x < min_x || cdata[i].x >= max_x || 
                   cdata[i].y < min_y || cdata[i].y >= max_y){
                    continue;
                }

                olderCellExpData *cptr = cell_exp_data + cdata[i].offset;
                for(uint32_t j = 0;j<cdata[i].gene_count;j++)
                {
                    if(gidmap.find((uint32_t)cptr[j].gene_id) == gidmap.end())
                    {
                        gene_ind.push_back(gid);
                        string str(gdata[cptr[j].gene_id].gene_name);
                        vec_gene.emplace_back(std::move(str));
                        gidmap.emplace((uint32_t)cptr[j].gene_id, gid++);
                    }
                    else
                    {
                        gene_ind.push_back(gidmap[(uint32_t)cptr[j].gene_id]);
                    }

                    count.push_back((uint32_t)cptr[j].count);
                    cell_ind.push_back(cid);
                }
                uniq_cell_id = cdata[i].x;
                uniq_cell_id = (uniq_cell_id << 32) | cdata[i].y;
                uniq_cells.emplace_back(uniq_cell_id);
                cid++;
            }
            free(cell_exp_data);
        } else {
            hid_t memtype = getMemtypeOfCellExpData();
            CellExpData* cell_exp_data = (CellExpData*)malloc(expression_num_ * sizeof(CellExpData));
            H5Dread(cell_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell_exp_data);
    
            unsigned long long uniq_cell_id;
            uint32_t gid = 0, cid = 0;
            map<uint32_t, uint32_t> gidmap;//新旧gid对照
            for(uint32_t i=0;i<cell_num_;i++)
            {
                if(cdata[i].x < min_x || cdata[i].x >= max_x || 
                   cdata[i].y < min_y || cdata[i].y >= max_y){
                    continue;
                }
    
                CellExpData *cptr = cell_exp_data + cdata[i].offset;
                for(uint32_t j = 0;j<cdata[i].gene_count;j++)
                {
                    if(gidmap.find(cptr[j].gene_id) == gidmap.end())
                    {
                        gene_ind.push_back(gid);
                        string str(gdata[cptr[j].gene_id].gene_name);
                        vec_gene.emplace_back(std::move(str));
                        gidmap.emplace(cptr[j].gene_id, gid++);
                    }
                    else
                    {
                        gene_ind.push_back(gidmap[cptr[j].gene_id]);
                    }
                    
                    count.push_back(cptr[j].count);
                    cell_ind.push_back(cid);
                }
                uniq_cell_id = cdata[i].x;
                uniq_cell_id = (uniq_cell_id << 32) | cdata[i].y;
                uniq_cells.emplace_back(uniq_cell_id);
                cid++;
            }
            free(cell_exp_data);
        }
    }
    else if(region.empty() && (!genelist.empty()))
    {
        set<string> gset;
        for(string &str : genelist)
        {
            gset.insert(str);
        }
        hid_t memtype = getMemtypeOfGeneExpData();
        GeneExpData* gene_exp_data = (GeneExpData*)malloc(expression_num_ * sizeof(GeneExpData));
        H5Dread(gene_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, gene_exp_data);

        unsigned long long uniq_cell_id;
        uint32_t gid = 0, cid = 0, oldcid = 0;
        map<uint32_t, uint32_t> cidmap;//新旧cid对照

        for(unsigned int gene_id = 0; gene_id < gene_num_; gene_id++)
        {
            string str(gdata[gene_id].gene_name);
            if(gset.find(str) != gset.end())
            {
                vec_gene.emplace_back(std::move(str));
                GeneExpData *gptr = gene_exp_data + gdata[gene_id].offset;
                for(uint32_t j = 0;j<gdata[gene_id].cell_count;j++)
                {
                    oldcid = gptr[j].cell_id;
                    if(cidmap.find(oldcid) == cidmap.end())
                    {
                        cell_ind.push_back(cid);
                        cidmap.emplace(oldcid, cid++);

                        uniq_cell_id = cdata[oldcid].x;
                        uniq_cell_id = (uniq_cell_id << 32) | cdata[oldcid].y;
                        uniq_cells.emplace_back(uniq_cell_id);
                    }
                    else
                    {
                        cell_ind.push_back(cidmap[oldcid]);
                    }

                    count.push_back(gptr[j].count);
                    gene_ind.push_back(gid);
                }
                gid++;
            }
        }

        free(gene_exp_data);
    }
    else if((!region.empty()) && (!genelist.empty()))
    {
        set<string> gset;
        for(string &str : genelist)
        {
            gset.insert(str);
        }
        hid_t memtype = getMemtypeOfGeneExpData();
        GeneExpData* gene_exp_data = (GeneExpData*)malloc(expression_num_ * sizeof(GeneExpData));
        H5Dread(gene_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, gene_exp_data);

        unsigned long long uniq_cell_id;
        uint32_t gid = 0, cid = 0, oldcid = 0;
        map<uint32_t, uint32_t> cidmap;//新旧cid对照
        for(unsigned int gene_id = 0; gene_id < gene_num_; gene_id++)
        {
            string str(gdata[gene_id].gene_name);
            if(gset.find(str) != gset.end())
            {
                vec_gene.emplace_back(std::move(str));
                GeneExpData *gptr = gene_exp_data + gdata[gene_id].offset;
                for(uint32_t j=0;j<gdata[gene_id].cell_count;j++)
                {
                    oldcid = gptr[j].cell_id;
                    if(cdata[oldcid].x < min_x || cdata[oldcid].x >= max_x || 
                        cdata[oldcid].y < min_y || cdata[oldcid].y >= max_y)
                    {
                        continue;
                    }

                    if(cidmap.find(oldcid) == cidmap.end())
                    {
                        cell_ind.push_back(cid);
                        cidmap.emplace(oldcid, cid++);

                        uniq_cell_id = cdata[oldcid].x;
                        uniq_cell_id = (uniq_cell_id << 32) | cdata[oldcid].y;
                        uniq_cells.emplace_back(uniq_cell_id);
                    }
                    else
                    {
                        cell_ind.push_back(cidmap[oldcid]);
                    }

                    
                    count.push_back(gptr[j].count);
                    gene_ind.push_back(gid);
                }
                gid++;
            }
        }
        free(gene_exp_data);
    }
    else if(region.empty() && genelist.empty())
    {
        if (isOldCellExpVersion) {
            hid_t memtype = getMemtypeOfOlderCellExpData();
            olderCellExpData* cell_exp_data = (olderCellExpData*)malloc(expression_num_ * sizeof(olderCellExpData));
            H5Dread(cell_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell_exp_data);

            for(uint32_t i=0;i<gene_num_;i++)
            {
                vec_gene.emplace_back(gene_array_[i].gene_name);
            }

            unsigned long long uniq_cell_id;
            for(uint32_t i=0;i<cell_num_;i++)
            {
                olderCellExpData *cptr = cell_exp_data + cdata[i].offset;
                for(uint32_t j = 0;j<cdata[i].gene_count;j++)
                {
                    gene_ind.push_back(cptr[j].gene_id);
                    count.push_back(cptr[j].count);
                    cell_ind.push_back(i);
                }
                uniq_cell_id = cdata[i].x;
                uniq_cell_id = (uniq_cell_id << 32) | cdata[i].y;
                uniq_cells.emplace_back(uniq_cell_id);
            }
            free(cell_exp_data);
        } else {
            hid_t memtype = getMemtypeOfCellExpData();
            CellExpData* cell_exp_data = (CellExpData*)malloc(expression_num_ * sizeof(CellExpData));
            H5Dread(cell_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell_exp_data);

            for(uint32_t i=0;i<gene_num_;i++)
            {
                vec_gene.emplace_back(gene_array_[i].gene_name);
            }

            unsigned long long uniq_cell_id;
            for(uint32_t i=0;i<cell_num_;i++)
            {
                CellExpData *cptr = cell_exp_data + cdata[i].offset;
                for(uint32_t j = 0;j<cdata[i].gene_count;j++)
                {
                    gene_ind.push_back(cptr[j].gene_id);
                    count.push_back(cptr[j].count);
                    cell_ind.push_back(i);
                }
                uniq_cell_id = cdata[i].x;
                uniq_cell_id = (uniq_cell_id << 32) | cdata[i].y;
                uniq_cells.emplace_back(uniq_cell_id);
            }
            free(cell_exp_data);
        }
    }
}


void CgefReader::getfiltereddata_exon(vector<int> &region, vector<string> &genelist,
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

    CellData *cdata = loadCell();
    GeneData *gdata = loadGene();

    if(genelist.empty()&& (!region.empty()))
    {
        if (isOldCellExpVersion) {
            hid_t memtype = getMemtypeOfOlderCellExpData();
            olderCellExpData* cell_exp_data = (olderCellExpData*)malloc(expression_num_ * sizeof(olderCellExpData));
            H5Dread(cell_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell_exp_data);

            uint16_t *pcell_exp_exon = (uint16_t *)malloc(2*expression_num_);
            hid_t did = H5Dopen(group_id_, "cellExpExon", H5P_DEFAULT);
            H5Dread(did, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, pcell_exp_exon);
            H5Dclose(did);

            unsigned long long uniq_cell_id;
            uint32_t gid = 0, cid = 0;
            map<uint32_t, uint32_t> gidmap;//新旧gid对照
            for(uint32_t i=0;i<cell_num_;i++)
            {
                if(cdata[i].x < min_x || cdata[i].x >= max_x || 
                   cdata[i].y < min_y || cdata[i].y >= max_y){
                    continue;
                }

                uint16_t *pexon = pcell_exp_exon + cdata[i].offset;
                olderCellExpData *cptr = cell_exp_data + cdata[i].offset;
                for(uint32_t j = 0;j<cdata[i].gene_count;j++)
                {
                    if(gidmap.find(cptr[j].gene_id) == gidmap.end())
                    {
                        gene_ind.push_back(gid);
                        string str(gdata[cptr[j].gene_id].gene_name);
                        vec_gene.emplace_back(std::move(str));
                        gidmap.emplace(cptr[j].gene_id, gid++);
                    }
                    else
                    {
                        gene_ind.push_back(gidmap[cptr[j].gene_id]);
                    }

                    exon.push_back(pexon[j]);
                    count.push_back(cptr[j].count);
                    cell_ind.push_back(cid);
                }
                uniq_cell_id = cdata[i].x;
                uniq_cell_id = (uniq_cell_id << 32) | cdata[i].y;
                uniq_cells.emplace_back(uniq_cell_id);
                cid++;
            }
            free(cell_exp_data);
            free(pcell_exp_exon);
        } else {
            hid_t memtype = getMemtypeOfCellExpData();
            CellExpData* cell_exp_data = (CellExpData*)malloc(expression_num_ * sizeof(CellExpData));
            H5Dread(cell_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell_exp_data);

            uint16_t *pcell_exp_exon = (uint16_t *)malloc(2*expression_num_);
            hid_t did = H5Dopen(group_id_, "cellExpExon", H5P_DEFAULT);
            H5Dread(did, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, pcell_exp_exon);
            H5Dclose(did);

            unsigned long long uniq_cell_id;
            uint32_t gid = 0, cid = 0;
            map<uint32_t, uint32_t> gidmap;//新旧gid对照
            for(uint32_t i=0;i<cell_num_;i++)
            {
                if(cdata[i].x < min_x || cdata[i].x >= max_x || 
                   cdata[i].y < min_y || cdata[i].y >= max_y){
                    continue;
                }

                uint16_t *pexon = pcell_exp_exon + cdata[i].offset;
                CellExpData *cptr = cell_exp_data + cdata[i].offset;
                for(uint32_t j = 0;j<cdata[i].gene_count;j++)
                {
                    if(gidmap.find(cptr[j].gene_id) == gidmap.end())
                    {
                        gene_ind.push_back(gid);
                        string str(gdata[cptr[j].gene_id].gene_name);
                        vec_gene.emplace_back(std::move(str));
                        gidmap.emplace(cptr[j].gene_id, gid++);
                    }
                    else
                    {
                        gene_ind.push_back(gidmap[cptr[j].gene_id]);
                    }

                    exon.push_back(pexon[j]);
                    count.push_back(cptr[j].count);
                    cell_ind.push_back(cid);
                }
                uniq_cell_id = cdata[i].x;
                uniq_cell_id = (uniq_cell_id << 32) | cdata[i].y;
                uniq_cells.emplace_back(uniq_cell_id);
                cid++;
            }
            free(cell_exp_data);
            free(pcell_exp_exon);
        }

    }
    else if(region.empty() && (!genelist.empty()))
    {
        set<string> gset;
        for(string &str : genelist)
        {
            gset.insert(str);
        }
        hid_t memtype = getMemtypeOfGeneExpData();
        GeneExpData* gene_exp_data = (GeneExpData*)malloc(expression_num_ * sizeof(GeneExpData));
        H5Dread(gene_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, gene_exp_data);

        uint16_t *pgene_exp_exon = (uint16_t *)malloc(2*expression_num_);
        hid_t did = H5Dopen(group_id_, "geneExpExon", H5P_DEFAULT);
        H5Dread(did, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, pgene_exp_exon);
        H5Dclose(did);

        unsigned long long uniq_cell_id;
        uint32_t gid = 0, cid = 0, oldcid = 0;
        map<uint32_t, uint32_t> cidmap;//新旧cid对照

        for(unsigned int gene_id = 0; gene_id < gene_num_; gene_id++)
        {
            string str(gdata[gene_id].gene_name);
            if(gset.find(str) != gset.end())
            {
                vec_gene.emplace_back(std::move(str));
                GeneExpData *gptr = gene_exp_data + gdata[gene_id].offset;
                uint16_t *pexon = pgene_exp_exon + gdata[gene_id].offset;
                for(uint32_t j = 0;j<gdata[gene_id].cell_count;j++)
                {
                    oldcid = gptr[j].cell_id;
                    if(cidmap.find(oldcid) == cidmap.end())
                    {
                        cell_ind.push_back(cid);
                        cidmap.emplace(oldcid, cid++);

                        uniq_cell_id = cdata[oldcid].x;
                        uniq_cell_id = (uniq_cell_id << 32) | cdata[oldcid].y;
                        uniq_cells.emplace_back(uniq_cell_id);
                    }
                    else
                    {
                        cell_ind.push_back(cidmap[oldcid]);
                    }

                    exon.push_back(pexon[j]);
                    count.push_back(gptr[j].count);
                    gene_ind.push_back(gid);
                }
                gid++;
            }
        }

        free(gene_exp_data);
        free(pgene_exp_exon);
    }
    else if((!region.empty()) && (!genelist.empty()))
    {
        set<string> gset;
        for(string &str : genelist)
        {
            gset.insert(str);
        }
        hid_t memtype = getMemtypeOfGeneExpData();
        GeneExpData* gene_exp_data = (GeneExpData*)malloc(expression_num_ * sizeof(GeneExpData));
        H5Dread(gene_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, gene_exp_data);

        uint16_t *pgene_exp_exon = (uint16_t *)malloc(2*expression_num_);
        hid_t did = H5Dopen(group_id_, "geneExpExon", H5P_DEFAULT);
        H5Dread(did, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, pgene_exp_exon);
        H5Dclose(did);

        unsigned long long uniq_cell_id;
        uint32_t gid = 0, cid = 0, oldcid = 0;
        map<uint32_t, uint32_t> cidmap;//新旧cid对照
        for(unsigned int gene_id = 0; gene_id < gene_num_; gene_id++)
        {
            string str(gdata[gene_id].gene_name);
            if(gset.find(str) != gset.end())
            {
                vec_gene.emplace_back(std::move(str));
                GeneExpData *gptr = gene_exp_data + gdata[gene_id].offset;
                uint16_t *pexon = pgene_exp_exon + gdata[gene_id].offset;
                for(uint32_t j=0;j<gdata[gene_id].cell_count;j++)
                {
                    oldcid = gptr[j].cell_id;
                    if(cdata[oldcid].x < min_x || cdata[oldcid].x >= max_x || 
                        cdata[oldcid].y < min_y || cdata[oldcid].y >= max_y)
                    {
                        continue;
                    }

                    if(cidmap.find(oldcid) == cidmap.end())
                    {
                        cell_ind.push_back(cid);
                        cidmap.emplace(oldcid, cid++);

                        uniq_cell_id = cdata[oldcid].x;
                        uniq_cell_id = (uniq_cell_id << 32) | cdata[oldcid].y;
                        uniq_cells.emplace_back(uniq_cell_id);
                    }
                    else
                    {
                        cell_ind.push_back(cidmap[oldcid]);
                    }

                    exon.push_back(pexon[j]);
                    count.push_back(gptr[j].count);
                    gene_ind.push_back(gid);
                }
                gid++;
            }
        }
        free(gene_exp_data);
        free(pgene_exp_exon);
    }
    else if(region.empty() && genelist.empty())
    {
        if (isOldCellExpVersion) {
            hid_t memtype = getMemtypeOfOlderCellExpData();
            olderCellExpData* cell_exp_data = (olderCellExpData*)malloc(expression_num_ * sizeof(olderCellExpData));
            H5Dread(cell_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell_exp_data);
    
            uint16_t *pcell_exp_exon = (uint16_t *)malloc(2*expression_num_);
            hid_t did = H5Dopen(group_id_, "cellExpExon", H5P_DEFAULT);
            H5Dread(did, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, pcell_exp_exon);
            H5Dclose(did);
    
            for(uint32_t i=0;i<gene_num_;i++)
            {
                vec_gene.emplace_back(gene_array_[i].gene_name);
            }
            
            unsigned long long uniq_cell_id;
            for(uint32_t i=0;i<cell_num_;i++)
            {
                olderCellExpData *cptr = cell_exp_data + cdata[i].offset;
                for(uint32_t j = 0;j<cdata[i].gene_count;j++)
                {
                    gene_ind.push_back(cptr[j].gene_id);
                    count.push_back(cptr[j].count);
                    exon.push_back(pcell_exp_exon[j]);
                    cell_ind.push_back(i);
                }
                uniq_cell_id = cdata[i].x;
                uniq_cell_id = (uniq_cell_id << 32) | cdata[i].y;
                uniq_cells.emplace_back(uniq_cell_id);
            }
            free(cell_exp_data);
            free(pcell_exp_exon);
        } else {
            hid_t memtype = getMemtypeOfCellExpData();
            CellExpData* cell_exp_data = (CellExpData*)malloc(expression_num_ * sizeof(CellExpData));
            H5Dread(cell_exp_dataset_id_, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, cell_exp_data);

            uint16_t *pcell_exp_exon = (uint16_t *)malloc(2*expression_num_);
            hid_t did = H5Dopen(group_id_, "cellExpExon", H5P_DEFAULT);
            H5Dread(did, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, pcell_exp_exon);
            H5Dclose(did);

            for(uint32_t i=0;i<gene_num_;i++)
            {
                vec_gene.emplace_back(gene_array_[i].gene_name);
            }

            unsigned long long uniq_cell_id;
            for(uint32_t i=0;i<cell_num_;i++)
            {
                CellExpData *cptr = cell_exp_data + cdata[i].offset;
                for(uint32_t j = 0;j<cdata[i].gene_count;j++)
                {
                    gene_ind.push_back(cptr[j].gene_id);
                    count.push_back(cptr[j].count);
                    exon.push_back(pcell_exp_exon[j]);
                    cell_ind.push_back(i);
                }
                uniq_cell_id = cdata[i].x;
                uniq_cell_id = (uniq_cell_id << 32) | cdata[i].y;
                uniq_cells.emplace_back(uniq_cell_id);
            }
            free(cell_exp_data);
            free(pcell_exp_exon);
        }

    }
}