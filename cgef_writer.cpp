#include "cgef_writer.h"
#include <stdlib.h>
#include <time.h> 
#include <random>

CgefWriter::CgefWriter( bool verbose) {
    str32_type_ = H5Tcopy(H5T_C_S1);
    H5Tset_size(str32_type_, 32);
    str64_type_ = H5Tcopy(H5T_C_S1);
    H5Tset_size(str64_type_, 64);
    verbose_ = verbose;

    // cerr << "create h5 file: " <<  output_cell_gef << endl;
    // hid_t fapl_id = H5Pcreate (H5P_FILE_ACCESS);
    // H5Pset_libver_bounds(fapl_id, H5F_LIBVER_V18, H5F_LIBVER_LATEST);
    // H5Pset_fclose_degree(fapl_id, H5F_CLOSE_STRONG);
    // file_id_ = H5Fcreate(output_cell_gef.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
    // group_id_ = H5Gcreate(file_id_, "/cellBin", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    // H5Pclose(fapl_id);
}

void CgefWriter::setOutput(const string& output_cell_gef)
{
    cerr << "create h5 file: " <<  output_cell_gef << endl;
    hid_t fapl_id = H5Pcreate (H5P_FILE_ACCESS);
    H5Pset_libver_bounds(fapl_id, H5F_LIBVER_V18, H5F_LIBVER_LATEST);
    H5Pset_fclose_degree(fapl_id, H5F_CLOSE_STRONG);
    file_id_ = H5Fcreate(output_cell_gef.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
    group_id_ = H5Gcreate(file_id_, "/cellBin", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Pclose(fapl_id);
}

void CgefWriter::setInput(const string& input_cell_gef)
{
    cerr << "open h5 file: " <<  input_cell_gef << endl;
    hid_t fapl_id = H5Pcreate (H5P_FILE_ACCESS);
    H5Pset_libver_bounds(fapl_id, H5F_LIBVER_V18, H5F_LIBVER_LATEST);
    H5Pset_fclose_degree(fapl_id, H5F_CLOSE_STRONG);
    file_id_ = H5Fopen(input_cell_gef.c_str(), H5F_ACC_RDWR, fapl_id);
    group_id_ = H5Gopen(file_id_, "/cellBin", H5P_DEFAULT);

    H5Pclose(fapl_id);
    openCellDataset();
    getAttr();
}

CgefWriter::~CgefWriter() {
    H5Tclose(str32_type_);
    H5Tclose(str64_type_);
    H5Gclose(group_id_);
    H5Fclose(file_id_);
}

void CgefWriter::openCellDataset() 
{
    unsigned long cprev = clock();
    hid_t cell_dataset_id = H5Dopen(group_id_, "cell", H5P_DEFAULT);
    if (cell_dataset_id < 0) {
        cerr << "failed open dataset: cell" << endl;
        reportErrorCode2File(errorCode::E_MISSINGFILEINFO, "failed open dataset: cell");
        exit(3);
    }

    hid_t s1_tid = H5Dget_type(cell_dataset_id);
    int nmemb = H5Tget_nmembers(s1_tid);
    if(nmemb < 9){
        cerr << "Please use geftools(>=0.6) to regenerate this cgef file." << endl;
        reportErrorCode2File(errorCode::E_LOWVERSION, 
                            "Please use geftools(>=0.6) to regenerate this cgef file.");
        exit(2);
    }

    hsize_t dims[1];
    hid_t cell_dataspace_id = H5Dget_space(cell_dataset_id);
    H5Sget_simple_extent_dims(cell_dataspace_id, dims, nullptr);
    cell_num_ = dims[0];

    hid_t memtype = getMemtypeOfCellData();
    m_cdataPtr = (CellData *) malloc(cell_num_ * sizeof(CellData));
    H5Dread(cell_dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_cdataPtr);


    hid_t attr = H5Aopen(cell_dataset_id, "minX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT32, &m_canvas[0]);

    attr = H5Aopen(cell_dataset_id, "minY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT32, &m_canvas[1]);

    attr = H5Aopen(cell_dataset_id, "maxX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT32, &m_canvas[2]);

    attr = H5Aopen(cell_dataset_id, "maxY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT32, &m_canvas[3]);

    H5Aclose(attr);

    H5Sclose(cell_dataspace_id);
    H5Dclose(cell_dataset_id);

    if (verbose_) printCpuTime(cprev, "openCellDataset");
}

void CgefWriter::getAttr()
{
    hid_t attr = H5Aopen(file_id_, "offsetX", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT32, &m_offsetX);

    attr = H5Aopen(file_id_, "offsetY", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT32, &m_offsetY);

    H5Aclose(attr);
}

void CgefWriter::storeCellBorder(short* borderPath, unsigned int cell_num) const {
    unsigned long cprev=clock();
    hsize_t dims[3];
    dims[0] = cell_num;
    dims[1] = BORDERCNT;
    dims[2] = 2;

    hid_t dataspace_id = H5Screate_simple(3, dims, nullptr);
    hid_t dataset_id = H5Dcreate(group_id_, "cellBorder", H5T_STD_I16LE, dataspace_id, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_STD_I16LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, borderPath);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    if(verbose_) printCpuTime(cprev, "storeCellBorder");
}

// void CgefWriter::storeCellBorder(short* borderPath, unsigned int cell_num) const {
//     unsigned long cprev=clock();
//     hsize_t dims[1] = {cell_num};
//     // dims[0] = cell_num;
//     // dims[1] = BORDERCNT;
//     // dims[2] = 2;

//     hid_t dataspace_id = H5Screate_simple(1, dims, nullptr);
//     hid_t dataset_id = H5Dcreate(group_id_, "cellBorder", H5T_STD_I16LE, dataspace_id, H5P_DEFAULT,
//                                  H5P_DEFAULT, H5P_DEFAULT);
//     H5Dwrite(dataset_id, H5T_STD_I16LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, borderPath);
//     H5Sclose(dataspace_id);
//     H5Dclose(dataset_id);
//     if(verbose_) printCpuTime(cprev, "storeCellBorder");
// }

void CgefWriter::storeCellBorder_cnt(vector<short> &borcnt) {
    unsigned long cprev=clock();
    hsize_t dims[1] = {borcnt.size()};


    hid_t dataspace_id = H5Screate_simple(1, dims, nullptr);
    hid_t dataset_id = H5Dcreate(group_id_, "cellBordercnt", H5T_STD_I16LE, dataspace_id, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_STD_I16LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, borcnt.data());
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    if(verbose_) printCpuTime(cprev, "storeCellBordercnt");
}


void CgefWriter::storeCellBorderWithAttr(short *borderPath, unsigned int cell_num, int *effective_rect) const {
    unsigned long cprev=clock();
    storeCellBorder(borderPath, cell_num);

    hid_t dataset_id = H5Dopen(group_id_, "cellBorder", H5P_DEFAULT);

    hsize_t dims_attr[1] = {1};
    hid_t attr;
    hid_t attr_dataspace = H5Screate_simple(1, dims_attr, nullptr);
    attr = H5Acreate(dataset_id, "minX", H5T_STD_I32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT, effective_rect);
    attr = H5Acreate(dataset_id, "minY", H5T_STD_I32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT, &effective_rect[1]);
    attr = H5Acreate(dataset_id, "maxX", H5T_STD_I32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT, &effective_rect[2]);
    attr = H5Acreate(dataset_id, "maxY", H5T_STD_I32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT, &effective_rect[3]);

    H5Aclose(attr);
    H5Sclose(attr_dataspace);
    H5Dclose(dataset_id);
    if(verbose_) printCpuTime(cprev, "storeCellBorderWithAttr");
}

void CgefWriter::storeCellExp() {
    unsigned long cprev=clock();
    hsize_t dims[1] = {cell_exp_list_.size()};

    hid_t memtype, filetype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(CellExpData));
    H5Tinsert(memtype, "geneID", HOFFSET(CellExpData, gene_id), H5T_NATIVE_UINT32);
    H5Tinsert(memtype, "count", HOFFSET(CellExpData, count), H5T_NATIVE_USHORT);

    filetype = H5Tcreate(H5T_COMPOUND, 6);
    H5Tinsert(filetype, "geneID", 0, H5T_STD_U32LE);
    H5Tinsert(filetype, "count", 4, H5T_STD_U16LE);

    hid_t dataspace_id = H5Screate_simple(1, dims, nullptr);
    hid_t dataset_id = H5Dcreate(group_id_, "cellExp", filetype, dataspace_id, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &cell_exp_list_[0]);

    hsize_t dims_attr[1] = {1};
    hid_t attr;
    hid_t attr_dataspace = H5Screate_simple(1, dims_attr, nullptr);
    attr = H5Acreate(dataset_id, "maxCount", H5T_STD_U16LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_USHORT, &max_mid_count_);

    H5Aclose(attr);
    H5Sclose(attr_dataspace);
    H5Tclose(memtype);
    H5Tclose(filetype);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    if(verbose_) printCpuTime(cprev, "storeCellExp");
}

void CgefWriter::storeCell(unsigned int block_num, unsigned int * block_index, const unsigned int *block_size) {
    unsigned long cprev=clock();
    m_cdataPtr = cell_list_.data();
    hsize_t dims[1] = {(hsize_t)cell_num_};

    hid_t memtype, filetype;
    memtype = getMemtypeOfCellData();

    filetype = H5Tcreate(H5T_COMPOUND, sizeof(CellData));
    H5Tinsert(filetype, "id", HOFFSET(CellData, id), H5T_STD_U32LE);
    H5Tinsert(filetype, "x", HOFFSET(CellData, x), H5T_STD_I32LE);
    H5Tinsert(filetype, "y", HOFFSET(CellData, y), H5T_STD_I32LE);
    H5Tinsert(filetype, "offset", HOFFSET(CellData, offset), H5T_STD_U32LE);
    H5Tinsert(filetype, "geneCount", HOFFSET(CellData, gene_count), H5T_STD_U16LE);
    H5Tinsert(filetype, "expCount", HOFFSET(CellData, exp_count), H5T_STD_U16LE);
    H5Tinsert(filetype, "dnbCount", HOFFSET(CellData, dnb_count), H5T_STD_U16LE);
    H5Tinsert(filetype, "area", HOFFSET(CellData, area), H5T_STD_U16LE);
    H5Tinsert(filetype, "cellTypeID", HOFFSET(CellData, cell_type_id), H5T_STD_U16LE);
    H5Tinsert(filetype, "clusterID", HOFFSET(CellData, cluster_id), H5T_STD_U16LE);
    //H5Tinsert(filetype, "incnt", HOFFSET(CellData, incnt), H5T_STD_U16LE);
    
    hid_t dataspace_id = H5Screate_simple(1, dims, nullptr);

    hid_t dpid = H5Pcreate (H5P_DATASET_CREATE);
    H5Pset_attr_phase_change(dpid, 0, 0);
    hid_t dataset_id = H5Dcreate(group_id_, "cell", filetype, dataspace_id, H5P_DEFAULT,
                                 dpid, H5P_DEFAULT);
    H5Dwrite(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &cell_list_[0]);

    // Create cell attribute
    // median
    auto * index = (unsigned int *) malloc(cell_num_ * sizeof(unsigned int));
    iota(index, index+cell_num_, 0);
    sort(index, index+cell_num_,[this](int a,int b){return cell_list_[a].area < cell_list_[b].area; });

    unsigned int mi = cell_num_ / 2;
    if (cell_num_ % 2 != 0) {
        cell_attr_.median_area = cell_list_[index[mi]].area;
    } else {
        cell_attr_.median_area = float(cell_list_[index[mi]].area + cell_list_[index[mi - 1]].area)/2;
    }

    iota(index, index+cell_num_, 0);
    sort(index, index+cell_num_,
         [this](int a,int b){return cell_list_[a].gene_count < cell_list_[b].gene_count; });
    if (cell_num_ % 2 != 0) {
        cell_attr_.median_gene_count = cell_list_[index[mi]].gene_count;
    } else {
        cell_attr_.median_gene_count = float(cell_list_[index[mi]].gene_count + cell_list_[index[mi - 1]].gene_count)/2;
    }

    iota(index, index+cell_num_, 0);
    sort(index, index+cell_num_,
         [this](int a,int b){return cell_list_[a].exp_count < cell_list_[b].exp_count; });
    if (cell_num_ % 2 != 0) {
        cell_attr_.median_exp_count = cell_list_[index[mi]].exp_count;
    } else {
        cell_attr_.median_exp_count = float(cell_list_[index[mi]].exp_count + cell_list_[index[mi - 1]].exp_count)/2;
    }

    iota(index, index+cell_num_, 0);
    sort(index, index+cell_num_,
         [this](int a,int b){return cell_list_[a].dnb_count < cell_list_[b].dnb_count; });
    if (cell_num_ % 2 != 0) {
        cell_attr_.median_dnb_count = cell_list_[index[mi]].dnb_count;
    } else {
        cell_attr_.median_dnb_count = float(cell_list_[index[mi]].dnb_count + cell_list_[index[mi - 1]].dnb_count)/2;
    }

    // average
    cell_attr_.average_gene_count = static_cast<float>(expression_num_)/  static_cast<float>(cell_num_);
    cell_attr_.average_exp_count = static_cast<float>(exp_count_sum_) /  static_cast<float>(cell_num_);
    cell_attr_.average_dnb_count = static_cast<float>(dnb_count_sum_) /  static_cast<float>(cell_num_);
    cell_attr_.average_area = static_cast<float>(area_sum_) /  static_cast<float>(cell_num_);

    hsize_t dimsAttr[1] = {1};
    hid_t attr;
    hid_t attr_dataspace = H5Screate_simple(1, dimsAttr, nullptr);
    attr = H5Acreate(dataset_id, "averageGeneCount", H5T_IEEE_F32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_FLOAT, &cell_attr_.average_gene_count);
    attr = H5Acreate(dataset_id, "averageExpCount", H5T_IEEE_F32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_FLOAT, &cell_attr_.average_exp_count);
    attr = H5Acreate(dataset_id, "averageDnbCount", H5T_IEEE_F32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_FLOAT, &cell_attr_.average_dnb_count);
    attr = H5Acreate(dataset_id, "averageArea", H5T_IEEE_F32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_FLOAT, &cell_attr_.average_area);
    attr = H5Acreate(dataset_id, "medianGeneCount", H5T_IEEE_F32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_FLOAT, &cell_attr_.median_gene_count);
    attr = H5Acreate(dataset_id, "medianExpCount", H5T_IEEE_F32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_FLOAT, &cell_attr_.median_exp_count);
    attr = H5Acreate(dataset_id, "medianDnbCount", H5T_IEEE_F32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_FLOAT, &cell_attr_.median_dnb_count);
    attr = H5Acreate(dataset_id, "medianArea", H5T_IEEE_F32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_FLOAT, &cell_attr_.median_area);
    attr = H5Acreate(dataset_id, "minX", H5T_STD_I32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT32, &cell_attr_.min_x);
    attr = H5Acreate(dataset_id, "maxX", H5T_STD_I32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT32, &cell_attr_.max_x);
    attr = H5Acreate(dataset_id, "minY", H5T_STD_I32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT32, &cell_attr_.min_y);
    attr = H5Acreate(dataset_id, "maxY", H5T_STD_I32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT32, &cell_attr_.max_y);
    attr = H5Acreate(dataset_id, "minGeneCount", H5T_STD_U16LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_USHORT, &cell_attr_.min_gene_count);
    attr = H5Acreate(dataset_id, "minExpCount", H5T_STD_U16LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_USHORT, &cell_attr_.min_exp_count);
    attr = H5Acreate(dataset_id, "minDnbCount", H5T_STD_U16LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_USHORT, &cell_attr_.min_dnb_count);
    attr = H5Acreate(dataset_id, "minArea", H5T_STD_U16LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_USHORT, &cell_attr_.min_area);
    attr = H5Acreate(dataset_id, "maxGeneCount", H5T_STD_U16LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_USHORT, &cell_attr_.max_gene_count);
    attr = H5Acreate(dataset_id, "maxExpCount", H5T_STD_U16LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_USHORT, &cell_attr_.max_exp_count);
    attr = H5Acreate(dataset_id, "maxDnbCount", H5T_STD_U16LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_USHORT, &cell_attr_.max_dnb_count);
    attr = H5Acreate(dataset_id, "maxArea", H5T_STD_U16LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_USHORT, &cell_attr_.max_area);
    // attr = H5Acreate(dataset_id, "col", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    // H5Awrite(attr, H5T_NATIVE_USHORT, m_x_len);
    // attr = H5Acreate(dataset_id, "row", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    // H5Awrite(attr, H5T_NATIVE_USHORT, m_y_len);

    // write block index
   // dimsAttr[0] = block_num + 1;
    // attr_dataspace = H5Screate_simple(1, dimsAttr, nullptr);
    // attr = H5Acreate(dataset_id, "blockIndex", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    // H5Awrite(attr, H5T_NATIVE_UINT32, block_index);

    // dimsAttr[0] = 4;
    // attr_dataspace = H5Screate_simple(1, dimsAttr, nullptr);
    // attr = H5Acreate(dataset_id, "blockSize", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    // H5Awrite(attr, H5T_NATIVE_UINT32, block_size);

    H5Aclose(attr);
    H5Tclose(memtype);
    H5Tclose(filetype);
    H5Sclose(attr_dataspace);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    free(index);

    storeBlkidx(block_num, block_index, block_size);
    if(verbose_) printCpuTime(cprev, "storeCell");
}

void CgefWriter::storeBlkidx(unsigned int block_num, unsigned int * block_index, const unsigned int *block_size)
{
    hsize_t dims[1] = {block_num + 1};
    hid_t ds = H5Screate_simple(1, dims, nullptr);
    hid_t dataset_id = H5Dcreate(group_id_, "blockIndex", H5T_STD_U32LE, ds, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, block_index);

    H5Sclose(ds);
    H5Dclose(dataset_id);

    dims[0] = 4;
    hid_t sds = H5Screate_simple(1, dims, nullptr);
    hid_t sdataset_id = H5Dcreate(group_id_, "blockSize", H5T_STD_U32LE, sds, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(sdataset_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, block_size);

    H5Sclose(sds);
    H5Dclose(sdataset_id);
}

void CgefWriter::storeCellLabel(vector<unsigned int> &vecdata)
{
    hsize_t dims[1] = {vecdata.size()};
    hid_t ds = H5Screate_simple(1, dims, nullptr);
    hid_t dataset_id = H5Dcreate(group_id_, "label", H5T_STD_U32LE, ds, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, vecdata.data());

    H5Sclose(ds);
    H5Dclose(dataset_id);
}

void CgefWriter::addDnbExp(vector<cv::Point> & dnb_coordinates,
                           map<unsigned long long int, pair<unsigned int, unsigned short>> & bin_gene_exp_map,
                           const DnbExpression *dnb_expression,
                           const cv::Point& center_point,
                           unsigned short area) {
    unsigned long long int bin_id;
    map<unsigned int, unsigned short> gene_count_in_cell;
    unsigned short gene_count = 0;
    unsigned short exp_count = 0;

    for (auto & dnb_coordinate : dnb_coordinates) {
        bin_id = static_cast<unsigned long long int>(dnb_coordinate.x);
        bin_id = bin_id << 32 | static_cast<unsigned int>(dnb_coordinate.y);

        auto iter = bin_gene_exp_map.find(bin_id);
        if(iter != bin_gene_exp_map.end()){
            pair<unsigned int, unsigned short> gene_info = iter->second;
            unsigned int end = gene_info.first + gene_info.second;
            for(unsigned int i = gene_info.first; i < end; i++){
                exp_count += dnb_expression[i].count;
                auto iter_gene = gene_count_in_cell.find(dnb_expression[i].gene_id);
                if(iter_gene != gene_count_in_cell.end()){
                    iter_gene->second += dnb_expression[i].count;
                } else{
                    gene_count_in_cell.insert(
                            map<unsigned int, unsigned short>::value_type(
                                    dnb_expression[i].gene_id, dnb_expression[i].count));
                    gene_count++;
                }
            }
        }
    }

    unsigned short cell_type_id = random_cell_type_num_ == 0 ? 0 : rand()%(random_cell_type_num_ + 1);

    CellData cell = {
            cell_num_,
            center_point.x,
            center_point.y,
            expression_num_, //offset
            gene_count,
            exp_count,
            dnb_coordinates.size(),
            area,
            cell_type_id
    };
    expression_num_ += gene_count;

    cell_attr_.min_x = cell.x < cell_attr_.min_x ? cell.x : cell_attr_.min_x;
    cell_attr_.max_x = cell.x > cell_attr_.max_x ? cell.x : cell_attr_.max_x;
    cell_attr_.min_y = cell.y < cell_attr_.min_y ? cell.y : cell_attr_.min_y;
    cell_attr_.max_y = cell.y > cell_attr_.max_y ? cell.y : cell_attr_.max_y;

    cell_attr_.min_area = area < cell_attr_.min_area ? area : cell_attr_.min_area;
    cell_attr_.max_area = area > cell_attr_.max_area ? area : cell_attr_.max_area;
    cell_attr_.min_gene_count = gene_count < cell_attr_.min_gene_count ? gene_count : cell_attr_.min_gene_count;
    cell_attr_.max_gene_count = gene_count > cell_attr_.max_gene_count ? gene_count : cell_attr_.max_gene_count;
    cell_attr_.min_exp_count = exp_count < cell_attr_.min_exp_count ? exp_count : cell_attr_.min_exp_count;
    cell_attr_.max_exp_count = exp_count > cell_attr_.max_exp_count ? exp_count : cell_attr_.max_exp_count;

    unsigned int dnb_count = dnb_coordinates.size();
    cell_attr_.min_dnb_count = dnb_count < cell_attr_.min_dnb_count ? dnb_count : cell_attr_.min_dnb_count;
    cell_attr_.max_dnb_count = dnb_count > cell_attr_.max_dnb_count ? dnb_count : cell_attr_.max_dnb_count;

    exp_count_sum_ += exp_count;
    dnb_count_sum_ += dnb_count;
    area_sum_ += area;

    cell_list_.emplace_back(cell);

    map<unsigned int, unsigned short> ::iterator iter_m;
    iter_m = gene_count_in_cell.begin();
    while(iter_m != gene_count_in_cell.end()) {
        unsigned int gene_id = iter_m->first;
        unsigned short count = iter_m->second;
        max_mid_count_ = count > max_mid_count_ ? count : max_mid_count_;

        // 用于生成geneExp dataset的数据
        GeneExpData gene_exp_tmp = {cell_num_, count};
        auto iter_gene_exp_map = gene_exp_map_.find(gene_id);
        if(iter_gene_exp_map != gene_exp_map_.end()){
            iter_gene_exp_map->second.emplace_back(gene_exp_tmp);
        } else{
            vector<GeneExpData> v_gene_exp_data;
            v_gene_exp_data.emplace_back(gene_exp_tmp);
            gene_exp_map_.insert(map<unsigned int, vector<GeneExpData>>::value_type(gene_id, v_gene_exp_data));
        }

        CellExpData cexp_tmp = {gene_id, count};
        cell_exp_list_.emplace_back(cexp_tmp);
        ++iter_m;
    }

    cell_num_ += 1;
}


void CgefWriter::storeAttr(CellBinAttr & cell_bin_attr) const {
    unsigned long cprev=clock();
    hsize_t dimsAttr[1] = {1};
    hid_t attr;
    hid_t attr_dataspace = H5Screate_simple(1, dimsAttr, nullptr);
    attr = H5Acreate(file_id_, "version", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT32, &cell_bin_attr.version);
    attr = H5Acreate(file_id_, "resolution", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT32, &cell_bin_attr.resolution);
    attr = H5Acreate(file_id_, "offsetX", H5T_STD_I32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT32, &cell_bin_attr.offsetX);
    attr = H5Acreate(file_id_, "offsetY", H5T_STD_I32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT32, &cell_bin_attr.offsetY);

    //Write createTime into cell bin gef
//    S32 time_str = getStrfTime();
//    attr = H5Acreate(file_id_, "createTime", str32_type_, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
//    H5Awrite(attr, str32_type_, &time_str);

    H5Aclose(attr);
    H5Sclose(attr_dataspace);

    hsize_t gef_dimsAttr[1] = {3};
    hid_t gef_dataspace_id = H5Screate_simple(1, gef_dimsAttr, nullptr);
    hid_t gef_attr = H5Acreate(file_id_, "geftool_ver", H5T_STD_U32LE, gef_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(gef_attr, H5T_NATIVE_UINT32, GEFVERSION);
    H5Sclose(gef_dataspace_id);
    H5Aclose(gef_attr);

    hsize_t kind_dims[1] = {1};
    hid_t k_did = H5Screate_simple(1, kind_dims, nullptr);
    hid_t k_attr = H5Acreate(file_id_, "omics", str32_type_, k_did, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(k_attr, str32_type_, cell_bin_attr.omics.c_str());
    H5Sclose(k_did);
    H5Aclose(k_attr);
    
    if(verbose_) printCpuTime(cprev, "storeAttr");
}

void CgefWriter::createGenedata(const vector<string> &gene_name_list)
{
    unsigned long cprev=clock();
    gene_num_ = gene_name_list.size();
    GeneData* gene_data_list;
    gene_data_list = static_cast<GeneData *>(malloc(gene_num_ * sizeof(GeneData)));

    unsigned int exp_count, min_exp_count = UINT32_MAX, max_exp_count = 0, offset = 0;
    unsigned int cell_count, min_cell_count = UINT32_MAX, max_cell_count = 0;
    unsigned short max_MID_count;

    vector<GeneExpData> gene_exp_list;
    gene_exp_list.reserve(expression_num_);
    for(unsigned int i = 0; i < gene_num_; i++){
        auto iter_gene_exp_map = gene_exp_map_.find(i);
        if(iter_gene_exp_map != gene_exp_map_.end()){
            vector<GeneExpData> tmp = iter_gene_exp_map->second;
            gene_exp_list.insert(gene_exp_list.end(), tmp.begin(), tmp.end());

            cell_count = static_cast<unsigned int>(tmp.size());
            max_MID_count = 0;
            exp_count = 0;
            for (auto gene_exp :tmp) {
                exp_count += gene_exp.count;
                max_MID_count = gene_exp.count > max_MID_count ? gene_exp.count : max_MID_count;
            }

            min_exp_count = min_exp_count < exp_count ? min_exp_count : exp_count;
            max_exp_count = max_exp_count > exp_count ? max_exp_count : exp_count;
            min_cell_count = min_cell_count < cell_count ? min_cell_count : cell_count;
            max_cell_count = max_cell_count > cell_count ? max_cell_count : cell_count;

            gene_data_list[i] = GeneData(
                    gene_name_list[i].c_str(),
                    offset,
                    static_cast<unsigned int>(tmp.size()),
                    exp_count,
                    max_MID_count);

            offset += tmp.size();
        }else{
            gene_data_list[i] = GeneData(
                    gene_name_list[i].c_str(),
                    offset,
                    0,
                    0,
                    0);
        }
    }
    storeGeneAndGeneExp(min_exp_count, max_exp_count, min_cell_count, max_cell_count, gene_data_list, gene_exp_list);
    free(gene_data_list);
    if(verbose_) printCpuTime(cprev, "createGenedata");
}

void CgefWriter::storeGeneAndGeneExp(unsigned int min_exp_count, unsigned int max_exp_count,
                                    unsigned int min_cell_count, unsigned int max_cell_count,
                                    GeneData* gene_data_list, vector<GeneExpData> &gene_exp_list) {
    hsize_t dims[1] = {gene_num_};

    // GeneData* gene_data_list;
    // gene_data_list = static_cast<GeneData *>(malloc(gene_num_ * sizeof(GeneData)));

    // unsigned int exp_count, min_exp_count = UINT32_MAX, max_exp_count = 0, offset = 0;
    // unsigned int cell_count, min_cell_count = UINT32_MAX, max_cell_count = 0;
    // unsigned short max_MID_count;

    // vector<GeneExpData> gene_exp_list;
    // gene_exp_list.reserve(expression_num_);
    // for(unsigned int i = 0; i < gene_num_; i++){
    //     auto iter_gene_exp_map = gene_exp_map_.find(i);
    //     if(iter_gene_exp_map != gene_exp_map_.end()){
    //         vector<GeneExpData> tmp = iter_gene_exp_map->second;
    //         gene_exp_list.insert(gene_exp_list.end(), tmp.begin(), tmp.end());

    //         cell_count = static_cast<unsigned int>(tmp.size());
    //         max_MID_count = 0;
    //         exp_count = 0;
    //         for (auto gene_exp :tmp) {
    //             exp_count += gene_exp.count;
    //             max_MID_count = gene_exp.count > max_MID_count ? gene_exp.count : max_MID_count;
    //         }

    //         min_exp_count = min_exp_count < exp_count ? min_exp_count : exp_count;
    //         max_exp_count = max_exp_count > exp_count ? max_exp_count : exp_count;
    //         min_cell_count = min_cell_count < cell_count ? min_cell_count : cell_count;
    //         max_cell_count = max_cell_count > cell_count ? max_cell_count : cell_count;

    //         gene_data_list[i] = GeneData(
    //                 gene_name_list[i].c_str(),
    //                 offset,
    //                 static_cast<unsigned int>(tmp.size()),
    //                 exp_count,
    //                 max_MID_count);

    //         offset += tmp.size();
    //     }else{
    //         gene_data_list[i] = GeneData(
    //                 gene_name_list[i].c_str(),
    //                 offset,
    //                 0,
    //                 0,
    //                 0);
    //     }
    // }

    hid_t memtype, filetype;
    memtype = getMemtypeOfGeneData();
    filetype = H5Tcreate(H5T_COMPOUND, 78);
    H5Tinsert(filetype, "geneName", 0, str64_type_);
    H5Tinsert(filetype, "offset", 64, H5T_STD_U32LE);
    H5Tinsert(filetype, "cellCount", 68, H5T_STD_U32LE);
    H5Tinsert(filetype, "expCount", 72, H5T_STD_U32LE);
    H5Tinsert(filetype, "maxMIDcount", 76, H5T_STD_U16LE);

    hid_t dataspace_id = H5Screate_simple(1, dims, nullptr);
    hid_t dataset_id = H5Dcreate(group_id_, "gene", filetype, dataspace_id, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, gene_data_list);

    hsize_t dims_attr[1] = {1};
    hid_t attr;
    hid_t attr_dataspace = H5Screate_simple(1, dims_attr, nullptr);
    attr = H5Acreate(dataset_id, "minExpCount", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT32, &min_exp_count);
    attr = H5Acreate(dataset_id, "maxExpCount", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT32, &max_exp_count);
    attr = H5Acreate(dataset_id, "minCellCount", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT32, &min_cell_count);
    attr = H5Acreate(dataset_id, "maxCellCount", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT32, &max_cell_count);

    memtype = getMemtypeOfGeneExpData();
    filetype = H5Tcreate(H5T_COMPOUND, 6);
    H5Tinsert(filetype, "cellID", 0, H5T_STD_U32LE);
    H5Tinsert(filetype, "count", 4, H5T_STD_U16LE);
    //H5Tinsert(filetype, "incnt", 6, H5T_STD_U16LE);

    hsize_t dims_exp[1] = {expression_num_};
    dataspace_id = H5Screate_simple(1, dims_exp, nullptr);
    dataset_id = H5Dcreate(group_id_, "geneExp", filetype, dataspace_id, H5P_DEFAULT,
                           H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &gene_exp_list[0]);

    attr = H5Acreate(dataset_id, "maxCount", H5T_STD_U16LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_USHORT, &max_mid_count_);

    H5Aclose(attr);
    H5Sclose(attr_dataspace);
    H5Tclose(memtype);
    H5Tclose(filetype);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
}

void CgefWriter::storeCellTypeList() {
    unsigned long cprev=clock();
    S32 cell_type = S32("default");
    cell_type_list_.emplace_back(cell_type);

    int i = 0;
    while(i < random_cell_type_num_){
        i++;
        cell_type = S32();
        sprintf(cell_type.value, "type%d", i);
        cell_type_list_.emplace_back(cell_type);
    }

    hsize_t dims[1];
    dims[0] = random_cell_type_num_ + 1;
    hid_t dataspace_id = H5Screate_simple(1, dims, nullptr);
    hid_t dataset_id = H5Dcreate(group_id_, "cellTypeList", str32_type_, dataspace_id, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, str32_type_, H5S_ALL, H5S_ALL, H5P_DEFAULT, &cell_type_list_[0]);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    if(verbose_) printCpuTime(cprev, "storeCellTypeList");
}

void CgefWriter::storeCellTypeList_N() {
    unsigned long cprev=clock();

    hsize_t dims[1];
    dims[0] = cell_type_list_.size();
    hid_t dataspace_id = H5Screate_simple(1, dims, nullptr);
    hid_t dataset_id = H5Dcreate(group_id_, "cellTypeList", str32_type_, dataspace_id, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, str32_type_, H5S_ALL, H5S_ALL, H5P_DEFAULT, &cell_type_list_[0]);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    if(verbose_) printCpuTime(cprev, "storeCellTypeList");
}

unsigned short CgefWriter::calcMaxCountOfGeneExp(vector<GeneExpData> &gene_exps) {
    unsigned max = 0;
    for (auto gene_exp :gene_exps) {
        max = gene_exp.count > max ? gene_exp.count : max;
    }
    return max;
}

int CgefWriter::write(BgefReader &common_bin_gef, Mask &mask) {
    map<unsigned long long int, pair<unsigned int, unsigned short>> bin_gene_exp_map;
    auto * dnb_exp_info = (DnbExpression *) malloc(common_bin_gef.getExpressionNum()  * sizeof(DnbExpression));
    common_bin_gef.getBinGeneExpMap(bin_gene_exp_map, dnb_exp_info);

    const vector<GefTools::Polygon>& polygons = mask.getPolygons();
// FILE *f = fopen("mask_raw.txt","wb");
// char pbuf[1024];

    unsigned long cprev=clock();
    for(unsigned int i = 0; i < mask.getCellNum(); i++){
        GefTools::Polygon p = polygons[i];
        cv::Rect roi = cv::Rect(p.getMinX(),p.getMinY(),p.getCols(),p.getRows());
        cv::Mat roi_mat = common_bin_gef.getWholeExpMatrix(roi);
        cv::Mat fill_points = p.getFillPolyMat();
        roi_mat = roi_mat.mul(fill_points);

        vector<cv::Point> non_zero_coordinates, non_zero_coordinates_offset;
        findNonZero(roi_mat,non_zero_coordinates);
        cv::Point offset = cv::Point(-p.getMinX(), -p.getMinY());
        offsetCoordinates(non_zero_coordinates, non_zero_coordinates_offset, offset);

        addDnbExp(
            non_zero_coordinates_offset,
            bin_gene_exp_map,
            dnb_exp_info,
            p.getCenter(),
            p.getAreaUshort());
        
        
        // const vector<Point> &vp = p.getBorder();
        // int l = 0;
        // for(const Point &po : vp)
        // {
        //     l += sprintf(pbuf+l,"%d-%d ",po.x, po.y);
        // }
        // l += sprintf(pbuf+l, "%d %d %d %d %d\n", vp.size(), p.getMinX(), p.getMinY(), p.getMaxX(), p.getMaxY());
        // fwrite(pbuf, 1, l, f);
    }
    //fclose(f);

    if(verbose_) printCpuTime(cprev, "addDnbExp");

    m_borderptr = static_cast<short *>(malloc(mask.getCellNum() * BORDERCNT * 2 * sizeof(short)));
    mask.getBorders(m_borderptr);

    ExpressionAttr expression_attr = common_bin_gef.getExpressionAttr();
    CellBinAttr cell_bin_attr = {
            /*.version = */1,
            /*.resolution = */expression_attr.resolution,
            /*.offsetX = */expression_attr.min_x,
            /*.offsetY = */expression_attr.min_y
    };

    storeAttr(cell_bin_attr);

    int effective_rect[4];
    mask.getEffectiveRectangle(effective_rect);
    storeCellBorderWithAttr(m_borderptr, mask.getCellNum(), effective_rect);
    storeCell(mask.getBlockNum(), mask.getBlockIndex(), mask.getBlockSize());
    storeCellExp();
    storeCellTypeList();

    vector<string> gene_name_list;
    gene_name_list.reserve(common_bin_gef.getGeneNum());
    common_bin_gef.getGeneNameList(gene_name_list);
    createGenedata(gene_name_list);

    free(dnb_exp_info);
    return 0;
}

unsigned short CgefWriter::getRandomCellTypeNum() const {
    return random_cell_type_num_;
}

void CgefWriter::setRandomCellTypeNum(unsigned short random_cell_type_num) {
    random_cell_type_num_ = random_cell_type_num;
}

bool CgefWriter::isVerbose() const {
    return verbose_;
}

void CgefWriter::setVerbose(bool verbose) {
    verbose_ = verbose;
}

int CgefWriter::addLevel_1()
{
    createBlktype();
    m_level_gid = H5Gcreate(group_id_, "level", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    int lev = 0;
    vector<int> vec_cellid;
    vector<block> vec_blk;
    vec_blk.emplace_back(0, cell_num_);

    for(int i=0;i<cell_num_;i++)
    {
        vec_cellid.emplace_back(i);
    }

    vector<int> vec_blk_idx;
    vec_blk_idx.emplace_back(0);
    int blknum[2]={1,1};
    writeCelldata(lev, blknum, vec_blk, vec_cellid, vec_blk_idx);
    lev++;
    hsize_t dims_attr[1] = {1};
    hid_t attr_dataspace = H5Screate_simple(1, dims_attr, nullptr);
    hid_t attr = H5Acreate(m_level_gid, "levelnum", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT, &lev);
    H5Aclose(attr);
    H5Sclose(attr_dataspace);
    
    H5Tclose(m_blk_memtype);
    H5Tclose(m_blk_filetype);
    H5Gclose(m_level_gid);
    return 0;
}

int CgefWriter::addLevel(int allocat, int cnum, float ratio, int *cansize, int *blknum)
{
    if(cansize[0] <= (m_canvas[0]+m_offsetX) &&
      cansize[2] >= (m_canvas[2]+m_offsetX) &&
      cansize[1] <= (m_canvas[1] + m_offsetY) && 
      cansize[3] >= (m_canvas[3] + m_offsetY))
    {
        m_canvas[0] = cansize[0];
        m_canvas[2] = cansize[2];
        printf("canvas ok\n");
    }
    else
    {
        printf("canvas too small\n");
        return 0;
    }

    m_x_len = cansize[2]-cansize[0];
    m_y_len = cansize[3]-cansize[1];
    
    m_blknum[0] = blknum[0];
    m_blknum[1] = blknum[1];

    m_allocat = allocat;
    createBlktype();
    m_level_gid = H5Gcreate(group_id_, "level", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for(int i=0;i<cell_num_;i++)
    {
        m_hash_cellid.emplace(i);
    }

    getblkcelldata_top(0, cnum);
    getblkcelldata(1, cnum);
    getblkcelldata(2, cnum);

    int cnt = 0, tmp = 0, lev = 3;
    long blkcnt = 0;
    while (1)
    {
        cnt = cell_num_*ratio;
        tmp = m_hash_cellid.size() - cnt;
        // blkcnt = pow(m_allocat,lev);
        // blkcnt*=blkcnt;
//printf("---%d %d %d\n", lev, cnt, tmp);
        if(tmp < 1000 || tmp < blkcnt) //最后一层数据太少，不再分层
        {
            getblkcelldata_bottom(lev);
            lev++;
            break;
        }
        else
        {
            getblkcelldata(lev, cnt);
        }
        lev++;
    }

    hsize_t dims_attr[1] = {1};
    hid_t attr_dataspace = H5Screate_simple(1, dims_attr, nullptr);
    hid_t attr = H5Acreate(m_level_gid, "levelnum", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT, &lev);
    H5Aclose(attr);
    H5Sclose(attr_dataspace);

    dims_attr[0] = 4;
    hid_t attr_ds = H5Screate_simple(1, dims_attr, nullptr);
    hid_t attr_di = H5Acreate(m_level_gid, "canvas", H5T_STD_I32LE, attr_ds, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr_di, H5T_NATIVE_INT, cansize);
    H5Sclose(attr_ds);
    H5Aclose(attr_di);
    
    
    H5Tclose(m_blk_memtype);
    H5Tclose(m_blk_filetype);
    H5Gclose(m_level_gid);
}

void CgefWriter::getblkcelldata_top(int lev, int cnt)
{
    vector<int> vec_cellid;
    vector<block> vec_blk;
    vec_blk.emplace_back(0, cnt);
    default_random_engine rand(time(NULL));
    uniform_int_distribution<int> rand1(0, m_hash_cellid.size()-1);
    set<int> set_tmp;
    int idx = 0;
    while (1)
    {
        idx = rand1(rand);
        if(set_tmp.emplace(idx).second)
        {
            vec_cellid.push_back(idx);
            m_hash_cellid.erase(idx); //记录被筛掉的cellid
        }
        if(set_tmp.size() >= cnt || m_hash_cellid.size() == 0)
        {
            break;
        }
    }

    vector<int> vec_blk_idx;
    vec_blk_idx.emplace_back(0);
    int blknum[2]={1,1};
    writeCelldata(lev, blknum, vec_blk, vec_cellid, vec_blk_idx);
}

void CgefWriter::getblkcelldata_bottom(int lev)
{
    if(m_hash_cellid.empty()) return;
    int x_num = pow(m_allocat,lev);//1<<lev;
    int y_num = x_num;
    if(x_num > m_blknum[0]) x_num = m_blknum[0];
    if(y_num > m_blknum[1]) y_num = m_blknum[1];
    int x_size = ceil(m_x_len*1.0/x_num);
    int y_size = ceil(m_y_len*1.0/y_num);

    vector<vector<int>> vec_vec_cellid;
    for(int i=0;i<x_num*y_num;i++)
    {
        vector<int> vectmp;
        vec_vec_cellid.emplace_back(std::move(vectmp));
    }
    unsigned int id = 0;
    for(auto itor = m_hash_cellid.begin();itor != m_hash_cellid.end();itor++)
    {
        CellData &cdata = m_cdataPtr[*itor];
        id = ((cdata.x+m_offsetX-m_canvas[0]) / x_size) + ((cdata.y+m_offsetY-m_canvas[2]) / y_size) *y_num;
        vec_vec_cellid[id].emplace_back(*itor);
    }

    vector<int> vec_blk_idx; //非空的块下标
    vector<int> vec_cellid;
    vector<block> vec_blk;
    int offset = 0, cnt = 0;
    for(int i=0;i<x_num*y_num;i++)
    {
        vector<int> &blkcellid = vec_vec_cellid[i];
        cnt = blkcellid.size();
        vec_blk.emplace_back(offset, cnt);
        offset += cnt;
        if(cnt)
        {
            vec_blk_idx.emplace_back(i);
        }
        vec_cellid.insert(vec_cellid.end(), blkcellid.begin(), blkcellid.end());
    }
    int blknum[2]={x_num, y_num};
    writeCelldata(lev, blknum, vec_blk, vec_cellid, vec_blk_idx);
}

void CgefWriter::getblkcelldata(int lev, int cnt)
{
    if(m_hash_cellid.empty()) return;
    int x_num = pow(m_allocat,lev);
    int y_num = x_num;
    if(x_num > m_blknum[0]) x_num = m_blknum[0];
    if(y_num > m_blknum[1]) y_num = m_blknum[1];
    int x_size = ceil(m_x_len*1.0/x_num);
    int y_size = ceil(m_y_len*1.0/y_num);

    vector<vector<int>> vec_vec_cellid;
    for(int i=0;i<x_num*y_num;i++)
    {
        vector<int> vectmp;
        vec_vec_cellid.emplace_back(std::move(vectmp));
    }
    unsigned int id = 0;
    for(auto itor = m_hash_cellid.begin();itor != m_hash_cellid.end();itor++)
    {
        CellData &cdata = m_cdataPtr[*itor];
        id = ((cdata.x+m_offsetX-m_canvas[0]) / x_size) + ((cdata.y+m_offsetY-m_canvas[2]) / y_size) *y_num;
        vec_vec_cellid[id].emplace_back(*itor);
    }

    vector<int> vec_blk_idx; //非空的块下标
    vector<int> vec_cellid;
    vector<block> vec_blk;
    int offset = 0;
    int idx = 0;
    int scnt = 0;
    int left = m_hash_cellid.size();
    for(int i=0;i<x_num*y_num;i++)
    {
        vector<int> &blkcellid = vec_vec_cellid[i];
        scnt = blkcellid.size()*cnt/left;
        default_random_engine rand(time(NULL));
        uniform_int_distribution<int> rand1(0, blkcellid.size()-1);
        vec_blk.emplace_back(offset, scnt);
        offset += scnt;
        set<int> set_tmp;

        if(scnt)
        {
            vec_blk_idx.emplace_back(i);
        }

        while(scnt)
        {
            idx = rand1(rand);
            if(set_tmp.emplace(idx).second)
            {
                vec_cellid.push_back(blkcellid[idx]);
                m_hash_cellid.erase(blkcellid[idx]); 
            }
            if(set_tmp.size() >= scnt)
            {
                break;
            }
        }
    }

    int blknum[2]={x_num, y_num};
    writeCelldata(lev, blknum, vec_blk, vec_cellid, vec_blk_idx);
}

void CgefWriter::createBlktype()
{
    m_blk_memtype = H5Tcreate(H5T_COMPOUND, sizeof(block));
    H5Tinsert(m_blk_memtype, "offset", HOFFSET(block, offset), H5T_NATIVE_UINT32);
    H5Tinsert(m_blk_memtype, "count", HOFFSET(block, count), H5T_NATIVE_UINT32);

    m_blk_filetype = H5Tcreate(H5T_COMPOUND, 8);
    H5Tinsert(m_blk_filetype, "offset", 0, H5T_STD_U32LE);
    H5Tinsert(m_blk_filetype, "count", 4, H5T_STD_U32LE);
}

void CgefWriter::writeCelldata(int lev, int *blknum, vector<block> &blk, vector<int> &vecid, vector<int> &vec_blk_idx)
{
    printf("%d %d %d\n", lev, vecid.size(), blk.size());

    char buf[32]={0};
    sprintf(buf, "L%d", lev);
    hid_t level_gid = H5Gcreate(m_level_gid, buf, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    hsize_t dimsAttr[1] = {2};
    hid_t attr_dataspace = H5Screate_simple(1, dimsAttr, nullptr);
    hid_t attr = H5Acreate(level_gid, "blknum", H5T_STD_U32LE, attr_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT32, blknum);
    H5Sclose(attr_dataspace);
    H5Aclose(attr);

    hsize_t blk_dims[1] = {blk.size()};
    hid_t blk_dspace_id = H5Screate_simple(1, blk_dims, nullptr);
    hid_t blk_dset_id = H5Dcreate(level_gid, "blk", m_blk_memtype, blk_dspace_id, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(blk_dset_id, m_blk_filetype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &blk[0]);

    H5Sclose(blk_dspace_id);
    H5Dclose(blk_dset_id);

    hsize_t id_dims[1] = {vecid.size()};
    hid_t id_dspace_id = H5Screate_simple(1, id_dims, nullptr);
    hid_t id_dset_id = H5Dcreate(level_gid, "cellid", H5T_NATIVE_UINT32, id_dspace_id, H5P_DEFAULT,
                                 H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(id_dset_id, H5T_STD_U32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vecid[0]);

    H5Sclose(id_dspace_id);
    H5Dclose(id_dset_id);

    hsize_t bid_dims[1] = {vec_blk_idx.size()};
    hid_t bid_dspace_id = H5Screate_simple(1, bid_dims, nullptr);
    hid_t bid_dset_id = H5Dcreate(level_gid, "noempty", H5T_NATIVE_UINT32, bid_dspace_id, H5P_DEFAULT,
                                H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(bid_dset_id, H5T_STD_U32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vec_blk_idx[0]);

    H5Sclose(bid_dspace_id);
    H5Dclose(bid_dset_id);
    H5Gclose(level_gid);
}

void CgefWriter::storeGeneExon(uint32_t minExon, uint32_t maxExon, uint32_t *geneExonPtr, uint16_t maxExpExon, vector<uint16_t> vec_exonExp)
{
    hsize_t dims[1] = {gene_num_};
    hid_t exon_sid = H5Screate_simple(1, dims, nullptr);
    hid_t exon_did = H5Dcreate(group_id_, "geneExon", H5T_STD_U32LE, exon_sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(exon_did, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, geneExonPtr);

    hsize_t dims_attr[1] = {1};
    hid_t attr_sid = H5Screate_simple(1, dims_attr, nullptr);
    hid_t attr = H5Acreate(exon_did, "minExon", H5T_STD_U32LE, attr_sid, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT32, &minExon);
    attr = H5Acreate(exon_did, "maxExon", H5T_STD_U32LE, attr_sid, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_UINT32, &maxExon);

    H5Aclose(attr);
    H5Sclose(exon_sid);
    H5Dclose(exon_did);

    dims[0] = vec_exonExp.size();
    hid_t exonExp_sid = H5Screate_simple(1, dims, nullptr);
    hid_t exonExp_did = H5Dcreate(group_id_, "geneExpExon", H5T_STD_U16LE, exonExp_sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(exonExp_did, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec_exonExp.data());

    attr = H5Acreate(exonExp_did, "maxExon", H5T_STD_U16LE, attr_sid, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_USHORT, &maxExpExon);
 
    H5Aclose(attr);
    H5Sclose(attr_sid);
    H5Sclose(exonExp_sid);
    H5Dclose(exonExp_did);
}

void CgefWriter::storeCellExon(uint16_t minExon, uint16_t maxExon, vector<uint16_t> vec_cellexon, uint16_t maxExpExon, vector<uint16_t> vec_cellexon_exp)
{
    hsize_t dims[1] = {(hsize_t)cell_num_};
    hid_t exon_sid = H5Screate_simple(1, dims, nullptr);
    hid_t exon_did = H5Dcreate(group_id_, "cellExon", H5T_STD_U16LE, exon_sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(exon_did, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec_cellexon.data());

    hsize_t dims_attr[1] = {1};
    hid_t attr_sid = H5Screate_simple(1, dims_attr, nullptr);
    hid_t attr = H5Acreate(exon_did, "minExon", H5T_STD_U16LE, attr_sid, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_USHORT, &minExon);
    attr = H5Acreate(exon_did, "maxExon", H5T_STD_U16LE, attr_sid, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_USHORT, &maxExon);

    H5Aclose(attr);
    H5Sclose(exon_sid);
    H5Dclose(exon_did);

    dims[0] = vec_cellexon_exp.size();
    hid_t exonExp_sid = H5Screate_simple(1, dims, nullptr);
    hid_t exonExp_did = H5Dcreate(group_id_, "cellExpExon", H5T_STD_U16LE, exonExp_sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(exonExp_did, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec_cellexon_exp.data());

    attr = H5Acreate(exonExp_did, "maxExon", H5T_STD_U16LE, attr_sid, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_USHORT, &maxExpExon);
 
    H5Aclose(attr);
    H5Sclose(attr_sid);
    H5Sclose(exonExp_sid);
    H5Dclose(exonExp_did);
}