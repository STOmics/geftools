/*
 * @Author: zhaozijian
 * @Date: 2022-02-10 14:53:03
 * @LastEditors: zhaozijian
 * @LastEditTime: 2022-04-29 15:35:32
 * @Description: file content
 */
#include "gef.h"
#include "utils.h"

hid_t getMemtypeOfGeneData(){
    hid_t memtype;
    // hid_t str32_type_ = H5Tcopy(H5T_C_S1);
    // H5Tset_size(str32_type_, 32);
    hid_t str64_type_ = H5Tcopy(H5T_C_S1);
    H5Tset_size(str64_type_, 64);
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(GeneData));
    H5Tinsert(memtype, "geneName", HOFFSET(GeneData, gene_name), str64_type_);
    H5Tinsert(memtype, "offset", HOFFSET(GeneData, offset), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "cellCount", HOFFSET(GeneData, cell_count), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "expCount", HOFFSET(GeneData, exp_count), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "maxMIDcount", HOFFSET(GeneData, max_mid_count), H5T_NATIVE_USHORT);
    return memtype;
}

hid_t getMemtypeOfGeneExpData() {
    hid_t memtype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(GeneExpData));
    H5Tinsert(memtype, "cellID", HOFFSET(GeneExpData, cell_id), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "count", HOFFSET(GeneExpData, count), H5T_NATIVE_USHORT);
    //H5Tinsert(memtype, "incnt", HOFFSET(GeneExpData, incnt), H5T_NATIVE_USHORT);
    return memtype;
}


hid_t getMemtypeOfCellData() {
    hid_t memtype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(CellData));
    H5Tinsert(memtype, "id", HOFFSET(CellData, id), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "x", HOFFSET(CellData, x), H5T_NATIVE_INT);
    H5Tinsert(memtype, "y", HOFFSET(CellData, y), H5T_NATIVE_INT);
    H5Tinsert(memtype, "offset", HOFFSET(CellData, offset), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "geneCount", HOFFSET(CellData, gene_count), H5T_NATIVE_USHORT);
    H5Tinsert(memtype, "expCount", HOFFSET(CellData, exp_count), H5T_NATIVE_USHORT);
    H5Tinsert(memtype, "dnbCount", HOFFSET(CellData, dnb_count), H5T_NATIVE_USHORT);
    H5Tinsert(memtype, "area", HOFFSET(CellData, area), H5T_NATIVE_USHORT);
    H5Tinsert(memtype, "cellTypeID", HOFFSET(CellData, cell_type_id), H5T_NATIVE_USHORT);
    H5Tinsert(memtype, "clusterID", HOFFSET(CellData, cluster_id), H5T_NATIVE_USHORT);
    //H5Tinsert(memtype, "incnt", HOFFSET(CellData, incnt), H5T_NATIVE_USHORT);
    
    return memtype;
}

hid_t getMemtypeOfCellExpData() {
    hid_t memtype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(CellExpData));
    H5Tinsert(memtype, "geneID", HOFFSET(CellExpData, gene_id), H5T_NATIVE_UINT);
    H5Tinsert(memtype, "count", HOFFSET(CellExpData, count), H5T_NATIVE_USHORT);
    return memtype;
}

hid_t getMemtypeOfOlderCellExpData() {
    hid_t memtype;
    memtype = H5Tcreate(H5T_COMPOUND, sizeof(olderCellExpData));
    H5Tinsert(memtype, "geneID", HOFFSET(olderCellExpData, gene_id), H5T_NATIVE_USHORT);
    H5Tinsert(memtype, "count", HOFFSET(olderCellExpData, count), H5T_NATIVE_USHORT);
    return memtype;
}

bool isOlderCellExpDataVersion(hid_t fileId)
{
    unsigned int geftoolVersion[3] ={0};
    hid_t attr = H5Aopen(fileId, "geftool_ver", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_UINT32, geftoolVersion);
    printf("version is %d.%d.%d ", geftoolVersion[0], geftoolVersion[1], geftoolVersion[2]);
    
    if (geftoolVersion[0] == 0 && geftoolVersion[1] <= 7) {
        if (geftoolVersion[1] == 7 && geftoolVersion[2] >= 6) {
            return false;
        } else {
            return true;
        }
    } else {
        return false;
    }
}
