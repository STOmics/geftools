/** @file cgef_reader.h
    @brief Declare a CgefReader class for reading cell bin GEFs.

    Created by huangzhibo on 2021/12/27.
*/

#ifndef GEFTOOLS_CGEF_READER_H
#define GEFTOOLS_CGEF_READER_H

#include <numeric>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include "hdf5.h"
#include "opencv2/opencv.hpp"
#include "utils.h"
#include "gef.h"

struct GEFTOOLS_API Region {
    unsigned int min_x;
    unsigned int max_x;
    unsigned int min_y;
    unsigned int max_y;
};

class GEFTOOLS_API CgefReader {
  private:
    hid_t file_id_;
    hid_t group_id_;
    hid_t str32_type_;
    hid_t cell_dataset_id_;
    hid_t cell_dataspace_id_;
    hid_t cell_exp_dataset_id_;
    hid_t cell_exp_dataspace_id_;
    hid_t gene_dataset_id_;
    hid_t gene_exp_dataset_id_;
    hid_t gene_exp_dataspace_id_;

    unsigned int gene_num_ = 0;
    unsigned int gene_num_current_ = 0;
    GeneData* gene_array_ = nullptr;
    GeneData* gene_array_current_ = nullptr;
    //如果 gene_id_to_index_[gene_id] >= 0，则 gene 在restrict_gene_的区域内。
    int* gene_id_to_index_ = nullptr;

    unsigned int cell_num_ = 0;
    unsigned int cell_num_current_ = 0;
    CellData* cell_array_ = nullptr;
    CellData* cell_array_current_ = nullptr;
    unsigned int* cell_id_array_current_ = nullptr;
    //如果 isInRegion[cell_id - start_cell_id] == true，则 cell 在restrict_region_的区域内。
    int* cell_id_to_index_ = nullptr;
    unsigned int start_cell_id = 0;
    unsigned int end_cell_id = UINT_MAX;

    unsigned long long expression_num_ = 0;
    unsigned int expression_num_current_ = 0;

    unsigned int block_size_[4]{};  ///< x_block_size, y_block_size, x_block_num, y_block_num
    unsigned int* block_index_;  ///< offset, count
    unordered_map<string,unsigned int> genename_to_id_;

    hid_t openCellDataset(hid_t group_id);
    hid_t openCellExpDataset(hid_t group_id);
    hid_t openGeneDataset(hid_t group_id);
    hid_t openGeneExpDataset(hid_t group_id);

    GeneData *loadGene(bool reload = false);

    CellData *loadCell(bool reload = false);

    Region region_{};  ///< min_x, max_x, min_y, max_y
    bool verbose_ = false;
    bool restrict_region_ = false;
    bool restrict_gene_ = false;

    // char *m_borderdataPtr = nullptr;
    // char *m_borderdata_currentPtr = nullptr;
    short *m_borderdataPtr_s = nullptr;
    //short *m_borderdata_currentPtr_s = nullptr;
    short *m_pborcnt = nullptr;

    unsigned int m_ver = 0;
    unsigned int m_resolution = 0;
    int offsetX = 0;
    int offsetY = 0;
    unsigned int m_ver_tool[3] ={0};
    int m_bordercnt = 0;
    bool m_bexon = false;
    bool isOldCellExpVersion = false;
  public:
    explicit CgefReader(const string &filename, bool verbose = false);
    ~CgefReader();

    unsigned int getGeneNum() const;
    unsigned int getCellNum() const;
    unsigned long long int getExpressionNum() const;

    string getGeneName(unsigned int gene_id);

    /**
     * @brief Gets gene name array, 32 bit for each string.
     * @param gene_names
     */
    void getGeneNames(char * gene_names);

    int getGeneId(string& gene_name);
    GeneData *getGene();
    GeneData getGene(unsigned int gene_id) const;
    CellData *getCell();
    CellData getCell(unsigned int cell_id) const;

    /**
     * @brief Gets gene name array, 32 bit for each string.
     * @param gene_list
     */
    void getGeneNameList(vector<string> & gene_list);

    /**
     * @brief Gets cell pos array.  store x,y in a number (unsigned long long int) : x << 32 | y
     * @param cell_name_list
     */
    void getCellNameList(unsigned long long int * cell_name_list);

    /**
     * @brief Use blocks that intersect the input region.
     *
     * Some member variables (e.g. cell_num_) of this class will be updated.
     * @param min_x
     * @param max_x
     * @param min_y
     * @param max_y
     */
    void restrictRegion(unsigned int min_x, unsigned int max_x, unsigned int min_y, unsigned int max_y);

    /**
     * @brief Restrict to a gene list.
     * @param gene_list
     * @param exclude
     */
    void restrictGene(vector<string> & gene_list, bool exclude = false);

    void freeRestriction();

    /**
     * @brief Update current gene array and gene number.
     * Delete genes that do not appear in the current cell array.
     */
    void updateGeneInfo();

    /**
     * @brief Gets indices for building csr_matrix.
     *
     * Examples:
     * @code
     * # Python
     * from scipy import sparse
     * sparse.csr_matrix((data, indices, indptr), shape=(M, N))
     * @endcode
     * @param indices  CSR format index array of the matrix. Cell id or Gene id array, same size as count.
     * @param indptr   CSR format index pointer array of the matrix.
     * @param count    CSR format data array of the matrix. Expression count.
     * @param order    Order of count, "gene" or "cell".
     * @return
     */
    int getSparseMatrixIndices(unsigned int * indices, unsigned int * indptr, unsigned int * count, const char * order);

    /**
     * @brief Gets indices for building csr_matrix.
     *
     * @param cell_ind     CSR format index array of the matrix. same size as count.
     * @param gene_ind     CSR format index array of the matrix. same size as count.
     * @param count        CSR format data array of the matrix. Expression count.
     */
    int getSparseMatrixIndices2(unsigned int * cell_ind, unsigned int * gene_ind, unsigned int * count);

    /**
     * @brief Gets cellId and count from the geneExp dataset.
     * @param cell_id  Cell id.
     * @param count    Gene expression count.
     */
    void getCellIdAndCount(unsigned int *cell_id, unsigned short *count) const;

    /**
     * @brief Gets geneId and count from the cellExp dataset.
     * @param gene_id  Gene id.
     * @param count    Gene expression count.
     */
    void getGeneIdAndCount(unsigned int *gene_id, unsigned short *count) const;

    /**
     *
     * @param gene_name
     * @param expressions
     * @return
     */
    unsigned int getExpressionCountByGene(string& gene_name, GeneExpData* expressions);

    unsigned int getExpressionCountByGeneId(unsigned int gene_id, GeneExpData* expressions);

    unsigned int getCellCount(string& gene_name);

    unsigned int getCellCount(unsigned int gene_id);

    unsigned short getGeneCount(unsigned int cell_id) const;

    GeneData getGeneDataByGeneId(unsigned int gene_id);

    /**
     * @brief Output gene expression info to gem format.
     *
     * GEM format:
     * # geneName  x  y  MIDcount  cellID
     * @param filename  The output filename. If it's "stdout", the contents will be printed to standard output.
     * @param gene_name_list List of genes to include.
     * @param force_genes Only warn about unknown subset genes, default: false.
     * @param exclude The genes in gene_name_list is to exclude, default: false.
     * @return Number of output entries.
     */
    unsigned int toGem(string & filename,
                       const vector<string> & gene_name_list = vector<string>(),
                       bool force_genes = false,
                       bool exclude = false);

    /**
     * @brief Determine whether the restrictRegion function is run to limit to a rectangular region.
     */
    bool isRestrictRegion() const;

    /**
     * @brief Determine whether the restrictRene function is run to limit to a gene list.
     */
    bool isRestrictGene() const;

    bool isInRegion(unsigned int cell_id);

    bool isVerbose() const;

    void setVerbose(bool verbose);

    void selectCells(unsigned int offset, unsigned int cell_count, CellData *cell) const;

    void selectCellExp(unsigned int offset, unsigned int exp_count, CellExpData *cell_exp_data) const;

    void selectOlderCellExp(unsigned int offset, unsigned int exp_count, olderCellExpData *cell_exp_data) const;

    void selectGeneExp(unsigned int offset, unsigned int cell_count, GeneExpData *gene_exp_data) const;

    void closeH5();

    // char* getCellBorders_char(bool ball, unsigned int cell_id);
    // short* getCellBorders_short(bool ball, unsigned int cell_id);
    int getCellBorders(vector<unsigned int> &cell_ind, vector<short> &border);
    void getCellBorders(vector<unsigned int> &cell_ind, vector<short> &border, vector<short> &borcnt);
    void getAttr();
    uint32_t* getGefVer()
    {
        return m_ver_tool;
    }

    void getfiltereddata(vector<int> &region, vector<string> &genelist,
                    vector<string> &vec_gene, vector<unsigned long long> &uniq_cells,
                    vector<unsigned int> &cell_ind, vector<unsigned int> &gene_ind, vector<unsigned int> &count);

    bool isContainExon(){return m_bexon;}
    void getfiltereddata_exon(vector<int> &region, vector<string> &genelist,
                    vector<string> &vec_gene, vector<unsigned long long> &uniq_cells,
                    vector<unsigned int> &cell_ind, vector<unsigned int> &gene_ind, 
                    vector<unsigned int> &count, vector<unsigned int> &exon);
};

#endif //GEFTOOLS_CGEF_READER_H
