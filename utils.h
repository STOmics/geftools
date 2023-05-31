/** @file utils.h
    @brief Collection of utility functions.

    Created by huangzhibo on 2021/12/17.
*/

#ifndef GEFTOOLS__UTILS_H_
#define GEFTOOLS__UTILS_H_

#include <string>
#include <ctime>
#include <iostream>
#include <fstream>
#include <array>
#include <unordered_map>
#include <iomanip>
#include <zlib.h>
#include "hdf5.h"
#include "opencv2/opencv.hpp"
#include "tiffio.h"

#ifdef _WIN32
    #define GEFTOOLS_API __declspec(dllexport)
#else 
    #define GEFTOOLS_API
#endif

using namespace std;
//using namespace cv;

const int READLEN = 256*1024;
const unsigned int GEFVERSION[3] = {0,7,18};
const int BORDERCNT = 32;

namespace errorCode {
    extern bool isInSAWFlow;
    const static char* E_INVALIDPARAM        = "SAW-A60001";
    const static char* E_FILEOPENERROR       = "SAW-A60002";
    const static char* E_PARSEFILEERROR      = "SAW-A60003";
    const static char* E_LOWVERSION          = "SAW-A60110";
    const static char* E_STEPERROR           = "SAW-A60111";
    const static char* E_FILEDATAERROR       = "SAW-A60120";
    const static char* E_MISSINGFILEINFO     = "SAW-A60121";
    const static char* E_FILEMISMATCH        = "SAW-A60122";
    const static char* E_CREATEFILEFAILED    = "SAW-A60130";
    const static char* E_ALLOCMEMORYFAILED   = "SAW-A60140";
    const static char* E_GENEEXPDIMDISMATCH  = "SAW-A60150";
}

void reportErrorCode2File(const char* errCode, const char* errMsg);

static union
{
    char c[4];
    unsigned long l;
}endian_test = { { 'l','?','?','b' } };
#define ENDIANNESS ((char)endian_test.l)


/**
 * \brief Define a String type with 32 byte
 *
 * It is compatible with HDF5 H5T_C_S1 type.
 *   strtype = H5Tcopy(H5T_C_S1);
 *   H5Tset_size(strtype, 32);
 */
struct S32 {
    S32() = default;

    explicit S32(const char *c) {
        int i = 0;
        while (c[i] != '\0') {
            value[i] = c[i];
            ++i;
        }
    }

    char value[32] = {0};
};

/**
 * @brief Split string based on a character delimiter.
 * @param s   The string to split.
 * @param delim  The delimiter, a character.
 * @return
 */
vector<string> split(const string &s, char delim);

/**
 * @brief Read lines from txt file.
 * @param filename
 * @return
 */
vector<string> readLines(const string &filename);

/**
 * @brief Print time.
 * @param prev   Previous time， ( the start time ).
 * @param message  Message of the step.
 * @return Current time, used as the next previous time.
 */
time_t printTime(time_t prev, string message);

/**
 * @brief Print cpu time.
 * @param prev   Previous clock time， ( the start time ).
 * @param message  Message of the step.
 * @return Current clock time, used as the next previous clock time.
 */
unsigned long printCpuTime(unsigned long prev, string message);

/**
 * @brief Gets the format string of now time.
 */
S32 getStrfTime();

/*!
 * \brief Copies one file to another path
 * \param src_file The source file path
 * \param dst_file The destination file path
 * \return true if the file was copied successfully
 */
bool copyFile(const string& src_file, const string& dst_file);

/**
 *
 * @param coordinates
 * @param new_coordinates
 * @param offset_point
 */
void offsetCoordinates(vector<cv::Point> & coordinates, vector<cv::Point> & new_coordinates, cv::Point & offset_point);


bool readline(gzFile f, std::string& line);



#define SPTOOLS_CSR_DEFINE_TEMPLATE(I, T) \
  template void csr_tocsc(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], I Bp[], I Bi[], T Bx[]);

#define SPTOOLS_CSR_EXTERN_TEMPLATE(I, T) \
  extern template void csr_tocsc(const I n_row, const I n_col, const I Ap[], const I Aj[], const T Ax[], I Bp[], I Bi[], T Bx[]);

/*
 * Compute B = A for CSR matrix A, CSC matrix B
 *
 * Also, with the appropriate arguments can also be used to:
 *   - compute B = A^t for CSR matrix A, CSR matrix B
 *   - compute B = A^t for CSC matrix A, CSC matrix B
 *   - convert CSC->CSR
 *
 * Input Arguments:
 *   I  n_row         - number of rows in A
 *   I  n_col         - number of columns in A
 *   I  Ap[n_row+1]   - row pointer
 *   I  Aj[nnz(A)]    - column indices
 *   T  Ax[nnz(A)]    - nonzeros
 *
 * Output Arguments:
 *   I  Bp[n_col+1] - column pointer
 *   I  Bi[nnz(A)]  - row indices
 *   T  Bx[nnz(A)]  - nonzeros
 *
 * Note:
 *   Output arrays Bp, Bi, Bx must be preallocated
 *
 * Note:
 *   Input:  column indices *are not* assumed to be in sorted order
 *   Output: row indices *will be* in sorted order
 *
 *   Complexity: Linear.  Specifically O(nnz(A) + max(n_row,n_col))
 *
 */
template <class I, class T>
void csr_tocsc(const I n_row,
               const I n_col,
               const I Ap[],
               const I Aj[],
               const T Ax[],
               I Bp[],
               I Bi[],
               T Bx[])
{
    const I nnz = Ap[n_row];

    //compute number of non-zero entries per column of A
    std::fill(Bp, Bp + n_col, 0);

    for (I n = 0; n < nnz; n++){
        Bp[Aj[n]]++;
    }

    //cumsum the nnz per column to get Bp[]
    for(I col = 0, cumsum = 0; col < n_col; col++){
        I temp  = Bp[col];
        Bp[col] = cumsum;
        cumsum += temp;
    }
    Bp[n_col] = nnz;

    for(I row = 0; row < n_row; row++){
        for(I jj = Ap[row]; jj < Ap[row+1]; jj++){
            I col  = Aj[jj];
            I dest = Bp[col];

            Bi[dest] = row;
            Bx[dest] = Ax[jj];

            Bp[col]++;
        }
    }

    for(I col = 0, last = 0; col <= n_col; col++){
        I temp  = Bp[col];
        Bp[col] = last;
        last    = temp;
    }
}

bool decideSuffix(string& filename, string suffix);

/**
 * @brief Verifies that the target object exists
 * @param loc_id
 * @param name
 * @return
 */
bool utils_hdf5_check_present(hid_t loc_id, const char *name);

bool is_bgef(string& filename);
bool is_cgef(string& filename);

bool h5AttrWrite(hid_t id, hid_t filetype, hid_t memtype, const char *name, int rank, hsize_t *dims, const void* value);
hid_t h5DatasetWrite(hid_t gid, hid_t filetype, hid_t memtype, const char *name, int rank, hsize_t *dims, const void *data);

size_t Rect_hash(const cv::Rect &rect);
bool Rectequal_to(const cv::Rect &a, const cv::Rect &b);

unsigned int parseResolutin(const string& filename);

long tifread(cv::Mat &img, const string &strtif);
#endif //GEFTOOLS__UTILS_H_
