#include "utils.h"

bool errorCode::isInSAWFlow = false;

void reportErrorCode2File(const char* errCode, const char* errMsg) {
    if (!errorCode::isInSAWFlow) return;
    fstream fs;
	fs.open("errcode.log", ios::app);
    S32 time_str = getStrfTime();
	if (!fs) {
		ofstream fout("errcode.log");
		if (fout) {
			fout << "[" << time_str.value << "]" << " " << errCode << ": " << errMsg << endl;
			fout.close();
		}
	} else {
		fs << "[" << time_str.value << "]" << " " << errCode << ": " << errMsg << endl;
        fs.close();
	}
}

S32 getStrfTime()
{
    time_t timep;
    time (&timep);
    S32 tmp;
    strftime(tmp.value, sizeof(tmp), "%Y-%m-%d %H:%M:%S",localtime(&timep) );
    return tmp;
};

time_t printTime(time_t prev, string message){
    time_t cur;
    time(&cur);
    std::cout << std::setw (30) << message;
    double cost;
    cost = difftime(cur, prev);
    printf(" - %.f sec\n",cost);
    return cur;
}

unsigned long printCpuTime(unsigned long prev, string message){
    unsigned long cur=clock();
    std::cout << std::setw (30) << message;
    printf (" - %.6f cpu sec\n", static_cast<double>(cur - prev) / CLOCKS_PER_SEC);
    prev=cur;
    return prev;
}

vector<string> split(const string &s, char delim) {
    vector<string> result;
    stringstream ss (s);
    string item;

    while (getline (ss, item, delim)) {
        result.emplace_back(item);
    }

    return result;
}

vector<string> readLines(const string &filename){
    vector<string> result;
    char data[1000] = {0};
    ifstream infile;
    infile.open(filename);

    while(infile.getline(data,1000)){
        result.emplace_back(data);
    }

    if(!infile.eof()) {
        cerr << "Error to read file : " << filename << endl;
        char errMsg[32]={0};
        sprintf(errMsg, "Error to read file : %s", filename);
        reportErrorCode2File(errorCode::E_PARSEFILEERROR, errMsg);
        exit(2);
    }

    infile.close();
    return result;
}

bool copyFile(const string& src_file, const string& dst_file) {
    ifstream fin(src_file, ios::binary);
    ofstream fout(dst_file, ios::binary);

    bool ret = true;

    while(!fin.eof()){
        char buf;
        fin.read(&buf, sizeof(char));

        if(fin.eof()) break;

        if (fout.bad())
        {
            ret = false;
            break;
        }
        fout.write(&buf, sizeof(char));
    }

    fout.close();
    fin.close();
    return ret;
}

void offsetCoordinates(vector<cv::Point> & coordinates, vector<cv::Point> & new_coordinates, cv::Point & offset_point){
    for (auto & coordinate: coordinates) {
        cv::Point p = cv::Point(coordinate.x - offset_point.x, coordinate.y - offset_point.y);
        new_coordinates.emplace_back(p);
    }
}

bool readline(gzFile f, std::string& line)
{
    const int GZ_LINE_LEN = 1024;
    char GZ_LINE_BUFF[GZ_LINE_LEN];
    if (gzgets(f, GZ_LINE_BUFF, GZ_LINE_LEN) == Z_NULL)
    {
        // end-of-file or error
        int         err;
        const char* msg = gzerror(f, &err);
        if (err != Z_OK)
        {
            std::cerr << "read gz file error, error_code: " << err << " error_msg: " << msg << std::endl;
        }
        // cout<<"eof"<<endl;
        return false;
    }
    line.assign(GZ_LINE_BUFF);
    return true;
}

bool decideSuffix(string& filename, string suffix){
//    char dot = '.';
//    char suff[10] = {0};
//    int j = 0;
    size_t c1 = filename.size();
    size_t c2 = suffix.size();
    do{
        c1--;
        c2--;
        if(filename[c1] != suffix[c2])
            return false;
    }while (c2>0);
    return true;
//    for(int i = 0; i<c; i++)
//    {
//        if(gname[i] == dot)
//            j = i;
//    }
//    int k = j;
//    j = c - j - 1;
//    for(int i = 0; i<j; i++)
//    {
//        suffix[i] = gname[k+i+1];
//    }
//    if (0==strcmp(suff,nsuff))
//        return true;
//    else
//        return false;
}

bool utils_hdf5_check_present(hid_t loc_id, const char *name) {
    htri_t bool_id;
    if ((bool_id = H5Lexists(loc_id, name, H5P_DEFAULT)) < 0 || !bool_id)        return false;
    if ((bool_id = H5Oexists_by_name(loc_id, name, H5P_DEFAULT)) < 0 || !bool_id)        return false;
    return true;
}

bool is_cgef(string &filename) {
    hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    bool is_c = false;
    if(H5Lexists(file_id, "cellBin", H5P_DEFAULT))
        is_c = true;
    H5Fclose(file_id);
    return is_c;
}

bool is_bgef(string &filename) {
    hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    bool is_b = false;
    if(H5Lexists(file_id, "geneExp", H5P_DEFAULT))
        is_b = true;
    H5Fclose(file_id);
    return is_b;
}



bool h5AttrWrite(hid_t id, hid_t filetype, hid_t memtype, const char *name, int rank, hsize_t *dims, const void* value)
{
    hid_t attr_s = H5Screate_simple(rank, dims, nullptr);
    hid_t attr_d = H5Acreate(id, name, filetype, attr_s, H5P_DEFAULT, H5P_DEFAULT);
    herr_t ret = H5Awrite(attr_d, memtype, value);
    if(ret < 0)
    {
        printf("%s write err\n", name);
        return false;
    }
    H5Sclose(attr_s);
    H5Aclose(attr_d);
    return true;
}

hid_t h5DatasetWrite(hid_t gid, hid_t filetype, hid_t memtype, const char *name, int rank, hsize_t *dims, const void *data)
{
    hid_t space_id = H5Screate_simple(rank, dims, nullptr);
    hid_t data_id = H5Dcreate(gid, name, filetype, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    herr_t ret = H5Dwrite(data_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    if(ret < 0)
    {
        printf("%s write err\n", name);
        return false;
    }
    H5Sclose(space_id);
    return data_id;
}

size_t Rect_hash(const cv::Rect &rect)
{
    unsigned long tmp = 0, tx = rect.x, ty = rect.y;
    tmp = (tx << 40) | (ty << 16) | (rect.width << 8) | rect.height;
    return tmp;
}

bool Rectequal_to(const cv::Rect &a, const cv::Rect &b)
{
    return (a.x == b.x) && (a.y == b.y) && (a.width == b.width) && (a.height == b.height);
}

unsigned int parseResolutin(const string& filename) {
    std::unordered_map<string, unsigned int> pitch(
        {{"CL1", 900},  {"N1", 900},   {"V3", 715},  {"K2", 715},
         {"S2", 715},   {"S1", 900},   {"F3", 715},  {"F1", 800},
         {"V1", 800},   {"DP84", 715}, {"DP8", 850}, {"FP2", 500},
         {"SS2", 500},  {"FP1", 600},  {"E1", 700},  {"DP40", 700},
         {"G1", 700},   {"A", 500},    {"B", 500},   {"C", 500},
         {"D", 500},    {"U", 715},    {"V", 715},   {"W", 715},
         {"X", 715},    {"Y", 500},    {"P1", 715},  {"SS84", 715},
         {"FP21", 500}, {"SS1", 600}});

    auto pos = filename.find_last_of('/');
    if (pos == std::string::npos)
        pos = -1;

    unsigned int result = 0;
    std::string chip_prefix = filename.substr(pos+1, 4);
    while (chip_prefix.size() >= 1){
        if (pitch.count(chip_prefix) != 0)
        {
            result = pitch[chip_prefix];
            break;
        }
        chip_prefix.pop_back();
    }

    return result;
}

long tifread(cv::Mat &img, const string &strtif)
{
    TIFF* tif = TIFFOpen(strtif.c_str(), "r");
    if (tif) 
    {
        uint32_t w, h, ch;
        size_t npixels;
        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
        TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &ch);
        //int nf = TIFFNumberOfDirectories(tif);
        //int k = TIFFScanlineSize(tif);
        npixels = w * h;
        img.create(h,w,CV_8UC1);
        uchar *ptr = img.data;
        for (int row = 0; row < h; row++)
        {
            TIFFReadScanline(tif, ptr, row);
            ptr += w;
        }
        TIFFClose(tif);
        printf("img row:%d col:%d\n", img.rows, img.cols);
        return npixels;
    }
    return 0;
}