#ifndef _BEDFILEREADER_H        
#define _BEDFILEREADER_H 

#include <fstream>  
#include <iostream> 


using namespace std;


#define MY_CHAR_BIT 8


class BedFileReader{

    public:
    /* ---------- Constructor ----------*/
    BedFileReader(string famName, string bimName, string bedName);

    /* ---------- Member functions ----------*/
    int findSnpIndex(string snpName);
    vector<int> readOneSnp(int snpIndex);
    void snp_index_func();
    void readAllSnp(string fileName);
    void readSomeSnp(string fileName, vector<string> snpNames);
    void close_bed();
    vector<float> calculatePRS(vector<string> snpList, vector<float> betaList);

    
    private:

    /* ---------- Data members ----------*/
    
    string famName_temp;
    string bimName_temp;
    string bedName_temp;
    
    // Streams for files
    ifstream fam;
    ifstream bim;
    ifstream bed;
    ofstream error;

    // Variables
    size_t m_line_counter; // Total number of individual 
    size_t m_size_of_esi; // Total number of bytes that hold all genotype info per one snp = one line length!!!
    size_t m_num_of_snps; // Total number of Snps
    map<string, int> snp_index;
    string snpName;
    vector<string> snpNames;
    int snpIndex;

    /* ---------- Member functions ----------*/
    void decode_byte(int* bits_val,size_t * individuals_counter, int* temp_snp_info0, int* temp_snp_info1);

    };

    #endif //_BEDFILEREADER_H        
