// header file
#include <fstream>  
#include <iostream> 
#include <stdio.h>
#include <bitset>
#include <vector>
#include <cstring>
#include <chrono>
#include <Rcpp.h>
#include "BedFileReader.h"

using namespace std;

/// @brief 
/// @param famName 
/// @param bimName 
/// @param bedName  

BedFileReader::BedFileReader(string famName, string bimName, string bedName){

    /* ---------- initialization ----------*/
    this->famName_temp = famName;
    this->bimName_temp = bimName;
    this->bedName_temp = bedName;
    this->m_line_counter = -1; // total number of individual 
    this->m_num_of_snps = -1; // total number of snps

    /* ---------- fam file reading ----------*/

    this->fam.open(famName);    // open file 

    // error messege
    if (!this->fam){
        error.open("reading_file.log");
        error << "CANT_OPEN_FAM_FILE4READ\n";
        error.close();

    }

    if (this->fam.is_open()){        // check if the file is open
        while(!this->fam.eof()){    // while end of the file
            string str;
            getline(this->fam, str);
            this->m_line_counter++; // read fam file to ge number of individuals
        }
        this->fam.close();  
        
    }
    
    cout << "Number of individuals: " << this->m_line_counter << endl;
    this->m_size_of_esi = (this->m_line_counter + 3)/4; //+3 for magic number & mode
    cout << "Number of bytes: " << this->m_size_of_esi << endl;

    /* ---------- bim file reading ----------*/

    this->bim.open(bimName);    

    // error messege
    if (!this->bim){
        error.open("reading_file.log", ios_base::app); // append error messege to log file
        error << "CANT_OPEN_BIM_FILE4READ\n";
        error.close();
    }

   if (this->bim.is_open()){
        while(!this->bim.eof()){
            string str;
            getline(this->bim, str);
            this->m_num_of_snps++;
        }
        this->bim.close();  
        
    }

    cout << "Number of snps: " << this->m_num_of_snps << endl;

    /* ---------- bed file reading ----------*/
        
    this->bed.open(bedName, ios::binary); // read it as binary

    //error messege
    if (!this->bed){
        error.open("reading_file.log", ios_base::app); 
        error << "CANT_OPEN_BED_FILE4READ\n";
        error.close();

    }

    char tmpbuf[3]; 
    memset(tmpbuf, '\0', sizeof(tmpbuf));
    this->bed.read(tmpbuf, sizeof(tmpbuf)); // three first bytes - permanent in bed file

}

// void BedFileReader::snp_index_func(){
    
//     this->bim.open(this->bimName_temp);

//     // error messege
//     if (!this->bim.is_open()){
//     error.open("reading_file.log", ios_base::app);
//     error << "BIM_FILE_READ_ERROR: cannot open the file\n";
//     error.close();
//     } 

//     size_t curLine = 0;
//     string line;

//     while(getline(this->bim, line)) { 
//         curLine++;
//         std::istringstream iss(line);

//         int column1;
//         string column2;

//         iss >> column1 >> column2; 
//         this->snp_index.insert({column2, curLine});
//     }

//     this->bim.close(); 
// }

void BedFileReader::snp_index_func(){
  this->bim.open(this->bimName_temp);
  // error messege
  if (!this->bim.is_open()){
  error.open("reading_file.log", ios_base::app);
  error << "BIM_FILE_READ_ERROR: cannot open the file\n";
  error.close();
  }
  size_t curLine = 0;
  string line;
  while(getline(this->bim, line)) {
    curLine++;
    std::istringstream iss(line);
    int column1;
    string column2;
    iss >> column1 >> column2;
    this->snp_index.insert({column2, curLine});
  }
  this->bim.close();
}

int BedFileReader::findSnpIndex(string snpName){
    
    this->bim.open(this->bimName_temp); //dictionary (map) : snpid(std::string) as key, snpidx as value

    // error messege
    if (!this->bim.is_open()){
    error.open("reading_file.log", ios_base::app);
    error << "BIM_FILE_READ_ERROR: cannot open the file\n";
    error.close();
    } 

    size_t curLine = 0;
    string line;
    
    while(getline(this->bim, line)) { 
        curLine++;
        
        if (line.find(snpName, 0) != string::npos) { //only 2nd col
            //cout << "found: " << this->snpName << " line: " << curLine << endl;
            this->snpIndex = curLine;
        }
    }
    this->bim.close();  

    //cout << this->snpIndex << endl;
    return this->snpIndex;
}

// int BedFileReader::findSnpIndex(string snpName){
//   snpIndex = this->snp_index.find(snpName)->second;
//   //cout << snpIndex << endl;
//   return snpIndex;
// }

vector<int> BedFileReader::readOneSnp(int snpIndex){

    this->snpIndex = snpIndex;
    /* ---------- initialization ----------*/
    int bits_val[MY_CHAR_BIT];
    vector<int> temp_snp_info0(m_line_counter, 0);
    vector<int> temp_snp_info1(m_line_counter, 0); 
    vector<char> encoded_snp_info(m_size_of_esi, 0); 
    vector<char> buff(m_size_of_esi, 0); 
    
    size_t individuals_counter = 0;
    size_t pos = (snpIndex-1) * this->m_size_of_esi + 3;
    size_t pos_cur = this->bed.tellg();

    // error messege
    if (!this->bed.is_open()){
    error.open("reading_file.log", ios_base::app);
    error << "1. BED_FILE_READ_ERROR: cannot open the file\n";
    error.close();
    } 

    if(!pos_cur){
        error.open("reading_file.log", ios_base::app); 
        error << "2. BED_FILE_READ_ERROR: cannot get the position\n";
        error.close();
    }

    if(!this->bed.seekg(pos, ios::beg)){
        error.open("reading_file.log", ios_base::app); 
        error << "3. BED_FILE_READ_ERROR: cannot get seek the position\n";
        error.close();
    }

    this->bed.read((char*) &buff[0],this->m_size_of_esi);

    for (size_t i = 0; i < this->m_size_of_esi; i++) {	
    //===============================================================					
    //=== This part converts Byte "buff[i]" to bits values "bits_val"
    //=== for example byte buff[0] = "w" ->  bits_val = 11101011

        memset(bits_val, 0, sizeof(bits_val));

        int k = MY_CHAR_BIT;  //8
            while (k > 0){
                -- k;
                bits_val[k] = (buff[i]&(1 << k) ? 1 : 0);
            }
    
        //=== here interpret Bit information "bits_val" to snps and count it - decode it
        decode_byte(bits_val, &individuals_counter, &temp_snp_info0[0], &temp_snp_info1[0]);
    
    }

  return temp_snp_info0;
}

void BedFileReader::readAllSnp(string fileName){

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    ofstream matrix;
    matrix.open(fileName, ios_base::app);
    for (int i= 1; i <= this->m_num_of_snps; i++){
        vector<int> a0 = readOneSnp(i); 
        for (int j= 0; j < this->m_line_counter; j++){
            matrix << a0[j] << " ";
        }
        matrix << "\n";
        
    }
    matrix.close();
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout << "Time take to read all snps is " << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << "[ms]" << endl;

}

// ADDED
void BedFileReader::readSomeSnp(string fileName, vector<string> snpNames){

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    ofstream matrix2;
    matrix2.open(fileName, ios_base::app);

    for (const string& snpName : snpNames) {
        cout << snpName << endl;
        int i = findSnpIndex(snpName);
        vector<int> a0 = readOneSnp(i);

        for (int j = 0; j < this->m_line_counter; j++) {
            matrix2 << a0[j] << " ";
        }
        matrix2 << "\n";
    }
    matrix2.close();
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout << "Time take to read some snps is " << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << "[ms]" << endl;
}

vector<float> BedFileReader::calculatePRS(vector<string> snpList, vector<float> betaList){

    vector<int> oneSnp; 
    vector<float> PRS(m_line_counter, 0);


    for (int i = 0; i < snpList.size(); i++){
        //snpIndex = findSnpIndex(snpList[i]);
        
        snpIndex = this->snp_index.find(snpList[i])->second;
        oneSnp = readOneSnp(snpIndex);
        for (int j = 0; j < oneSnp.size(); j++){
            PRS[j] += oneSnp[j] * betaList[i];
        }
    }
    return PRS;
}


void BedFileReader::decode_byte(int* bits_val, size_t * individuals_counter, int* temp_snp_info0, int* temp_snp_info1){
	//int flag = 0;
	for (int i = 0; i < 4; ++i)
	{
		if(bits_val[i*2] == 0 && bits_val[i*2+1] == 0)
		{
			*individuals_counter += 1;
			if (*individuals_counter > m_line_counter)
				return;
			temp_snp_info0[*individuals_counter - 1] = 2;
			temp_snp_info1[*individuals_counter - 1] = 0;
			//flag = 1 ;//Homozegote 1 for example GG      write 20 ; 00 - will point to [0] letter   + 2 to [0]

		}
		else if(bits_val[i*2] == 1 && bits_val[i*2+1] == 1)
		{
			*individuals_counter += 1;
			if (*individuals_counter > m_line_counter)
				return;
			temp_snp_info0[*individuals_counter - 1] = 0;
			temp_snp_info1[*individuals_counter - 1] = 2;
			//flag = 2 ;//Homozegote 2 for example AA      write 02 ; 11 - will point to [1] letter   + 2 to [1]
		}
		else if(bits_val[i*2] == 1 && bits_val[i*2+1] == 0)
		{
			*individuals_counter += 1;		
			if (*individuals_counter > m_line_counter)
				return;
			temp_snp_info0[*individuals_counter - 1] = 9;
			temp_snp_info1[*individuals_counter - 1] = 9;
			//flag = 3 ; //Missing value                   nothing to add - write 99 ;
		}
		else if(bits_val[i*2] == 0 && bits_val[i*2+1] == 1)
		{
			*individuals_counter += 1;
			if (*individuals_counter > m_line_counter)
				return;
			temp_snp_info0[*individuals_counter - 1] = 1;
			temp_snp_info1[*individuals_counter - 1] = 1;

			//flag = 4 ; //Heterozegote for example AG  or GA     write 11 ; 01 - will point to [0] and [1] letter   +1 +1 to [1]
		}
		// else
		// 	flag = 5 ; //Error
	}
}

void BedFileReader::close_bed(){
  if(this->bed.is_open()){
    this->bed.close();
  }
}
   
// [[Rcpp::export]]
/// create an external pointer to a Uniform object
RcppExport SEXP BedFileReader__new(SEXP famName_, SEXP bimName_, SEXP bedName_) {
  // convert inputs to appropriate C++ types
  string famName = Rcpp::as<string>(famName_);
  string bimName = Rcpp::as<string>(bimName_);
  string bedName = Rcpp::as<string>(bedName_);

  // create a pointer to an Uniform object and wrap it as an external pointer
  Rcpp::XPtr<BedFileReader> reader( new BedFileReader(famName, bimName, bedName), true );

  // return the external pointer to the R side
  return reader;
}

RcppExport SEXP BedFileReader__snp_index_func(SEXP xp){
    Rcpp::XPtr<BedFileReader> reader(xp);
    reader->snp_index_func();
}

RcppExport SEXP BedFileReader__readOneSnp(SEXP xp, SEXP snpIndex_){
    Rcpp::XPtr<BedFileReader> reader(xp);

    int snpIndex = Rcpp::as<int>(snpIndex_);

    vector<int> oneSnp = reader->readOneSnp(snpIndex);
    
    return Rcpp::wrap(oneSnp);
}

RcppExport SEXP BedFileReader__calculatePRS(SEXP xp, SEXP snpList_, SEXP betaList_){
    Rcpp::XPtr<BedFileReader> reader(xp);

    vector<string> snpList = Rcpp::as< vector<string> >(snpList_);
    vector<float> betaList = Rcpp::as< std::vector<float> >(betaList_);

    vector<float> PRS = reader->calculatePRS(snpList, betaList);

    return Rcpp::wrap(PRS);
}

RcppExport SEXP BedFileReader__readAllSnp(SEXP xp, SEXP filename_){
    Rcpp::XPtr<BedFileReader> reader(xp);

    string filename = Rcpp::as<string>(filename_);

    reader->readAllSnp(filename);
}

//ADDED
RcppExport SEXP BedFileReader__readSomeSnp(SEXP xp, SEXP filename_, SEXP snpNames_){
    Rcpp::XPtr<BedFileReader> reader(xp);

    string filename = Rcpp::as<string>(filename_);
    vector<string> snpNames = Rcpp::as<vector<string>>(snpNames_);

    reader->readSomeSnp(filename, snpNames);
}

RcppExport SEXP BedFileReader__findSnpIndex(SEXP xp, SEXP snpName_){
    Rcpp::XPtr<BedFileReader> reader(xp);

    string snpName = Rcpp::as<string>(snpName_);

    int snpIndex = reader->findSnpIndex(snpName);
    return Rcpp::wrap(snpIndex);
}


