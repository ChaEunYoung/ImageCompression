#define _CRT_SECURE_NO_WARNINGS
#include<iostream>
#include<fstream>
#include<vector>
#include<iomanip>
#include<string>
#include<stdlib.h>
#include<bitset>
#include"RLC.h"
#include"DPCM.h"
using namespace std;

#define header_size 54
#define MV_NUMBER 16
#define MACRO_NUMBER 8
#define P_VALUE 15
#define W 352 //44
#define H 288 //36 (44*36 = 1584)
#define PI 3.141592
#define FILE_PATH "C:\\Users\\user\\Source\\Repos\\Video Compression\\Video Compression\\Video\\"
#define OUT_FILE_PATH "C:\\Users\\user\\Source\\Repos\\Video Compression\\Video Compression\\outFile\\"
//#define width 256
//#define height 256

int Qtable[8 * 8] = {
16, 11, 10, 16, 24, 40, 51, 61,
12, 12, 14, 19, 26, 58, 60, 55,
14, 13, 16, 24, 40, 57, 69, 56,
14, 17, 22, 29, 51, 87, 80, 63,
18, 22, 37, 56, 68, 109, 103, 77,
24, 35, 55, 64, 81, 104, 113, 92,
49, 64, 78, 87, 103, 121, 120, 101,
72, 92, 95, 98, 112, 100, 103, 99 };

struct node
{
	int characters;
	unsigned int frequency;
	string code;
	node * leftChild;
	node * rightChild;
};
int frame_num = 0;
vector<node> nodeArray;
vector<node> dpcmHuffmanCode;
vector<node> rlcHuffmanCode;
vector<node> symbolArry;
string totalVideoCode;
string copy;

void showHeaderInfo(char* input) {
	cout << "헤더 정보 :" << endl;
	cout << "filesize: " << *(int*)(input + 2) << endl;
	cout << "비트맵 시작 위치: " << *(int*)(input + 10) << endl;
	cout << "header_size: " << *(int*)(input + 14) << endl;
	cout << "Width: " << *(int*)(input + 18) << endl;
	cout << "height: " << *(int*)(input + 22) << endl;
	cout << "픽셀당 비트수: " << *(short*)(input + 28) << endl;
	cout << endl;
}
void showValue(unsigned char** input, int w, int h) {
	cout << endl;
	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			cout << setw(5) << (unsigned int)input[i][j];
		}
		cout << endl;
	}
}
void showValue(double** input, int w, int h) {
	cout << endl;
	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			cout << setw(10) << input[i][j];
		}
		cout << endl;
	}
}
void showValue(double* input, int w, int h) {
	cout << endl;
	for (int i = 0; i < h; i++) {
		for (int j = 0; j < w; j++) {
			cout << setw(10) << input[i*8+j];
		}
		cout << endl;
	}
}

string decToBinary(float input) {
	string result;
	//cout << input << "\t";
	int num = input;
	if (num == 0) { result.insert(0, "0"); return result; }
	if (input < 0) {
		num = num * -1;
	}
	while (num >= 1) {
		int binary = num % 2;
		if (binary == 0) result.insert(0, "0");
		else result.insert(0, "1");
		num = num / 2;
	}

	if (input < 0) {
		for (int i = 0; i < result.size(); i++) {
			if (result[i] == '0') result[i] = '1';
			else result[i] = '0';
		}

	}

	return result;
}

string* char_to_string(char* input, int data_size) {
	string result;
	for (int k = 0; k < data_size; k++) {
		string temp_string;
		char temp_char = input[k];
		for (int i = 0; i < 8; i++) {
			if ((temp_char & 0x80) == 0x80) {
				temp_string.append(1, '1');
			}
			else temp_string.append(1, '0');
			temp_char = temp_char << 1;
		}
		result += temp_string;
	}
	return &result;
}

void char_to_string(char* input, int data_size, string& out) {
	string result;
	for (int k = 0; k < data_size; k++) {
		string temp_string;
		char temp_char = input[k];
		for (int i = 0; i < 8; i++) {
			if ((temp_char & 0x80) == 0x80) {
				temp_string.append(1, '1');
			}
			else temp_string.append(1, '0');
			temp_char = temp_char << 1;
		}
		result += temp_string;
	}
	out = result;
}

string int_to_string(int* input, int data_size) {
	string result;
	for (int k = 0; k < data_size; k++) {
		string temp_string;
		int  temp_char = input[k];
		for (int i = 0; i < 8 * sizeof(int); i++) {
			if ((temp_char & 0x80000000) == 0x80000000) {
				temp_string.append(1, '1');
			}
			else temp_string.append(1, '0');
			temp_char = temp_char << 1;
		}
		result += temp_string;
	}
	return result;
}

string short_to_string(vector<short> input, int data_size) {
	string result;
	for (int k = 0; k < data_size; k++) {
		string temp_string;
		short  temp_char = input[k];
		for (int i = 0; i < 8 * sizeof(short); i++) {
			if ((temp_char & 0x8000) == 0x8000) {
				temp_string.append(1, '1');
			}
			else temp_string.append(1, '0');
			temp_char = temp_char << 1;
		}
		result += temp_string;
	}
	return result;
}

int String_to_int(string str8char)
{
	int tempc = 0;
	for (int i = 0; i < 8*sizeof(int); i++)
	{
		if (str8char[i] == '0')
			tempc += 0;
		else
			tempc += 1;
		if (i != 8 * sizeof(int)-1)
			tempc = tempc << 1;
	}
	return tempc;
}

char String_to_char(string str8char)
{
	char tempc = 0;
	for (int i = 0; i < 8; i++)
	{
		if (str8char[i] == '0')
			tempc += 0;
		else
			tempc += 1;
		if (i != 8 - 1)
			tempc = tempc << 1;
	}
	return tempc;
}

short String_to_short(string str8char)
{
	short tempc = 0;
	for (int i = 0; i < 8* sizeof(short); i++)
	{
		if (str8char[i] == '0')
			tempc += 0;
		else
			tempc += 1;
		if (i != 8* sizeof(short) - 1)
			tempc = tempc << 1;
	}
	return tempc;
}

void depthFirstSearch(node * tempRoot, string s, bool rlcDpcm)
{
	node * root1 = new node;
	root1 = tempRoot;

	root1->code = s;
	if (root1 == NULL)
	{

	}
	else if (root1->leftChild == NULL && root1->rightChild == NULL)
	{
		cout << "\t" << root1->characters << "\t" << root1->code << "\t" << "\t" << root1->frequency << endl;
		if (rlcDpcm) { // true-> dpcm
			node temp;
			temp.code = root1->code;
			temp.characters = root1->characters;
			dpcmHuffmanCode.push_back(temp);
		}
		else {          //rlc
			node temp;
			temp.code = root1->code;
			temp.characters = root1->characters;
			rlcHuffmanCode.push_back(temp); // rlc 자리
		}
	}
	else
	{
		root1->leftChild->code = s.append("0");
		s.erase(s.end() - 1);
		root1->rightChild->code = s.append("1");
		s.erase(s.end() - 1);

		depthFirstSearch(root1->leftChild, s.append("0"), rlcDpcm);
		s.erase(s.end() - 1);
		depthFirstSearch(root1->rightChild, s.append("1"), rlcDpcm);
		s.erase(s.end() - 1);
	}

}

node extractMin()
{	
	unsigned int min = UINT_MAX;
	vector<node>::iterator iter, position;
	for (iter = nodeArray.begin(); iter != nodeArray.end(); iter++)
	{
		if (min > (*iter).frequency)
		{
			position = iter;
			min = (*iter).frequency;
		}
	}

	node tempNode = (*position);
	nodeArray.erase(position);

	return tempNode;
}

node getHuffmanTree()
{
	while (!nodeArray.empty())
	{
		node * tempNode = new node;
		node * tempNode1 = new node;
		node * tempNode2 = new node;
		*tempNode1 = extractMin();
		*tempNode2 = extractMin();

		tempNode->leftChild = tempNode1;
		tempNode->rightChild = tempNode2;
		tempNode->frequency = tempNode1->frequency + tempNode2->frequency;
		nodeArray.push_back(*tempNode);

		//Root Node만 남았으므로 Huffman Tree 완성
		if (nodeArray.size() == 1)  break;
	}
	
	
	return nodeArray[0];
}

void getHuffmanCode(dpcm* input, int input_size) //dpcm
{
	//vector<node> temp_dpcmHuffmanCode;
	int size = 11; //데이터에 사용되는 문자의 갯수 입력
	unsigned int tempInt;
	int *fre = new int[11];
	for (int i = 0; i < 11; i++) fre[i] = 0;

	for (int i = 0; i < input_size; i++) {
		for (int j = 0; j < 11; j++) {
			if (input[i].size == j) {
				fre[j]++;
			}
		}
	}
	
	
	//각 문자별 빈도 수 노드 생성
	for (int i = 0; i < size; i++)
	{
		node tempNode;
		tempNode.characters = i;
		tempNode.frequency = fre[i];
		tempNode.leftChild = NULL;
		tempNode.rightChild = NULL;
		nodeArray.push_back(tempNode);
	}

	
	//Huffman Tree 생성
	node root = getHuffmanTree();
	
	
	cout << endl;
	cout << endl;
	cout << "\t" << "size" << "\t" << "Code" << "\t" << "\t" << "frequancy" << endl;
	cout << "\t" << "------------------------------------" << endl;
	
	
	//Huffman Conding Table 생성
	depthFirstSearch(&root, "", true);
	
}

void getHuffmanCode(vector<node> input)  //rlc
{
	
	int size = input.size(); //데이터에 사용되는 문자의 갯수 입력
	unsigned int tempInt;

	//각 문자별 빈도 수 노드 생성
	for (int i = 0; i < size; i++)
	{
		node tempNode;
		tempNode.characters = input[i].characters;
		tempNode.frequency = input[i].frequency;
		tempNode.leftChild = NULL;
		tempNode.rightChild = NULL;
		nodeArray.push_back(tempNode);
	}

	//Huffman Tree 생성
	node root = getHuffmanTree();
	
	cout << endl;
	cout << endl;
	cout << "\t" << "symbol1" << "\t" << "Code" << "\t" << "\t" << "frequancy" << endl;
	cout << "\t" << "------------------------------------" << endl;
	
	//Huffman Conding Table 생성
	depthFirstSearch(&root, "",false);
}

string getDcCode(dpcm* dpcmInput, int dcArray_size) { //dcCode 생성
	string dcCode;
	for (int i = 0; i < dcArray_size; i++) {

		for (int j = 0; j < dpcmHuffmanCode.size(); j++) {
			if (dpcmInput[i].size == dpcmHuffmanCode[j].characters) {
				dcCode += dpcmHuffmanCode[j].code;
				break;
			}
		}

		dcCode += decToBinary(dpcmInput[i].amplitude);
	}
	
	//for(int i =0; i< dpcmHuffmanCode.size();i++)cout << dpcmHuffmanCode[i].characters << endl;
	//getchar();
	
	while (!nodeArray.empty()) nodeArray.pop_back();
	return dcCode;
}

string DPCM(vector<double> input) { // 입력: dcArray
	string dcCode;
	dpcm* difference = new dpcm[input.size()];
	difference[0] = dpcm(input[0]);
	for (int i = 0; i < input.size() - 1; i++) {
		difference[i + 1] = dpcm(input.at(i + 1) - input.at(i));
	}

	getHuffmanCode(difference, input.size()); // 코드 백터 생성
	return getDcCode(difference, input.size()); //DPCM 코드화 -> huffman+amplitude
}

string getAcCode(rlc* rlcInput,vector<double> acArray) {
	string acCode;
	for (int i = 0; i < acArray.size() / 2; i++) {

		for (int j = 0; j < rlcHuffmanCode.size(); j++) {
			if (rlcInput[i].getsymbol() == rlcHuffmanCode[j].characters) {
				acCode += rlcHuffmanCode[j].code; // sym1에 대한 코드
				break;
			}
		}
		acCode += decToBinary(rlcInput[i].amplitude); // 그때의 amp
	}
	//while (!rlcHuffmanCode.empty()) rlcHuffmanCode.pop_back();
	while (!nodeArray.empty()) nodeArray.pop_back();
	while (!symbolArry.empty()) symbolArry.pop_back();
	return acCode;
}

void getsymbolArray(rlc* input,int acarray_size) {
	for (int i = 0; i < acarray_size / 2; i++) {

		for (int k = 0; k < 1515; k++) {
			if (input[i].getsymbol() == k) {

				bool newsymbol = true;
				for (int j = 0; j < symbolArry.size(); j++) {
					if (input[i].getsymbol() == symbolArry[j].characters) {
						symbolArry[j].frequency++;
						newsymbol = false;
						break; break;
					}
				}
				if (newsymbol) {
					node temp;
					temp.characters = k;
					temp.frequency = 1;
					symbolArry.push_back(temp);
					break;
				}
			}
		}
	}
}

string RLC_symbol(vector<double> input) { //acarray-> accode
	string temp_acCode;
	vector<double> temp_input = input;
	rlc* difference = new rlc[input.size() / 2];
	for (int i = 0, ii = 0; i < input.size(); ii++, i += 2) {
		difference[ii] = rlc(input[i], input[i + 1]); // len, amp
	}
	getsymbolArray(difference, input.size());
	getHuffmanCode(symbolArry);
	temp_acCode= getAcCode(difference, temp_input);
	return temp_acCode;
}

vector<double> RLC(double* input) {
	vector<double> rlc;
	for (int i = 0; i < 63; i++) {
		int count = 0;
		if (input[i] == 0) {
			while (input[i] == 0) {
				count++;
				if (count == 15) {
					break;
				}
				if (i == 62) {
					rlc.push_back(0);
					rlc.push_back(0);

					return rlc;
				}
				i++;
			}
		}
		rlc.push_back(count);
		rlc.push_back(input[i]);
	}
	return rlc;
}

void zigzagScan(double** q_input, double* output, bool a) {
	//Q_count++;
	double* zigzag = new double[63];
	//double* zigzag = new double[63];
	vector<double> temp;
	int index = 0;
	int i = 0, j = 1;
	while (i != 7 || j != 7) {
		if (j > i) {
			int temp_i = i;
			int temp_j = j;
			while (temp_i != j) {
				zigzag[index] = q_input[i][j]; index++;
				i++; j--;
			}
			zigzag[index] = q_input[i][j]; index++;
			if (i < 7) { i++; }
			else { j++; }
		}
		else if (j < i) {
			int temp_i = i;
			int temp_j = j;
			while (temp_i != j) {
				zigzag[index] = q_input[i][j]; index++;
				i--; j++;
			}
			zigzag[index] = q_input[i][j]; index++;
			if (j < 7) { j++; }
			else { i++; }
		}
	}
	zigzag[index] = q_input[i][j];
	for (int i = 0; i < 63; i++) {
		output[i] = zigzag[i];
	}
	/*
	if (a) {
		temp = RLC(zigzag);
		for (int i = 0; i < temp.size(); i++) {
			acarray.push_back(temp.at(i));
			//cout << temp[i] << "\t";
		}

	}
	else {
		iz = zigzag;
	}
	*/
	delete[] zigzag;
}

void quantization(double** input, double **output, double outDC) {

	for (int i = 0; i < MACRO_NUMBER; i++) {
		double temp = 0;
		for (int j = 0; j < MACRO_NUMBER; j++) {

			temp = ((double)input[i][j] / (double)(Qtable[i * MACRO_NUMBER + j]));
			output[i][j] = round(temp);
			if (i == 0 && j == 0) {
				outDC = round(temp);
				
			}
		}
	}
}

void Div_ImgBlock_MV(unsigned char ** data, unsigned char ** block, int m, int n)
{
	for (int i = 0; i < MV_NUMBER; i++)
	{
		for (int j = 0; j < MV_NUMBER; j++)
		{
			block[i][j] = data[m + i][n + j]; // mn은 실제 이미지의 좌표
		}
	}
}

void Div_ImgBlock_MV(double ** data, double** block, int m, int n)
{
	for (int i = 0; i < MACRO_NUMBER; i++)
	{
		for (int j = 0; j < MACRO_NUMBER; j++)
		{
			block[i][j] = data[m + i][n + j]; // mn은 실제 이미지의 좌표
		}
	}
}

void IDCT(double** dct, unsigned char** output, int w, int h) {
	int width = w; int height = h;
	int N = 8, M = 8;;
	int image_size = width * height;
	int mcrNb_x = width / N; // 32
	int mcrNb_y = height / M; // 32
	double sum = 0, temp = 0;
	int  u, n, m, v;

	for (int mcr_y = 0; mcr_y < mcrNb_y; mcr_y++) {
		for (int mcr_x = 0; mcr_x < mcrNb_x; mcr_x++) {
			for (int i = 0; i < N; i++) {
				n = mcr_y * N + i;
				for (int j = 0; j < M; j++) {
					m = mcr_x * M + j;
					sum = 0;
					for (int k = 0; k < N; k++) {
						v = mcr_y * N + k;
						for (int l = 0; l < M; l++) {
							u = mcr_x * M + l;
							double theta_x = (double)(2.*j + 1)*l*PI / (2.*M);
							double theta_y = (double)(2.*i + 1)*k*PI / (2.*N);
							double ck_y = (k) ? 0.5 : sqrt((double)1.0 / (double)N);
							double ck_x = (l) ? 0.5 : sqrt((double)1.0 / (double)N);
							sum += ck_x * ck_y*(double)cos(theta_x)*cos(theta_y)*dct[v][u];
						}

					}
					if (sum > 255)  output[n][m] = 255;
					else if (sum < 0)  output[n][m] = 0;
					else output[n][m] = sum;
				}
			}
		}
	}
}

void DCT(unsigned char** input, double** dct, int w, int h) {
	int width = w; int height = h;
	int N = 8, M = 8;;
	int image_size = width * height;
	int mcrNb_x = width / N; // 32
	int mcrNb_y = height / M; // 32
	double sum = 0;
	int  u, n, m, v;

	for (int mcr_y = 0; mcr_y < mcrNb_y; mcr_y++) {
		for (int mcr_x = 0; mcr_x < mcrNb_x; mcr_x++) {
			for (int k = 0; k < N; k++) {
				v = mcr_y * N + k;
				for (int l = 0; l < M; l++) {
					u = mcr_x * M + l;
					sum = 0;
					for (int i = 0; i < N; i++) {
						n = mcr_y * N + i;
						for (int j = 0; j < M; j++) {
							m = mcr_x * M + j;
							double theta_x = (double)(2.*j + 1)*l*PI / (2.*M);
							double theta_y = (double)(2.*i + 1)*k*PI / (2.*N);
							sum += (double)cos(theta_x)*cos(theta_y)*input[n][m];
						}
					}
					double ck_y = (k) ? 0.5 : sqrt((double)1.0 / (double)N);
					double ck_x = (l) ? 0.5 : sqrt((double)1.0 / (double)N);
					dct[v][u] = ck_x * ck_y * sum;
				}
			}
		}
	}
}

void get_Full_data(unsigned char* image, char* image_header) {
	string inputName = "Video_Compression";
	ifstream inputFile;
	char* header = new char[header_size];
	unsigned char* Full_R = new unsigned char[W*H * 25];
	int image_Count = 0;

	for (int i = 1; i <= 25; i++) {
		char* temp_header = new char[header_size];
		unsigned char* temp_RGB = new unsigned char[W*H * 3];
		string buf = to_string(i);
		if (i < 10) buf = "0" + buf;

		inputFile.open(FILE_PATH + inputName + buf + ".bmp", ios::binary);

		if (inputFile.is_open() == NULL) {
			cout << "File open error!" << endl;
			exit(0);
		}
		image_Count = W * H*(i - 1);

		inputFile.read((char*)temp_header, header_size);
		inputFile.read((char*)temp_RGB, W*H * 3);
		for (int n = 0, nn = 0; n < H*W * 3; nn++, n += 3) {
			Full_R[nn + image_Count] = temp_RGB[n];
		}

		delete[] temp_RGB;
		if (i != 25) delete[] temp_header;
		else {
			showHeaderInfo(temp_header);
			for (int n = 0; n < header_size; n++) {
				image_header[n] = temp_header[n];
			}
			delete[] temp_header;
		}
		inputFile.clear();
		inputFile.close();
	}
	for (int i = 0; i < W*H * 25; i++) {
		image[i] = Full_R[i];
	}


	delete[] header;
	delete[] Full_R;

	cout << "complete read data!" << endl;
}

string JPEG_encoding(unsigned char** input) {
	vector<double> dcArray;
	vector<double> acArray;
	vector<double> temp_ac;
	string dcCode;
	string acCode;
	string total_Code;
	string acHuffCode;
	string dcHuffCode;
	string temp_header;
	string dcHuff;
	string acHuff;
	short x = 0, y = 0;
	int count = 0;
	
	
	unsigned char** inputImage = new unsigned char*[H];
	int* sizeInfo = new int[6]; // 전체, dc크기, dchuff크기, achuff크기, dc심볼수, ac 심볼수
	double** dct = new double*[H];
	double** temp_block = new double*[MACRO_NUMBER];
	double** temp_qout = new double*[MACRO_NUMBER];
	double* temp_zout = new double[63];


	for (int i = 0; i < H; i++) {
		inputImage[i] = new unsigned char[W];
		dct[i] = new double[W];
	}
	for (int i = 0; i < MACRO_NUMBER; i++) {
		temp_block[i] = new double[MACRO_NUMBER];
		temp_qout[i] = new double[MACRO_NUMBER];
	}

	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			inputImage[i][j] = input[i][j];
		}
	}

	DCT(inputImage, dct, W, H);

	double temp_DC = 1;
	for (int i = 0; i < H / MACRO_NUMBER; i++) // MV_NUMBER 마크로블록의 크기 결국 256/8 = 32
	{
		for (int j = 0; j < W / MACRO_NUMBER; j++) // MV_NUMBER 마크로블록의 크기 결국 256 / 8 = 32
		{
			x = j * MACRO_NUMBER; // 실제 이미지의 좌표 마크로 블록 첫번째 인덱스
			y = i * MACRO_NUMBER; // 실제 이미지의 좌표 마크로 블록 첫번째 인덱스
			Div_ImgBlock_MV(dct, temp_block, y, x);
			quantization(temp_block, temp_qout, temp_DC);
			dcArray.push_back(temp_qout[0][0]);
	
			zigzagScan(temp_qout, temp_zout, true);
			temp_ac = RLC(temp_zout);
			for (int k = 0; k < temp_ac.size(); k++) {
				acArray.push_back(temp_ac[k]);
			}
		}
	}
	dcCode = DPCM(dcArray);
	acCode = RLC_symbol(acArray);
	
	total_Code = dcCode + acCode;
	
	const int dcSymNum = dpcmHuffmanCode.size();
	const int acSymNum = rlcHuffmanCode.size();
	char* dc_SymbolArray = new char[dcSymNum];
	char* ac_SymbolArray = new char[acSymNum];
	char* dc_CodeLen = new char[dcSymNum];
	char* ac_CodeLen = new char[acSymNum];


	/*
	허프만 코드 저장법:
	1. 먼저 허프만 코드를 배열에 담는다.
	2. 8비트를 모두 string으로 바꾼다 (char to string)
	3. 허프만 코드 사이즈를 맨 처음 알려준 후 읽을 때는 그 크기 만큼 읽는다.
	@하나의 프레임이 갖는 데이터(헤더는 첫 번째 프레임만)
	허프만 코드를 포함한 전체 데이터 크기(int)
	dc허프만 코드의 크기 (int), ac허프만 코드 크기
	dc심볼갯수 ac 심볼갯수
	dc 심볼 어레이
	dc 코드 len 어레이
	dc 코드
	 ac 도 동일
	데이터의 크기 
	*/


	for (int i = 0; i < dcSymNum; i++) {
		dcHuff += dpcmHuffmanCode[i].code;
		dc_SymbolArray[i] = (char)dpcmHuffmanCode[i].characters;
		dc_CodeLen[i] = dpcmHuffmanCode[i].code.size();
		
	}
	
	for (int ii = 0; ii < acSymNum; ii++) {
		acHuff += rlcHuffmanCode[ii].code;
		ac_SymbolArray[ii] = rlcHuffmanCode[ii].characters;
		ac_CodeLen[ii] = rlcHuffmanCode[ii].code.size();
	}
	
	string temp1;
	string temp2;
	string temp3;
	string temp4;
	
	char_to_string(dc_SymbolArray, dcSymNum, temp1);
	char_to_string(dc_CodeLen, dpcmHuffmanCode.size(), temp2);
	dcHuffCode = (temp1)+ (temp2)+dcHuff;
	

	char_to_string(ac_SymbolArray, rlcHuffmanCode.size(), temp3);
	char_to_string(ac_CodeLen, rlcHuffmanCode.size(), temp4);
	acHuffCode = (temp3) + (temp4) + acHuff;
	
	


	
	 //전체 크기
	sizeInfo[0] = sizeInfo[1] + sizeInfo[2] + sizeInfo[3] + sizeInfo[4] + sizeof(int) * 8 * 6;
	sizeInfo[1] = dcHuffCode.size(); // dc허프만 메모리
	sizeInfo[2] = acHuffCode.size(); // ac허프만 메모리
	sizeInfo[3] = dpcmHuffmanCode.size(); // dc 심볼수
	sizeInfo[4] = rlcHuffmanCode.size(); // ac 심볼수
	sizeInfo[5] = dcCode.size(); // ac 심볼수

	temp_header = int_to_string(sizeInfo, 6);
	
	
	total_Code = temp_header + dcHuffCode + acHuffCode + dcCode + acCode;
	
	
	

	for (int i = 0; i < H; i++) {
		delete[] dct[i];
		delete[] inputImage[i];
	}
	
	
	for (int i = 0; i < MACRO_NUMBER; i++) {
		delete[] temp_block[i];
		delete[] temp_qout[i];
	}
	delete[] dct;
	delete[] inputImage;
	delete[] temp_block;
	delete[] temp_qout;
	
	delete[] dc_SymbolArray;
	delete[] ac_SymbolArray;
	delete[] dc_CodeLen;
	delete[] ac_CodeLen;
	delete[] sizeInfo;
	delete[] temp_zout;
	
	
	
	while (!rlcHuffmanCode.empty()) rlcHuffmanCode.pop_back();
	while (!dpcmHuffmanCode.empty()) dpcmHuffmanCode.pop_back();
	cout << "complete JPEG Encoding" << endl;

	return total_Code;
}

//--------------------------------------JPEG_DECODE-------------------------------------------------------------
vector<node> dpcmHuffmanCode_docode;
vector<node> rlcHuffmanCode_decode;
string copy2;
unsigned char * outCopy = new unsigned char[W*H];

void getDcHuffmanCode(string code, int symNum) {
	int iter = 0;
	const int symN = symNum;
	char* CodeLen = new char[symN];
	for (int i = 0; i < symN; i++) {
		node temp;
		iter = i*8;
		
		temp.characters = (int)String_to_char(code.substr(iter, 8)); //charactor
		dpcmHuffmanCode_docode.push_back(temp);
		//cout << dpcmHuffmanCode_docode[i].characters << endl;

	}
	for (int i = 0; i < symN; i++) {
		iter += 8;
		CodeLen[i] = String_to_char(code.substr(iter, 8));
	}
	
	iter += 8;
	for (int i = 0; i < symN; i++) {
		dpcmHuffmanCode_docode[i].code = code.substr(iter, CodeLen[i]);
		iter += CodeLen[i];
	}

	
	/*
	cout << endl << "*디코드한 dpcm size " << endl;
	for (int i = 0; i < dpcmHuffmanCode_docode.size(); i++) {
		cout << "size: " << dpcmHuffmanCode_docode.at(i).characters << " ,code: " << dpcmHuffmanCode_docode.at(i).code << endl;
	}
	getchar();
	*/
	delete[] CodeLen;
}

void getAcHuffmanCode(string code, int symNum) {
	int iter = 0;
	const int symN = symNum;
	char* CodeLen = new char[symN];
	
	for (int i = 0; i < symN; i++) {
		node temp;
		iter = i * 8;

		temp.characters = (unsigned char)String_to_char(code.substr(iter, 8)); //charactor
		rlcHuffmanCode_decode.push_back(temp);
		//cout << rlcHuffmanCode_decode[i].characters << endl;
	}
	
	for (int i = 0; i < symN; i++) {
		iter += 8;
		CodeLen[i] = (int)String_to_char(code.substr(iter, 8));
	}
	
	iter += 8;
	for (int i = 0; i < symN; i++) {
		rlcHuffmanCode_decode[i].code = code.substr(iter, CodeLen[i]);
		iter += CodeLen[i];
	}
	
	/*
	cout << endl << "*디코드한 dpcm size " << endl;
	for (int i = 0; i < rlcHuffmanCode_decode.size(); i++) {
	cout << "size: " << rlcHuffmanCode_decode.at(i).characters << " ,code: " << rlcHuffmanCode_decode.at(i).code << endl;
	}
	getchar();
	*/
	
	delete[] CodeLen;
}

double BinaryToDec(string input)
{
	bool minus = false;

	if (input[0] == '0') { //음수
		minus = true;
		for (int i = 0; i < input.size(); i++) {
			if (input[i] == '0') input[i] = '1';
			else input[i] = '0';
		}
	} // 1의 보수 취해주고
	  // 꼭 -1 곱해주기!

	char result = 0;
	int size = input.size();
	for (int i = 0; i < size; i++) {

		if (input[i] == '1') {
			result += 1;
		}
		else {
			//result <<= 1;
		}
		if (i != size - 1) {
			result <<= 1;
		}
	}


	double resultresult = result;

	if (minus) { //음수
		resultresult = -resultresult;
	}  
	return resultresult;
}

void dc_decode(vector<double> &input, vector<double>& DC) {
	//cout << input.size() << endl;
	double temp = input[0] + input[1]; 
	DC.push_back(input[0]);
	DC.push_back(temp);

	for (int i = 1; i < input.size() - 1; i++) {
		temp = temp + input[i + 1];
		DC.push_back(temp);
	}
}

void getDcArr(string& input, vector<double>& DCArr) {
	
	string temp;
	string temp_dpcm;

	for (int i = 0; i < input.size(); i++) { // i++를 의심해봐

		temp.push_back(input.at(i)); //i 만큼 읽고
		for (int j = 0; j < dpcmHuffmanCode_docode.size(); j++) {

			if (temp == dpcmHuffmanCode_docode[j].code) {
				i++;
				// i가 알려준 사이즈 만큼 읽고
				if (dpcmHuffmanCode_docode[j].characters == 0) {
					temp_dpcm = input.substr(i, 1);

					//float temp_char = BinaryToDec(temp_dpcm);
					float temp_char = 0;
					DCArr.push_back(temp_char);

				}
				else {
					temp_dpcm = input.substr(i, dpcmHuffmanCode_docode[j].characters);
					char temp_char = BinaryToDec(temp_dpcm);
					DCArr.push_back((double)temp_char);
					i = i + dpcmHuffmanCode_docode[j].characters - 1;
				}
				while (!temp.empty()) temp.pop_back();
				while (!temp_dpcm.empty())temp_dpcm.pop_back();
				break;
			}
		}
		//if (input.size() - i < 8) { cout <<"stop at: "<< i << endl; break; } //남은 비트가 비트가 안되면 종료
	}

	//cout << "디코드한 dpcm: " << dcArray_decode.size() << endl;
	//for (int i = 0; i < dcArray_decode.size(); i++)cout << dcArray_decode[i] << "\t";
	//cout << endl<< endl;
	//for (int i = 0; i < dcArray_decode.size(); i++)cout << difference1[i].amplitude << "\t";
	
}

void getImage(vector<double> input, vector<double>& dc, unsigned char** outImage) {
	//cout << dc.size();
	//getchar();
	double** decodedata = new double*[H];
	unsigned char ** reImage = new unsigned char*[H];
	double* iz = new double[63];
	int DC_count = 0;

	for (int i = 0; i < H; i++) {
		decodedata[i] = new double[W];
		reImage[i] = new unsigned char[W];
	}

	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			decodedata[i][j] = 0;
		}
	}

	double** Qinput = new double*[8];
	for (int i = 0; i < 8; i++) Qinput[i] = new double[8];
	int count = 0;
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 8; j++) {
			Qinput[i][j] = count;
			count++;
		}
	}

	zigzagScan(Qinput,iz, false);
	
	vector<double> temp;
	vector<double> showdata;

	int macro_num = 0;
	int macro_num2 = 0;
	for (int macro_h = 0; macro_h < H; macro_h += 8) {
		for (int macro_w = 0; macro_w < W; macro_w += 8) {

			if (DC_count == (352/8 * 288/8) -1) {
				break; break;
			}

			double* result = new double[64];
			double* q_result = new double[64];

			for (int i = 0; i < 64; i++) {
				result[i] = 0;
				q_result[i] = 0;
			}

			for (int i = macro_num; i < macro_num + 63; i++) {
				temp.push_back(input[i]);
			}
			//for (int k = macro_num2; k < macro_num2 + 64; k++) {
			//	showdata.push_back(quantout[k]);
			//}

			// i_zigzag
			for (int i = 0; i < 63; i++) {
				int index = iz[i];
				result[index] = temp[i];
			}



			result[0] = dc[DC_count];
			


			for (int i = 0; i < 64; i++) q_result[i] = result[i] * (double)(Qtable[i]); // qualtization;

			for (int i = 0; i < 8; i++) {
				for (int j = 0; j < 8; j++) {

					decodedata[macro_h + i][macro_w + j] = q_result[i * 8 + j];
				}
			}

			delete[] result;
			delete[] q_result;

			while (!temp.empty()) temp.pop_back();
			//while (!showdata.empty()) showdata.pop_back();

			DC_count++;
			macro_num += 63;
			macro_num2 += 64;

		}
	}
	cout << "i_zigzag complete" << endl;
	cout << "i_quantization complete" << endl;
	IDCT(decodedata, reImage, W, H);
	cout << "IDCT complete" << endl;
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			outImage[i][j] = reImage[i][j];
		}
	}
	for (int i = 0; i < H; i++) {
		delete[] decodedata[i];
		delete[] reImage[i];
	}

	delete[] decodedata;
	delete[] reImage;
	delete[] iz;
}

vector<double> ac_decode(vector<double> input) { //AC
								
	vector<double> temp_arr;
	vector<double> total_arr;
	
	//cout << input.size() << endl;
	int index = 0;
	for (int i = 0; i < input.size(); i += 2) {

		if (input[i] == 0 && input[i + 1] == 0) {//00;
			for (int j = index; j < 63; j++) {
				temp_arr.push_back(0);
			}

			for (int v = 0; v < 63; v++) {
				total_arr.push_back(temp_arr[v]);
			}

			while (!temp_arr.empty()) temp_arr.pop_back();
			index = 0;
		}
		else if (input[i] == 0 && input[i + 1] != 0) { //amplitude only
			temp_arr.push_back(input[i + 1]);
			index++;
		}
		else if (input[i] != 0) { //zero amplitude
			if (input[i] == 15) { //zero 15;
				for (int zero = 0; zero < input[i]; zero++) {
					temp_arr.push_back(0);
					index++;
				}
			}
			else {                  //zero amplitude
				for (int zero = 0; zero < input[i]; zero++) {
					temp_arr.push_back(0);
					index++;
				}
				temp_arr.push_back(input[i + 1]);
				index++;
			}
		}
		if (index == 63) {
			index = 0;
			for (int v = 0; v < 63; v++) {
				total_arr.push_back(temp_arr[v]);
			}
			while (!temp_arr.empty()) temp_arr.pop_back();
		}
	}
	
	return total_arr;
}

void getAcArr(string& input, vector<double>& output) {
	string temp;
	string temp_rlc;

	for (int i = 0; i < input.size(); i++) {
		temp.push_back(input.at(i)); //i 만큼 읽고

		for (int j = 0; j < rlcHuffmanCode_decode.size(); j++) {
			if (temp == rlcHuffmanCode_decode[j].code) {
				i++;
				double zero_len;
				double ampl_size;

				if (rlcHuffmanCode_decode[j].characters == 0) {
					zero_len = 0;
					ampl_size = 1;
					//temp_rlc = input.substr(i, ampl_size);

					//float temp_char = BinaryToDec(temp_dpcm);
					//float temp_char = 0;
					output.push_back(0);
					output.push_back(0);
				}
				else {
					if (rlcHuffmanCode_decode[j].characters > 9) {
						zero_len = rlcHuffmanCode_decode[j].characters / 10;
						ampl_size = rlcHuffmanCode_decode[j].characters % 10;
					}
					else {
						zero_len = 0;
						ampl_size = rlcHuffmanCode_decode[j].characters;
					}

					temp_rlc = input.substr(i, ampl_size);
					char temp_char = BinaryToDec(temp_rlc);
					output.push_back(zero_len);
					if (ampl_size == 0) {
						ampl_size = 1;
						output.push_back(0);
					}
					else { output.push_back((double)temp_char); }

					i = i + ampl_size - 1;
				}

				while (!temp.empty()) temp.pop_back();
				while (!temp_rlc.empty())temp_rlc.pop_back();
				break;
			}
		}
		if (input.size() - i < 8) break; //남은 비트가 비트가 안되면 종료
	}


	//outdata = i_zigzagScan(ac_decode(acArray_decode));
}

void JPEG_decoding(string& input, unsigned char** out) {
	cout << "start JPEG decoding" << endl;

	int* temp_header = new int[5];
	int temp_index = 0;
	vector<double> DC;
	vector<double> dcArray_decode; // 36*44 = 1584;
	vector<double> acArray_decode;

	string tempH;
	string DPCM_huff_data;
	string RLC_huff_data;
	string dcCode_decode;
	string acCode_decode;

	unsigned char** outdata = new unsigned char*[H];
	for (int i = 0; i < H; i++) {
		outdata[i] = new unsigned char[W];
	}

	// 허프만 정보 받아오기
	for (int i = 0; i < 6; i++) {
		tempH = (input).substr(i * sizeof(int) * 8, sizeof(int) * 8);
		temp_header[i] = String_to_int(tempH);
		temp_index += sizeof(int) * 8;
		//cout << "info: " << temp_header[i] << endl;
	}


	const int temp_data_size = temp_header[0];
	const int temp_DcHuffM = temp_header[1];
	const int temp_AcHuffM = temp_header[2];
	const int temp_DcSymNum = temp_header[3];
	const int temp_AcSymNum = temp_header[4];
	const int temp_DCsize = temp_header[5];
	int temp_ACsize = temp_data_size - temp_DcHuffM - temp_AcHuffM - (sizeof(int) * 8 * 6) - temp_DCsize;

	// dc 허프만 추출
	DPCM_huff_data = (input).substr(temp_index, temp_DcHuffM);
	RLC_huff_data = (input).substr(temp_index + temp_DcHuffM, temp_AcHuffM);


	getDcHuffmanCode(DPCM_huff_data, temp_DcSymNum);
	getAcHuffmanCode(RLC_huff_data, temp_AcSymNum);

	dcCode_decode = (input).substr(temp_index + temp_DcHuffM + temp_AcHuffM, temp_DCsize);
	acCode_decode = (input).substr(temp_index + temp_DcHuffM + temp_AcHuffM + temp_DCsize, temp_ACsize);


	// Dc
	getDcArr(dcCode_decode, dcArray_decode);

	dc_decode(dcArray_decode, DC);


	//AC
	getAcArr(acCode_decode, acArray_decode); // acArray_decode 생성

	getImage(ac_decode(acArray_decode), DC, outdata);

	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			out[i][j] = outdata[i][j];
		}
	}
	
	for (int i = 0; i < H; i++) {
		delete[] outdata[i];
	}
	delete[] outdata;
	
	//delete[] temp_header;
	//cout << "df" << endl;
	while (!dpcmHuffmanCode_docode.empty()) dpcmHuffmanCode_docode.pop_back();
	while (!rlcHuffmanCode_decode.empty()) rlcHuffmanCode_decode.pop_back();

}

//---------------------------------------VIDEO_ENCODING---------------------------------------------------------

void Get_MotionVector(vector<short>& moVec, unsigned char ** r_frame, unsigned char** t_frame, int width, int height)
{
	moVec.clear(); // 각 프레임마다 clear
	short mv_x = 0, mv_y = 0, x = 0, y = 0;
	unsigned char ** cp_block_r = new unsigned char *[MV_NUMBER];
	for (int i = 0; i < MV_NUMBER; i++)
		cp_block_r[i] = new unsigned char[MV_NUMBER];
	unsigned char ** cp_block_t = new unsigned char *[MV_NUMBER];
	for (int i = 0; i < MV_NUMBER; i++)
		cp_block_t[i] = new unsigned char[MV_NUMBER];


	for (int i = 0; i < height / MV_NUMBER; i++) // MV_NUMBER 마크로블록의 크기 결국 256/8 = 32
	{
		for (int j = 0; j < width / MV_NUMBER; j++) // MV_NUMBER 마크로블록의 크기 결국 256 / 8 = 32
		{
			x = j * MV_NUMBER; // 실제 이미지의 좌표 마크로 블록 첫번째 인덱스
			y = i * MV_NUMBER; // 실제 이미지의 좌표 마크로 블록 첫번째 인덱스
			int diff = INT_MAX; //최대값을 찾기 위한 인자 
			short temp_mv_x = 0;
			short temp_mv_y = 0;
			for (mv_y = -1 * P_VALUE; mv_y < P_VALUE; mv_y++) //P_VALUE 검색 범위
			{
				for (mv_x = -1 * P_VALUE; mv_x < P_VALUE; mv_x++)//P_VALUE 검색 범위
				{
					if (((y + mv_y) >= 0) && ((y + mv_y) <= (height - MV_NUMBER)) && ((x + mv_x) >= 0) && ((x + mv_x) <= (width - MV_NUMBER)))
					{
						int temp = 0;
						Div_ImgBlock_MV(r_frame, cp_block_r, y + mv_y, x + mv_x); // r의 8by8 레퍼런스 프레임의 마크로 블록만 움직임 
						Div_ImgBlock_MV(t_frame, cp_block_t, y, x); // t의 8 by 8; 타겟 프레임의 마크로 블록은 움직이지 않는다 
						for (int k = 0; k < MV_NUMBER; k++)
						{
							for (int l = 0; l < MV_NUMBER; l++)
							{
								temp += abs((int)cp_block_r[k][l] - (int)cp_block_t[k][l]);
							}
						}
						if (temp < diff)
						{
							diff = temp;
							temp_mv_x = mv_x;
							temp_mv_y = mv_y;
						}
					}
				}
			}
			
			moVec.push_back(temp_mv_y); // 최종 선택된 mv_y  -> 최종적으로 한 프레임의 1024개의 모션 백터가 찾아짐;
			moVec.push_back(temp_mv_x); // 최종 선택된 mv_x
		}
	}

	for (int i = 0; i < MV_NUMBER; i++)
		delete[] cp_block_t[i];
	delete[] cp_block_t;
	for (int i = 0; i < MV_NUMBER; i++)
		delete[] cp_block_r[i];
	delete[] cp_block_r;
}

void MC_predic(vector<short>& moVec, unsigned char ** r_frame, unsigned char ** t_frame, unsigned char** pred, unsigned char** difference,int width,int height) {
	// 마크로 블록을 모션백터 만큼 움직이면 타겟프레임와 유사해 진다구?!?

	short mv_x = 0, mv_y = 0, x = 0, y = 0;
	unsigned char ** p_prime = new unsigned char *[H];
	for (int i = 0; i < H; i++)
		p_prime[i] = new unsigned char[W];
	unsigned char ** cp_block_r = new unsigned char *[MV_NUMBER];
	for (int i = 0; i < MV_NUMBER; i++)
		cp_block_r[i] = new unsigned char[MV_NUMBER];
	unsigned char ** cp_block_t = new unsigned char *[MV_NUMBER];
	for (int i = 0; i < MV_NUMBER; i++)
		cp_block_t[i] = new unsigned char[MV_NUMBER];


	for (int index = 0; index < moVec.size(); index += 2) {
		int move_i = moVec[index];
		int move_j = moVec[index + 1];
		
	}


	for (int index = 0; index < moVec.size(); index+=2) {
		int move_i = moVec[index];
		int move_j = moVec[index+1];

		for (int mcr_y = 0; mcr_y < height / MV_NUMBER; mcr_y++) { // MV_NUMBER 마크로블록의 크기 결국 256/8 = 32
			for (int mcr_x = 0; mcr_x < width / MV_NUMBER; mcr_x++) { // MV_NUMBER 마크로블록의 크기 결국 256 / 8 = 32
				y = mcr_y * MV_NUMBER; // 실제 이미지의 좌표 마크로 블록 첫번째 인덱스
				x = mcr_x * MV_NUMBER; // 실제 이미지의 좌표 마크로 블록 첫번째 인덱스

				for (int i = 0; i < 0 +MV_NUMBER; i++) {
					for (int j = 0; j <0 + MV_NUMBER; j++) {
						if (((y + i + move_i) >= 0) && ((y + i + move_i) < (H)) && ((x + j + move_j) >= 0) && ((x + j + move_j) < (W))) {
							p_prime[y + i][x + j] = r_frame[y + i + move_i][x + j + move_j];
						}
						//else {
						//	p_prime[y + i][x + j] = 0;
						//}
					}
				}

			}
		}
	}

	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			pred[i][j] = p_prime[i][j];
		}
	}


	short temp = 0;
	short min = 9999;
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			temp = -(t_frame[i][j]- p_prime[i][j]);
			//cout << temp << "\t";
			if (temp > 255) difference[i][j] = 255;
			else if (temp < 0) difference[i][j]= 0;
			else difference[i][j] = temp;
			//difference[i][j] = temp;
			/*
			temp = t_frame[i][j] - temp_image[i][j];
			if (temp < min) {
				min = temp;
				cout << min << "\t";
			}*/
		}
	}

	for (int i = 0; i < H; i++)
		delete[] p_prime[i];
	delete[] p_prime;
	for (int i = 0; i < MV_NUMBER; i++)
		delete[] cp_block_t[i];
	delete[] cp_block_t;
	for (int i = 0; i < MV_NUMBER; i++)
		delete[] cp_block_r[i];
	delete[] cp_block_r;

}

//unsigned char 데이터를 MV_NUMBER*MV_NUMBER 블록에 (m,n)에서부터 MV_NUMBER, MV_NUMBER범위를 넣어주는 함수

void encoding() {
	string outName = "Image_Out";
	string outName1 = "I_Out";
	string outName2 = "P_Out";
	string header_;
	ifstream outputFile;

	bool check_Iframe = true;
	int count = 0;
	string outdata;

	char* header = new char[54];
	unsigned char* Fulldata = new unsigned char[W*H * 25];
	unsigned char** I_Frame = new unsigned char*[H];
	
	for (int i = 0; i < H; i++) {
		I_Frame[i] = new unsigned char[W];
	}
	get_Full_data(Fulldata, header);

	char_to_string(header, header_size, header_);
	outdata += header_;
	unsigned char** temp_Iframe_decode = new unsigned char*[H];
	unsigned char** outRGB = new unsigned char*[H];
	unsigned char** R_frame = new unsigned char*[H];
	for (int i = 0; i < H; i++) {
		R_frame[i] = new unsigned char[W];
		temp_Iframe_decode[i] = new unsigned char[W];
		outRGB[i] = new unsigned char[W * 3];
	}

	unsigned char** T_frame = new unsigned char*[H];
	unsigned char** diff = new unsigned char*[H];
	unsigned char** pred = new unsigned char*[H];
	unsigned char** diff_decode = new unsigned char*[H];
	//unsigned char** p_frame_decode = new unsigned char*[H];



	for (int i = 0; i < H; i++) {
		T_frame[i] = new unsigned char[W];
		diff[i] = new unsigned char[W];
		diff_decode[i] = new unsigned char[W];
		pred[i] = new unsigned char[W];
		//p_frame_decode[i] = new unsigned char[W];
	}


	unsigned char** temp_Iframe = new unsigned char*[H];
	for (int i = 0; i < H; i++) {
		temp_Iframe[i] = new unsigned char[W];
	}

	
	//main roof
	for (int image_Num = 1; image_Num <= 25; image_Num++) {
		
		int start_Index = (image_Num - 1)*W*H;
		int count_index = 0;
		string temp_i;
		string diff_jpeg;
		



		if (check_Iframe) { // Iframe 
			//if (image_Num ==1) { // Iframe	
			count = 0;
			check_Iframe = false;
			cout << "----------------Frame num: "<< image_Num <<"---------------------" <<endl;

			for (int i = 0; i < H; i++) {
				for (int j = 0; j < W; j++) {
					temp_Iframe[i][j] = Fulldata[start_Index + count_index];
					count_index++;
				}
			}

			temp_i=JPEG_encoding(temp_Iframe);
			cout << "@@현재 프레임의 압축률:" << (double)temp_i.size() / (double)(W*H);

			// 사이즈 먼저 저장
			int size[2] = { diff_jpeg.size() , 0 };
			string size2 = int_to_string(size,2);
			outdata += size2;
			// 데이터 저장
			outdata += temp_i;
			//
			JPEG_decoding(temp_i, temp_Iframe_decode);
			
			for (int i = 0; i < H; i++) {
				for (int j = 0; j < W; j++) {
					R_frame[i][j] = temp_Iframe_decode[i][j];
				}
			}

			
			for (int i = 0; i < H; i++) {
				for (int j = 0,jj=0; j < W * 3;jj++, j+=3) {
					outRGB[i][j] = R_frame[i][jj];
					outRGB[i][j+1] = R_frame[i][jj];
					outRGB[i][j+2] = R_frame[i][jj];
				}
			}
			/*
			string buf = to_string(image_Num);
			if (image_Num < 10) buf = "0" + buf;
			ofstream outputFile;
			outputFile.open(OUT_FILE_PATH+ outName+buf+".bmp",ios::binary);
			outputFile.write((char*)header, header_size);
			for (int i = 0; i < H; i++) {
				outputFile.write((char*)outRGB[i], W * 3);
			}
			cout << "complete save " << image_Num << "th Iframe!" << endl;
			*/
		}

		else {					//Pframe
			/*
			decode_Iframe과 현재의 Pframe를 통해 MV를 구한다.
			*/
			cout << "----------------Frame num: " << image_Num << "---------------------" << endl;

			count++;
			vector<short> temp_moVec;
			string movString;
			

			for (int i = 0; i < H; i++) {
				for (int j = 0; j < W; j++) {
					T_frame[i][j] = Fulldata[start_Index + count_index];
					
					//cout << (int)R_frame[i][j] << "\t";
					count_index++;
				}
			}

			Get_MotionVector(temp_moVec, R_frame, T_frame, W, H);
			movString = short_to_string(temp_moVec, temp_moVec.size());


			MC_predic(temp_moVec, R_frame, T_frame, pred ,diff,W,H);
			// short to sting 으로 모션 백터 저장

			diff_jpeg = JPEG_encoding(diff);



			// 사이즈 저장
			int size[2] = { diff_jpeg.size()+movString.size(),movString.size() };
			string size2 = int_to_string(size, 2);
			outdata += size2;
			// 모션벡터 저장
			outdata += movString;
			// 데이터 저장
			outdata += diff_jpeg;


			cout << "@@현재 프레임의 압축률:" << (double)diff_jpeg.size() / (double)(W*H);
			JPEG_decoding(diff_jpeg, diff_decode);
			
			int tempdata = 0;
			for (int i = 0; i < H; i++) {
				for (int j = 0; j < W; j++) {
					R_frame[i][j] = pred[i][j] - diff_decode[i][j];
				}
			}
			
			// 저장
			for (int i = 0; i < H; i++) {
				for (int j = 0, jj = 0; j < W * 3; jj++, j += 3) {
					outRGB[i][j] = diff[i][jj];
					outRGB[i][j + 1] = diff[i][jj];
					outRGB[i][j + 2] = diff[i][jj];
				}
			}
			
			string aa = "diff";
			string buf2 = to_string(image_Num);
			if (image_Num < 10) buf2 = "0" + buf2;
			ofstream outputFile2;
			outputFile2.open(OUT_FILE_PATH + aa + buf2 + ".bmp", ios::binary);
			outputFile2.write((char*)header, header_size);
			for (int i = 0; i < H; i++) {
				outputFile2.write((char*)outRGB[i], W * 3);
			}
			
		

			if (count == 3) { 
				check_Iframe = true;
			}
		}
		
	}
	cout << "전체 압축률: " <<  (double)outdata.size()/(double)(25 * H*W)  << endl;
	copy2 = outdata;

	for (int i = 0; i < H; i++) {
		delete[] I_Frame[i];
		delete[] temp_Iframe_decode[i];
		delete[] outRGB[i];
		delete[] R_frame[i];
		delete[] T_frame[i];
		delete[] diff[i];
		delete[] pred[i];
		delete[] diff_decode[i];
	}
	
	delete[] I_Frame;
	delete[] temp_Iframe_decode;
	delete[] outRGB;
	delete[] R_frame;
	delete[] T_frame;
	delete[] diff;
	delete[] I_Frame;
	delete[] pred;
	delete[] diff_decode;
	delete[] Fulldata;
}

//---------------------------------------VIDEO_DECODING---------------------------------------------------------

void decoding() {
	string outName = "Image_Out";
	string outName1 = "I_Out";
	string outName2 = "P_Out";
	ifstream outputFile;

	bool check_Iframe = true;
	int count = 0;
	string inputdata = copy2;

	char* header = new char[54];
	
	for (int i = 0; i < 54; i++) {
		int index = i * 8;
		string temp_string = inputdata.substr(i,8);
		header[i] = String_to_char(temp_string);
	}

	unsigned char** temp_Iframe_decode = new unsigned char*[H];
	unsigned char** outRGB = new unsigned char*[H];
	unsigned char** R_frame = new unsigned char*[H];	
	unsigned char** P_frame = new unsigned char*[H];

	for (int i = 0; i < H; i++) {
		P_frame[i] = new unsigned char[W];
		R_frame[i] = new unsigned char[W];
		temp_Iframe_decode[i] = new unsigned char[W];
		outRGB[i] = new unsigned char[W * 3];
	}

	unsigned char** T_frame = new unsigned char*[H];
	unsigned char** diff = new unsigned char*[H];
	unsigned char** pred = new unsigned char*[H];
	unsigned char** diff_decode = new unsigned char*[H];
	//unsigned char** p_frame_decode = new unsigned char*[H];



	for (int i = 0; i < H; i++) {
		T_frame[i] = new unsigned char[W];
		diff[i] = new unsigned char[W];
		diff_decode[i] = new unsigned char[W];
		pred[i] = new unsigned char[W];
		//p_frame_decode[i] = new unsigned char[W];
	}


	unsigned char** temp_Iframe = new unsigned char*[H];
	for (int i = 0; i < H; i++) {
		temp_Iframe[i] = new unsigned char[W];
	}

	int start_Index = header_size;


	//main roof
	for (int image_Num = 1; image_Num <= 25; image_Num++) {
		int * sizeInfo = new int[2];
		int datalen;
		int dataSize;
		int vecSize;
		string dataSize_;
		string vecSize_;
		vector<short> temp_moVec;
		

		//데이터 정보 읽어오기
		dataSize_ = inputdata.substr(start_Index, sizeof(int) * 8);
		start_Index += sizeof(int) * 8;
		dataSize = String_to_int(dataSize_);

		vecSize_ = inputdata.substr(start_Index, sizeof(int) * 8);
		start_Index += sizeof(int) * 8;
		vecSize = String_to_int(vecSize_);
		
		datalen = dataSize - vecSize;




		

		
		

		if (check_Iframe) { // Iframe 

			string vecdata = inputdata.substr(start_Index, vecSize);
			start_Index += vecSize;

			string getdata = inputdata.substr(start_Index, datalen);
			start_Index += datalen;

			JPEG_decoding(getdata, R_frame);
			
			string buf = to_string(image_Num);
			if (image_Num < 10) buf = "0" + buf;
			ofstream outputFile_de;
			outputFile_de.open(OUT_FILE_PATH+ outName+buf+".bmp",ios::binary);
			outputFile_de.write((char*)header, header_size);
			for (int i = 0; i < H; i++) {
				outputFile_de.write((char*)outRGB[i], W * 3);
			}
			cout << "complete save " << image_Num << "th Iframe!" << endl;
			
		}

		else {					//Pframe
			count++;
			vector<short> temp_moVec;
			string vecdata = inputdata.substr(start_Index, vecSize);
			start_Index += vecSize;

			for (int i = 0; i < (W / 8)*(H / 8); i++) {
				int vecIndex = i * sizeof(short) * 8;
				string temp_mo = vecdata.substr(vecIndex, sizeof(short) * 8);
				temp_moVec[i] = String_to_short(temp_mo);
			}


			string getdata = inputdata.substr(start_Index, datalen);
			start_Index += datalen;

			JPEG_decoding(getdata,P_frame)

			for (int i = 0; i < H; i++) {
				for (int j = 0; j < W; j++) {
					T_frame[i][j] = Fulldata[start_Index + count_index];

					//cout << (int)R_frame[i][j] << "\t";
					count_index++;
				}
			}

			Get_MotionVector(temp_moVec, R_frame, T_frame, W, H);
			MC_predic(temp_moVec, R_frame, T_frame, pred, diff, W, H);


			diff_jpeg = JPEG_encoding(diff); //  이거 저장@@@
											 // 사이즈 저장
			int size[1] = { diff_jpeg.size() };
			string size2 = int_to_string(size, 1);
			outdata += size2;
			// 데이터 저장
			outdata += diff_jpeg;


			cout << "@@현재 프레임의 압축률:" << (double)diff_jpeg.size() / (double)(W*H);
			JPEG_decoding(diff_jpeg, diff_decode);

			int tempdata = 0;
			for (int i = 0; i < H; i++) {
				for (int j = 0; j < W; j++) {
					R_frame[i][j] = pred[i][j] - diff_decode[i][j];
				}
			}

			// 저장
			for (int i = 0; i < H; i++) {
				for (int j = 0, jj = 0; j < W * 3; jj++, j += 3) {
					outRGB[i][j] = diff[i][jj];
					outRGB[i][j + 1] = diff[i][jj];
					outRGB[i][j + 2] = diff[i][jj];
				}
			}

			string aa = "diff";
			string buf2 = to_string(image_Num);
			if (image_Num < 10) buf2 = "0" + buf2;
			ofstream outputFile2;
			outputFile2.open(OUT_FILE_PATH + aa + buf2 + ".bmp", ios::binary);
			outputFile2.write((char*)header, header_size);
			for (int i = 0; i < H; i++) {
				outputFile2.write((char*)outRGB[i], W * 3);
			}

			if (count == 3) {
				check_Iframe = true;
			}
		}

	}


	cout << "전체 압축률: " << (double)outdata.size() / (double)(25 * H*W) << endl;

	
	for (int i = 0; i < H; i++) {
		delete[] I_Frame[i];
		delete[] temp_Iframe_decode[i];
		delete[] outRGB[i];
		delete[] R_frame[i];
		delete[] T_frame[i];
		delete[] diff[i];
		delete[] pred[i];
		delete[] diff_decode[i];
	}

	delete[] I_Frame;
	delete[] temp_Iframe_decode;
	delete[] outRGB;
	delete[] R_frame;
	delete[] T_frame;
	delete[] diff;
	delete[] I_Frame;
	delete[] pred;
	delete[] diff_decode;
	delete[] Fulldata;
	cout << "COMPLETE ENCODING" << endl;
}



int main()
{
	encoding();
	decoding();

	return 0;
}


