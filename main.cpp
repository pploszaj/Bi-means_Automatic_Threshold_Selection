#include <iostream>

class ThresholdSelection {
public:
    int numRows;
    int numCols;
    int minVal;
    int maxVal;
    int BiGaussThrVal;
    int histHeight;
    int maxHeight;
    int* histAry;
    int* GaussAry;
    int* bestFitGaussAry;
    char** Graph

    ThresholdSelection(int rows, int cols, int min, int max) : numRows(rows), numCols(cols),minVal(min), maxVal(max),BiGaussThrVal(0), histHeight(0), maxHeight(0){
        histAry = new int[maxVal + 1];
        for (int i = 0; i <= maxVal; ++i) histAry[i] = 0;

        gaussAry = new int[maxVal + 1];
        for (int i = 0; i <= maxVal; ++i) gaussAry[i] = 0;

        bestFitGaussAry = new int[maxVal + 1];
        for (int i = 0; i <= maxVal; ++i) bestFitGaussAry[i] = 0;

        graph = new char*[maxVal + 1];
        for (int i = 0; i <= maxVal; ++i) {
            graph[i] = new char[histHeight + 1];
            for (int j = 0; j <= histHeight; ++j) {
                graph[i][j] = ' ';
            }
        }
    }

    int loadHist(const std::string& inFile) {
        std::ifstream input(inFile);
        if(!input.is_open()){
            std::cerr < "Error opening input file: " << inFile << std::endl;
            return -1;
        }

        int maxHistValue = 0;
        int value;
        int index = 0;
        while (input >> value) {
            if (index > maxVal) {
                std::cerr << "Histogram data exceeds expected range." << std::endl;
                break;
            }
            histAry[index++] = value;
            if (value > maxHistValue) {
                maxHistValue = value;
            }
        }
        histHeight = maxHistValue; // Update histHeight to reflect the loaded histogram's maximum value
        input.close();
        return maxHistValue;
    }


};