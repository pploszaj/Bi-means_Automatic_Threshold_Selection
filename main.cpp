#include <iostream>
#include <fstream>
#include <string>

class ThresholdSelection {
public:
    int numRows, numCols, minVal, maxVal;
    int biGaussThrVal; // the auto selected threshold value by the Bi-Gaussian method
    int histHeight; // The largest hist[i] in the input histogram
    int maxHeight; // The largest hist[i] within a given range of the histogram. Initialize to 0.
    int* histAry; // a 1D integer array (size of maxVal + 1) to store the histogram. It needs to be dynamically allocated at run time; initialize to zero.
    int* gaussAry;  // a 1D integer array (size of maxVal + 1) to store the “modified” Gaussian curve values. It needs to be dynamically allocated at run time. initialize to zero.
    int* bestFitGaussAry; // to store the best biGaussian curves. Initialize to zero.
    char** graph; // a 2-D char array size of maxVal+1 by histHeight+1, initialize to blank

    ThresholdSelection(int rows, int cols, int min, int max) : numRows(rows), numCols(cols),minVal(min), maxVal(max),biGaussThrVal(0), histHeight(0), maxHeight(0){
//        histAry = new int[maxVal + 1];
//        for (int i = 0; i < maxVal; i++) histAry[i] = 0;


        bestFitGaussAry = new int[maxVal + 1];
        for (int i = 0; i < maxVal; i++) bestFitGaussAry[i] = 0;

        graph = new char*[maxVal + 1];
        for (int i = 0; i < maxVal; i++) {
            graph[i] = new char[histHeight + 1];
            for (int j = 0; j < histHeight; j++) {
                graph[i][j] = ' ';
            }
        }
    }

    int loadHist(int* ary, std::ifstream& inFile) {
        if(!inFile){
            std::cerr << "File is not open" << std::endl;
            return -1;
        }

        int maxHistValue = 0;
        int histVal, count;
        while (inFile >> histVal >> count) {
            if (histVal < 0 || histVal > maxVal) {
                std::cerr << "Histogram value " << histVal << " is out of range." << std::endl;
                break;
            }
            ary[histVal] = count;
            if (count > maxHistValue) {
                maxHistValue = count;
            }
        }
        histHeight = maxHistValue; // Update histHeight
        return maxHistValue;
    }

    void dispHist(std::ofstream& outFile) {
        outFile << numRows << " " << numCols << " " << minVal << " " << maxVal << std::endl;
        for (int i = minVal; i <= maxVal; ++i) {
            outFile << i << " (" << histAry[i] << "):";
            for (int j = 0; j < histAry[i]; ++j) {
                outFile << "+";
            }
            outFile << std::endl;
        }
    }

    void copyArys(int* ary1, int* ary2, int size) {
        for (int i = 0; i < size; ++i) {
            ary2[i] = ary1[i];
        }
    }

    void plotHist(){
        for (int i = 0; i <= maxVal; ++i) {
            int count = histAry[i];
            // Plot '+' for each count at the ith position, from bottom up
            for (int j = 0; j < count && j <= histHeight; ++j) {
                // Assuming the bottom of the graph is the 0th row, and we fill upwards
                graph[i][histHeight - j] = '+';
            }
        }
    }

    void setZero(int* ary, int size){
        for(let i = 0; i < size; i++){
            ary[i] = 0;
        }
    }

    int biGaussian(int* gaussAry, std::ofstream& debugFile){
        debugFile << "Entering biGaussian method" << std::endl;
        double sum1, sum2, total, minSumDiff;
        int offSet = (maxVal - minVal) / 10;
        int dividePt = offSet;
        int bestThr = dividePt;
        minSumDiff = 99999.0;

        while(dividePt < (maxVal - offSet)){
            setZero(gaussAry, maxVal + 1);
            sum1 = fitGauss(0, dividePt, debugFile);
            sum2 = fitGauss(dividePt, debugFile);
            total = sum1 + sum2;
            if(total < minSumDiff) {
                minSumDiff = total;
                bestThr = dividePt;
                copyArys(gaussAry, bestFitGaussAry);
            }

            debugFile << "In biGaussian (): dividePt = " << dividePt <<", sum1= "<< sum1 <<", sum2= "<< sum2 <<", total= "<< total <<", minSumDiff = "<< minSumDiff <<"and bestThr= "<< bestThr << std::endl;
            dividePt++;
        }

        debugFile << "leaving biGaussian method, minSumDiff = " << minSumDiff << "bestThr is " << bestThr << std::endl;
        return bestThr;
    }

    double fitGauss(int leftIndex, int rightIndex, int* histAry, int* gaussAry, int maxHeight, int** graph, std::ofstream& debugFile){
        debugFile << "“Entering fitGauss method" << std::endl;
        double mean, var, gVal;
        double sum = 0.0;
        mean = computeMean(leftIndex, rightIndex, mean, maxHeight, histAry, debugFile);
        var = computeVar(leftIndex, rightIndex, maxHeight, histAry, debugFile);
        int index = leftIndex;
        while(index <= rightIndex){
            gVal = modifiedGauss(index, mean, var, maxHeight);
            sum += abs(gVal - histAry[index]);
            gaussAry[index] = gVal;
            index++;
        }

        debugFile << "leaving fitGauss method, sum is: " << sum << std::endl;

        return sum;
    }

    double computeMean(){}

    double computeVar(){}

    void modifiedGauss(){}

    void printGraph(std::ofstream& deBugFile){
        if (graph == nullptr) {
            deBugFile << "Graph is not initialized." << std::endl;
            return;
        }

        deBugFile << "In main (), below is the Graph after plotting the histogram onto Graph" << std::endl;
        for (int j = histHeight; j >= 0; --j) {
            for (int i = 0; i <= maxVal; ++i) {
                deBugFile << graph[i][j];
            }
            deBugFile << std::endl;
        }
    }

};

int main(int argc, char* argv[]){
    if(argc < 5){
        std::cerr << "Incorrect number of arguments" << std::endl;
        return 1;
    }

    //step 0
    std::ifstream inFile1(argv[1]);
    std::ofstream outFile1(argv[2]);
    std::ofstream outFile2(argv[3]);
    std::ofstream deBugFile(argv[4]);

    if (!inFile1 || !outFile1 || !outFile2 || !deBugFile) {
        std::cerr << "File open error." << std::endl;
        return 1;
    }

    //step 1
    int numRows, numCols, minVal, maxVal;
    inFile1 >> numRows >> numCols >> minVal >> maxVal;

    ThresholdSelection ts(numRows, numCols, minVal, maxVal);

    int* histAry = new int[maxVal + 1];
    for (int i = 0; i < maxVal; i++) histAry[i] = 0;
    ts.histHeight = ts.loadHist(histAry, inFile1);

    //step 2
    outFile1 << "in main(), below is the input histogram" << std::endl;
    ts.dispHist(outFile1);

    //step 3
    ts.plotHist();
    ts.printGraph(deBugFile);

    //step 4
    int* gaussAry = new int[maxVal + 1];
    for (int i = 0; i < maxVal; i++) gaussAry[i] = 0;
    ts.biGaussThrVal = ts.biGaussian(gaussAry, deBugFile);

}