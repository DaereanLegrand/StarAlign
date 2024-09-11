#include <climits>
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>

using std::string;
using std::vector;

enum Direction { DIAG, UP, LEFT };

void 
storeBRCA1Sequences(const std::string& filename, std::vector<std::string>& forwardSeqs, std::vector<std::string>& reverseSeqs) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file: " + filename);
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.find("F:") != std::string::npos) {
            size_t pos = line.find("5′-");
            size_t end = line.find("-3′");
            if (pos != std::string::npos && end != std::string::npos && pos < end) {
                forwardSeqs.push_back(line.substr(pos + 4, end - pos - 4));
            }
        } else if (line.find("R:") != std::string::npos) {
            size_t pos = line.find("5′-");
            size_t end = line.find("-3′");
            if (pos != std::string::npos && end != std::string::npos && pos < end) {
                reverseSeqs.push_back(line.substr(pos + 4, end - pos - 4));
            }
        }
    }

    if (forwardSeqs.empty() && reverseSeqs.empty()) {
        throw std::runtime_error("No valid sequences found in file: " + filename);
    }
}

int 
getAlignmentCount(const std::string& a, const std::string& b, std::vector<std::vector<Direction>>& directionTable) {
    int aSize = a.length() + 1;
    int bSize = b.length() + 1;

    std::vector<std::vector<int>> scoreTable(aSize, std::vector<int>(bSize, 0));
    directionTable.resize(aSize, std::vector<Direction>(bSize, DIAG));

    for (int i = 0; i < aSize; i++) {
        scoreTable[i][0] = i * (-2);
        directionTable[i][0] = UP;
    }

    for (int j = 0; j < bSize; j++) {
        scoreTable[0][j] = j * (-2);
        directionTable[0][j] = LEFT;
    }

    int x = aSize - 1;
    int y = bSize - 1;

    for (int line = 1; line <= (x + y - 1); line++) { 
        int start_j = std::max(0, line - x); 
        int count = std::min(line, std::min((y - start_j), x)); 

        for(int k = 0; k < count; k++) {
            int i = std::min(x, line) - k - 1;
            int j = start_j + k;

            i += 1;
            j += 1;

            if (i == 0 || j == 0) continue;

            int matchScore = (a[i - 1] == b[j - 1]) ? 1 : -1;
            int diag = scoreTable[i - 1][j - 1] + matchScore;
            int up = scoreTable[i - 1][j] - 2;
            int left = scoreTable[i][j - 1] - 2;

            int maxScore = std::max({diag, up, left});
            scoreTable[i][j] = maxScore;

            if (maxScore == diag) {
                directionTable[i][j] = DIAG;
            } else if (maxScore == up) {
                directionTable[i][j] = UP;
            } else {
                directionTable[i][j] = LEFT;
            }
        }
    }

    return scoreTable[aSize - 1][bSize - 1];
}

void
traceback(const std::string& a, const std::string& b, const std::vector<std::vector<Direction>>& directionTable, string& alignedA, string& alignedB) {
    int i = a.length();
    int j = b.length();

    while (i > 0 || j > 0) {
        if (directionTable[i][j] == DIAG) {
            alignedA = a[i - 1] + alignedA;
            alignedB = b[j - 1] + alignedB;
            i--;
            j--;
        } else if (directionTable[i][j] == UP) {
            alignedA = a[i - 1] + alignedA;
            alignedB = '-' + alignedB;
            i--;
        } else { // LEFT
            alignedA = '-' + alignedA;
            alignedB = b[j - 1] + alignedB;
            j--;
        }
    }

    while (i > 0) {
        alignedA = a[i - 1] + alignedA;
        alignedB = '-' + alignedB;
        i--;
    }
    while (j > 0) {
        alignedA = '-' + alignedA;
        alignedB = b[j - 1] + alignedB;
        j--;
    }
}

void
starAlign(vector<string> sequences)
{
    vector<vector<int>> score(sequences.size(), vector<int>(sequences.size(), 0));
    vector<int> sum(sequences.size(), 0);

    for (int i = 0; i < sequences.size(); i++) {
        for (int j = 0; j < sequences.size(); j++) {
            if (i != j) {
                vector<vector<Direction>> directionTable;
                score[i][j] = getAlignmentCount(sequences[i], sequences[j], directionTable);
            } else {
                score[i][j] = 0;
            }
            sum[i] += score[i][j];
        }
    }

    for (int i = 0; i < sequences.size(); i++) {
        for (int j = 0; j < sequences.size(); j++) {
            printf("%3d\t", score[i][j]);
        }
        printf("\n");
    }

    printf("\n");

    int cmax = 0;
    for (int i = 1; i < sum.size(); i++) {
        cmax = (sum[i] > sum[cmax]) ? i : cmax;
    }

    printf("Sequence with maximum sum: %d %s\n", cmax, sequences[cmax].c_str());
    printf("\n");

    vector<string> aligned;
    aligned.push_back(sequences[cmax]);
    for (int i = 0; i < sequences.size(); i++) {
        if (i != cmax) {
            vector<vector<Direction>> directionTable;
            getAlignmentCount(sequences[cmax], sequences[i], directionTable);

            string alignedA, alignedB;
            traceback(sequences[cmax], sequences[i], directionTable, alignedA, alignedB);
            aligned.push_back(alignedB);

            printf("Alignment between %s and %s:\n", sequences[cmax].c_str(), sequences[i].c_str());
            printf("%s\n", sequences[cmax].c_str());
            printf("%s\n", alignedB.c_str());
        }
    }
    printf("\n");

    for (auto i: aligned) {
        printf("%s\n", i.c_str());
    }
}

int
main(int argc, char *argv[]) {
    vector<string> sequences = {
        "ATTGCCATT",
        "ATGGCCATT",
        "ATCCAATTTT",
        "ATCTTCTT",
        "ACTGACC"
    };

    starAlign(sequences);

    if (argc < 2) return 0;

    std::vector<std::string> forwardSequences;
    std::vector<std::string> reverseSequences;

    storeBRCA1Sequences(argv[1], forwardSequences, reverseSequences);

    starAlign(forwardSequences);
    starAlign(reverseSequences);

    return 0;
}
