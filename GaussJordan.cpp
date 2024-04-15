#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <algorithm>

using namespace std;

// Đọc ma trận từ bàn phím
void readMatrix(vector<vector<double>>& matrix) {
    int rows, cols;
    cin >> rows >> cols;  // Đọc số hàng và số cột của ma trận bo sung
    matrix.resize(rows, vector<double>(cols));
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            cin >> matrix[i][j];  // Đọc giá trị từng phần tử của ma trận
        }
    }
}

// In ma trận ra màn hình
void printMatrix(const vector<vector<double>>& matrix) {
    for (const vector<double>& row : matrix) {
        for (double element : row) {
            cout << setw(8) << element << " ";
        }
        cout << endl;
    }
}

// Tìm phần tử giải và vị trí của nó
void findPivot(const vector<vector<double>>& matrix, vector<int>& index_row, vector<int>& index_column) {
    double pivot_element = 0.0;
    int row_pivot_element = -1;
    int col_pivot_element = -1;
    for (int row = 0; row < matrix.size(); row++) {
        if (find(index_row.begin(), index_row.end(), row) != index_row.end()) {
            continue;  // Bỏ qua hàng đã có phần tử giải
        }
        for (int col = 0; col < matrix[row].size() - 1; col++) {
            if (find(index_column.begin(), index_column.end(), col) != index_column.end()) {
                continue;  // Bỏ qua cột đã có phần tử giải
            }
            double abs_value = abs(matrix[row][col]);
            if (abs_value > pivot_element) {
                pivot_element = abs_value;
                row_pivot_element = row;
                col_pivot_element = col;
            }
        }
    }
    if (row_pivot_element != -1 && col_pivot_element != -1) {
        index_row.push_back(row_pivot_element);
        index_column.push_back(col_pivot_element);
        cout << "Phần tử giải: " << matrix[row_pivot_element][col_pivot_element] << endl;
        cout << "Vị trí: (" << row_pivot_element + 1 << ", " << col_pivot_element + 1 << ")" << endl;
        cout << endl;
    }
}

// Khử ma trận bằng phương pháp Gauss-Jordan
void gaussJordan(vector<vector<double>>& matrix, vector<int>& index_row, vector<int>& index_column) {
    findPivot(matrix, index_row, index_column);
    if (index_row.size() > 0 && index_column.size() > 0) {
        int pivot_row = index_row.back();
        int pivot_col = index_column.back();
        double pivot_element = matrix[pivot_row][pivot_col];
        for (int col = 0; col < matrix[pivot_row].size(); col++) {
            matrix[pivot_row][col] /= pivot_element;
        }
        for (int row = 0; row < matrix.size(); row++) {
            if (row != pivot_row) {
                double multiplier = matrix[row][pivot_col];
                for (int col = 0; col < matrix[row].size(); col++) {
                    matrix[row][col] -= multiplier * matrix[pivot_row][col];
                }
            }
        }
    }
}

// Chuẩn hóa ma trận
void normalizeMatrix(vector<vector<double>>& matrix, const vector<int>& index_row, const vector<int>& index_column) {
    for (int i = 0; i < index_row.size(); i++) {
        int row = index_row[i];
        int col = index_column[i];
        matrix[row][col] = 1.0;
    }
    cout << "Chuẩn hóa ma trận:" << endl;
    printMatrix(matrix);
    cout << "===========================" << endl;
}

// Hiển thị nghiệm của hệ phương trình

void displaySolutions(const vector<vector<double>>& matrix, const vector<int>& index_row, const vector<int>& index_column) {
    vector<vector<double>> result(matrix[0].size() - 1, vector<double>(matrix[0].size(), 0.0));
    for (int col = 0; col < matrix[0].size() - 1; col++) {
        if (find(index_column.begin(), index_column.end(), col) != index_column.end()) {
            int index = distance(index_column.begin(), find(index_column.begin(), index_column.end(), col));
            int pivot_row = index_row[index];
            result[col][0] = matrix[pivot_row].back();
            for (int i = 0; i < matrix[0].size() - 1; i++) {
                if (find(index_column.begin(), index_column.end(), i) == index_column.end()) {
                    result[col][i + 1] = -matrix[pivot_row][i];
                }
            }
        }
        else {
            result[col][col + 1] = 1.0;
        }
    }
    cout << "Nghiệm của hệ phương trình:" << endl;
    for (int i = 0; i < result.size(); i++) {
        cout << "X" << i + 1 << " = ";
        for (int k = 0; k < result[0].size(); k++) {
            if (k == 0) {
                cout << result[i][k];
            }
            else {
                cout << " + " << result[i][k] << ".X" << k;
            }
        }
        cout << endl;
    }
}

// Kiểm tra hạng của hệ phương trình
void checkRank(const vector<vector<double>>& matrix, const vector<int>& index_row, const vector<int>& index_column) {
    int rank1 = 0;
    int rank2 = 0;
    for (const vector<double>& row : matrix) {
        if (any_of(row.begin(), row.end() - 1, [](double element) { return element != 0.0; })) {
            rank1++;
        }
        if (any_of(row.begin(), row.end(), [](double element) { return element != 0.0; })) {
            rank2++;
        }
    }
    if (rank1 < rank2) {
        cout << "Hệ phương trình vô nghiệm" << endl;
    }
    else if (rank1 < matrix[0].size() - 1) {
        cout << "Hệ phương trình có vô số nghiệm" << endl;
        displaySolutions(matrix, index_row, index_column);
    }
    else {
        cout << "Hệ phương trình có nghiệm duy nhất" << endl;
        displaySolutions(matrix, index_row, index_column);
    }
}


int main() {
    vector<vector<double>> matrix;
    readMatrix(matrix);
    
    cout << "Ma trận mở rộng vừa nhập là:" << endl;
    printMatrix(matrix);
    cout << "===========================" << endl;
    
    vector<int> index_row;
    vector<int> index_column;
    
    while (index_row.size() < matrix.size() && index_column.size() < matrix[0].size() - 1) {
        gaussJordan(matrix, index_row, index_column);
    }
    
    normalizeMatrix(matrix, index_row, index_column);
    checkRank(matrix, index_row, index_column);
    
    cout << "===========================" << endl;
    
    return 0;
}
