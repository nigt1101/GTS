#include <iostream>
#include<fstream>
#include <vector>

using namespace std;

void gauss(vector<vector<double>>& A, vector<vector<double>>& B, vector<vector<double>>& x) {
    int m = A.size();
    int n = A[0].size();
    int p = B[0].size();
    x.resize(m, vector<double>(p));

    for (int i = 0; i < n; i++) {
        // Chọn hàng i làm hàng chính
        for (int j = i + 1; j < n; j++) {
            // Trừ đi hàng i nhân với hệ số thích hợp để loại bỏ phần tử A[j][i]
            double t = A[j][i] / A[i][i];
            for (int k = 0; k < n; k++) {
                A[j][k] -= t * A[i][k];
            }
            for (int l = 0; l < p; l++) {
                B[j][l] -= t * B[i][l];
            }
        }
    }

    // Giải hệ phương trình tam giác trên
    for (int i = n - 1; i >= 0; i--) {
        for (int l = 0; l < p; l++) {
            x[i][l] = B[i][l];
        }
        for (int j = i + 1; j < n; j++) {
            for (int l = 0; l < p; l++) {
                x[i][l] -= A[i][j] * x[j][l];
            }
        }
        for (int l = 0; l < p; l++) {
            x[i][l] /= A[i][i];
        }
    }
}

int main() {
    //mo file
	ifstream file;
	/*chu y thay doi dia chi cua file
	co the tim dia chi cua file bang cach vao properties cua file roi copy phan location
	thay doi '\' thanh '\\'
	*/
	file.open("C:\\Users\\dell\\Downloads\\4_Seidel.txt"); 
	if(file.fail() == true){
		cout << "\nKhong ton tai file";	
	}

    int m, n, p;
    file >> m;
    file >> n;
    file >> p;

    // Nhập ma trận A
    vector<vector<double>> A(m, vector<double>(n));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            file >> A[i][j];
        }
    }

    // Nhập ma trận B
    vector<vector<double>> B(n, vector<double>(p));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < p; j++) {
            file >> B[i][j];
        }
    }

    // Giải hệ phương trình
    vector<vector<double>> x;
    gauss(A, B, x);

    // In ra nghiệm
    for (int i = 0; i < m; i++) {
        for (int l = 0; l < p; l++) {
            cout << "x" << i + 1 << " = " << x[i][l] << endl;
        }
    }

    return 0;
}
