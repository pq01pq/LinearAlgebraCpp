#include "linalg.h"

using namespace std;
using namespace linalg;

int main()
{
	try {
		cout << "행 축약" << endl;
		Matrixx matrix1(3, 6);
		matrix1 <<
			0, 3, -6, 6, 4, -5,
			3, -7, 8, -5, 8, 9,
			3, -9, 12, -9, 6, 15;
		matrix1.reduce();
		cout << matrix1 << endl;

		cout << "역행렬" << endl;
		Matrixx matrix2(3, 3);
		matrix2 <<
			1, 2, 0,
			0, 0, 1,
			1, 1, 0;
		Matrixx matrix2inverse = matrix2.inverse();
		cout << matrix2inverse << endl;

		cout << matrix2 * matrix2inverse << endl;
		cout << matrix2inverse * matrix2 << endl;

		cout << "블록행렬 곱" << endl;
		Matrixx a(4, 4);
		a <<
			1, -7, 8, 4,
			6, 3, -5, 0,
			-2, -1, 0, 6,
			11, -4, 4, 1;
		Matrixx b(4, 3);
		b <<
			7, -2, 5,
			8, 12, 12,
			-3, -1, 4,
			0, 4, 1;

		Matrixx a00 = a.block(0, 0, 2, 2);	Matrixx a01 = a.block(0, 2, 2, 2);
		Matrixx a10 = a.block(2, 0, 2, 2);	Matrixx a11 = a.block(2, 2, 2, 2);

		Matrixx b00 = b.block(0, 0, 2, 2);	Matrixx b01 = b.block(0, 2, 2, 1);
		Matrixx b10 = b.block(2, 0, 2, 2);	Matrixx b11 = b.block(2, 2, 2, 1);
		
		// '&' == 좌우 append, '|' == 상하 append (우선순위 : & > | )
		Matrixx blockMultiply =
			(a00*b00 + a01*b10) & (a00*b01 + a01*b11) |
			(a10*b00 + a11*b10) & (a10*b01 + a11*b11);
		
		cout << a * b << endl;
		cout << blockMultiply << endl;
	}
	catch (const logic_error& e) {
		cout << e.what() << endl;
	}
}