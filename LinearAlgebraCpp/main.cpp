#include "linalg.h"

using namespace std;
using namespace linalg;

int main()
{
	try {
		Matrixx matrix1(3, 6);
		matrix1 <<
			0, 3, -6, 6, 4, -5,
			3, -7, 8, -5, 8, 9,
			3, -9, 12, -9, 6, 15;
		matrix1.reduce();
		cout << matrix1 << endl;

		Matrixx matrix2(3, 3);
		matrix2 <<
			1, 2, 0,
			0, 0, 1,
			1, 1, 0;
		Matrixx matrix2inverse = matrix2.inverse();
		cout << matrix2inverse << endl;

		cout << matrix2 * matrix2inverse << endl;
		cout << matrix2inverse * matrix2 << endl;
	}
	catch (const logic_error& e) {
		cout << e.what() << endl;
	}
}