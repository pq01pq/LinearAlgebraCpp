#include "linalg.h"

using namespace std;
using namespace linalg;

int main()
{
	try {
		Matrixx matrix(3, 6);
		matrix <<
			0, 3, -6, 6, 4, -5,
			3, -7, 8, -5, 8, 9,
			3, -9, 12, -9, 6, 15;
		matrix.reduce();
		cout << matrix << endl;

		Matrixx matrix2(3, 3);
		matrix2 <<
			1, 2, 0,
			0, 0, 1,
			1, 1, 0;
		Matrixx matrix3 = matrix2.inverse();

		cout << matrix2 * matrix3 << endl;
		cout << matrix3 * matrix2 << endl;

		Matrixx matrix4(4, 4);
		matrix4 <<
			1, 0, 0, -1,
			0, -1, 1, 0,
			-1, -1, 1, -1,
			0, 1, 1, -1;
		cout << matrix4.inverse() << endl;
	}
	catch (const logic_error& e) {
		cout << e.what() << endl;
	}
}