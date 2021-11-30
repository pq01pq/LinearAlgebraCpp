#include "linalg.h"

using namespace std;
using namespace linalg;

int main()
{
	Matrixx matrix(3, 3);
	matrix <<
		1, 2, 3,
		4, 5, 6,
		7, 8, 9;
	cout << matrix;

	Vectorr vector(5);
	vector << 1, 2, 3, 4, 5;
	cout << vector;

	Matrixx matrix2 = matrix;
	cout << matrix2;

}