#include "linalg.h"

using namespace std;
using namespace linalg;

class In;
class Out;

class In {
public:
	In() = default;
	~In() {
		
	}
	void setParent(Out& out) {
		parent = &out;
	}
	Out* parent;
};

class Out {
public:
	Out()
	{
		in = new In();
		in->setParent(*this);
	}
	~Out() {}
	void removeChild() {
		delete in;
	}
	In* in;
};



int main()
{
	/*int a[3][3] = {
		{1, 1, 1},
		{1, 1, 1},
		{1, 1, 1}
	};
	cout << a[0] << endl;
	cout << &(a[0][0]) << endl;
	cout << &(a[0][1]) << endl;
	cout << &(a[0][2]) << endl;
	cout << endl;
	cout << a[1] << endl;
	cout << &(a[1][0]) << endl;
	cout << &(a[1][1]) << endl;
	cout << &(a[1][2]) << endl;
	cout << endl;
	cout << a[2] << endl;
	cout << &(a[2][0]) << endl;
	cout << &(a[2][1]) << endl;
	cout << &(a[2][2]) << endl;*/

	/*Out out;
	cout << &out << endl;
	cout << out.in->parent << endl;
	out.removeChild();
	cout << &out << endl;*/
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