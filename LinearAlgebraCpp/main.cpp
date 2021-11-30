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