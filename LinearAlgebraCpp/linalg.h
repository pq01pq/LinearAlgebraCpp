#pragma once
#include <iostream>
#include <string>
#include <stdexcept>

namespace linalg {
	class Matrixx;
	class Roww;
	class Vectorr;
	class MatrixxCell;

	class Matrixx {
		friend class Roww;
		friend class Vectorr;
		friend class MatrixxCell;
	public:
		Matrixx(int rowLength, int colLength);
		Matrixx(const Matrixx& copyMatrix);
		~Matrixx();

	private:
		Roww* rows;
		int rowLength, colLength;
	};

	class Roww {
		friend class Matrixx;
		friend class Vectorr;
		friend class MatrixxCell;
	public:
		Roww() = default;
		Roww(int colLength);
		~Roww();

	private:
		MatrixxCell* cells;
		int colLength;
	};

	class MatrixxCell {
		friend class Matrixx;
		friend class Roww;
		friend class Vectorr;
	public:

	private:
		double value;
	};
}