#include "linalg.h"

using namespace std;

namespace linalg {
	Allocator::Allocator(Allocatable& target, const int sequence)
		: target(target), sequence(sequence)
	{
	}
	Allocator::~Allocator()
	{
	}
	Allocator& Allocator::operator,(const double value)
	{
		target.allocate(sequence, value);
		sequence++;
		return *this;
	}






	Matrixx::Matrixx(int height, int width)
		: height(height), width(width)
	{
		rows = new Roww[height];
		for (int row = 0; row < height; row++) {
			rows[row].init(width);
		}
	}
	Matrixx::Matrixx(const Matrixx& copyMatrix)
		: height(copyMatrix.getHeight()), width(copyMatrix.getWidth())
	{
		rows = new Roww[height];
		for (int row = 0; row < height; row++) {
			rows[row].init(width);
			for (int col = 0; col < width; col++) {
				rows[row][col] = copyMatrix[row][col];
			}
		}
	}
	Matrixx::~Matrixx()
	{
		delete[] rows;
	}

	Roww& Matrixx::operator[](int row)
	{
		return rows[row];
	}
	const Roww& Matrixx::operator[](int row) const
	{
		return rows[row];
	}

	double& Matrixx::operator()(int row, int col) const
	{
		return rows[row][col];
	}

	Allocator& Matrixx::operator<<(const double value)
	{
		this->allocate(0, value);
		return *(new Allocator(*this, 1));
	}
	void Matrixx::allocate(const int sequence, const double value)
	{
		rows[sequence / width][sequence % width] = value;
	}
	
	const int Matrixx::getHeight() const
	{
		return height;
	}
	const int Matrixx::getWidth() const
	{
		return width;
	}

	const std::string Matrixx::str() const
	{
		std::string matrixString = "(" + to_string(height) + " x " + to_string(width) + " Matrix)\n";
		for (int row = 0; row < height; row++) {
			matrixString += rows[row].str();
		}
		return matrixString;
	}






	Roww::Roww(int width)
		: width(width)
	{
		cells = new double[width];
	}
	Roww::Roww(const Roww& copyRow)
		: Roww(copyRow.getWidth())
	{
		for (int col = 0; col < width; col++) {
			cells[col] = copyRow[col];
		}
	}
	linalg::Roww::~Roww()
	{
		delete[] cells;
	}
	void Roww::init(int width)
	{
		this->width = width;
		cells = new double[width];
	}

	double& Roww::operator[](int col)
	{
		return cells[col];
	}
	const double& Roww::operator[](int col) const
	{
		return cells[col];
	}

	Allocator& Roww::operator<<(const double value)
	{
		this->allocate(0, value);
		return *(new Allocator(*this, 1));
	}
	void Roww::allocate(const int sequence, const double value)
	{
		cells[sequence] = value;
	}

	const int Roww::getWidth() const
	{
		return width;
	}

	const std::string Roww::str() const
	{
		std::string rowString = "[\t";
		for (int col = 0; col < width; col++) {
			rowString += to_string(cells[col]) + "\t";
		}
		rowString += "]\n";
		return rowString;
	}

	
	
	






	linalg::Vectorr::Vectorr(int height)
		: height(height)
	{
		cells = new double[height];
	}
	Vectorr::Vectorr(const Vectorr& copyVector)
		: Vectorr(copyVector.getHeight())
	{
		for (int row = 0; row < height; row++) {
			cells[row] = copyVector[row];
		}
	}
	linalg::Vectorr::~Vectorr()
	{
		delete[] cells;
	}
	double& Vectorr::operator[](int row)
	{
		return cells[row];
	}
	const double& Vectorr::operator[](int row) const
	{
		return cells[row];
	}
	Allocator& Vectorr::operator<<(const double value)
	{
		this->allocate(0, value);
		return *(new Allocator(*this, 1));
	}
	void Vectorr::allocate(const int sequence, const double value)
	{
		cells[sequence] = value;
	}

	const int Vectorr::getHeight() const
	{
		return height;
	}

	const std::string Vectorr::str() const
	{
		std::string vectorString = "(" + to_string(height) + " row Vector)\n";
		for (int row = 0; row < height; row++) {
			vectorString += "[\t" + to_string(cells[row]) + "\t]\n";
		}
		return vectorString;
	}
	
	void Vectorr::init(int height)
	{
		this->height = height;
		cells = new double[height];
	}






	//linalg::Celll::Celll(double value)
	//	: value(value)
	//{
	//}
	///*Celll::Celll(const Celll& copyCell)
	//	: value(copyCell.get())
	//{
	//}*/
	//linalg::Celll::~Celll()
	//{
	//}

	//Celll::operator double() const
	//{
	//	return value;
	//}

	///*Celll& Celll::operator=(const Celll& rightCell)
	//{
	//	if (this == &rightCell) {
	//		return *this;
	//	}

	//	value = rightCell;
	//	return *this;
	//}*/

	//void Celll::set(double value)
	//{
	//	this->value = value;
	//}
	//const double Celll::get() const
	//{
	//	return value;
	//}

}