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

	Matrixx& Matrixx::operator=(const Matrixx& rightMatrix)
	{
		if (this == &rightMatrix) {
			return *this;
		}

		Matrixx copyMatrix(rightMatrix);
		swap(*this, copyMatrix);
		return *this;
	}
	void swap(Matrixx& leftMatrix, Matrixx& rightMatrix) noexcept
	{
		std::swap(leftMatrix.rows, rightMatrix.rows);
		std::swap(leftMatrix.height, rightMatrix.height);
		std::swap(leftMatrix.width, rightMatrix.width);
	}
	
	
	const int Matrixx::getHeight() const
	{
		return height;
	}
	const int Matrixx::getWidth() const
	{
		return width;
	}

	Roww& Matrixx::getRow(int row) const
	{
		Roww copyRow(row);
		for (int col = 0; col < width; col++) {
			copyRow[col] = rows[row][col];
		}
		return copyRow;
	}
	Vectorr& Matrixx::getColumn(int col) const
	{
		Vectorr copyVector(width);
		for (int row = 0; row < height; row++) {
			copyVector[row] = rows[row][col];
		}
		return copyVector;
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
		entries = new double[width];
	}
	Roww::Roww(const Roww& copyRow)
		: Roww(copyRow.getWidth())
	{
		for (int col = 0; col < width; col++) {
			entries[col] = copyRow[col];
		}
	}
	linalg::Roww::~Roww()
	{
		delete[] entries;
	}
	void Roww::init(int width)
	{
		this->width = width;
		entries = new double[width];
	}

	double& Roww::operator[](int col)
	{
		return entries[col];
	}
	const double& Roww::operator[](int col) const
	{
		return entries[col];
	}

	Allocator& Roww::operator<<(const double value)
	{
		this->allocate(0, value);
		return *(new Allocator(*this, 1));
	}
	void Roww::allocate(const int sequence, const double value)
	{
		entries[sequence] = value;
	}

	Roww& Roww::operator=(const Roww& rightRow)
	{
		if (this == &rightRow) {
			return *this;
		}

		Roww copyRow(rightRow);
		swap(*this, copyRow);
		return *this;
	}
	void swap(Roww& leftRow, Roww& rightRow) noexcept
	{
		std::swap(leftRow.entries, rightRow.entries);
		std::swap(leftRow.width, rightRow.width);
	}

	const int Roww::getWidth() const
	{
		return width;
	}

	const std::string Roww::str() const
	{
		std::string rowString = "[\t";
		for (int col = 0; col < width; col++) {
			rowString += to_string(entries[col]) + "\t";
		}
		rowString += "]\n";
		return rowString;
	}

	
	
	






	linalg::Vectorr::Vectorr(int height)
		: height(height)
	{
		entries = new double[height];
	}
	Vectorr::Vectorr(const Vectorr& copyVector)
		: Vectorr(copyVector.getHeight())
	{
		for (int row = 0; row < height; row++) {
			entries[row] = copyVector[row];
		}
	}
	linalg::Vectorr::~Vectorr()
	{
		delete[] entries;
	}
	double& Vectorr::operator[](int row)
	{
		return entries[row];
	}
	const double& Vectorr::operator[](int row) const
	{
		return entries[row];
	}
	Allocator& Vectorr::operator<<(const double value)
	{
		this->allocate(0, value);
		return *(new Allocator(*this, 1));
	}
	void Vectorr::allocate(const int sequence, const double value)
	{
		entries[sequence] = value;
	}

	Vectorr& Vectorr::operator=(const Vectorr& rightVector)
	{
		if (this == &rightVector) {
			return *this;
		}

		Vectorr copyRow(rightVector);
		swap(*this, copyRow);
		return *this;
	}
	void swap(Vectorr& leftVector, Vectorr& rightVector) noexcept
	{
		std::swap(leftVector.entries, rightVector.entries);
		std::swap(leftVector.height, rightVector.height);
	}

	const int Vectorr::getHeight() const
	{
		return height;
	}

	const std::string Vectorr::str() const
	{
		std::string vectorString = "(" + to_string(height) + " row Vector)\n";
		for (int row = 0; row < height; row++) {
			vectorString += "[\t" + to_string(entries[row]) + "\t]\n";
		}
		return vectorString;
	}
	
	void Vectorr::init(int height)
	{
		this->height = height;
		entries = new double[height];
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