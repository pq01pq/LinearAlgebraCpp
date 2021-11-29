#include "linalg.h"

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
		: height(copyMatrix.height), width(copyMatrix.width)
	{
		rows = new Roww[height];
		for (int row = 0; row < height; row++) {
			rows[row].init(width);
			for (int col = 0; col < width; col++) {
				rows[row][col] = copyMatrix[row][col];
			}
		}
	}
	Matrixx::Matrixx(const Roww& copyRow)
		: height(1), width(copyRow.width)
	{
		rows = new Roww[1];
		rows[0].init(width);
		for (int col = 0; col < height; col++) {
			rows[0][col] = copyRow[col];
		}
	}
	Matrixx::Matrixx(const Vectorr& copyVector)
		: height(copyVector.height), width(1)
	{
		rows = new Roww[height];
		for (int row = 0; row < height; row++) {
			rows[row].init(1);
			rows[row][0] = copyVector[row];
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

	Matrixx Matrixx::operator+() const
	{
		return Matrixx(*this);
	}
	Matrixx Matrixx::operator-() const
	{
		Matrixx negativeMatrix(*this);
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < width; col++) {
				negativeMatrix[row][col] = preventNegativeZero(-rows[row][col]);
			}
		}
		return negativeMatrix;
	}

	Allocator& Matrixx::operator<<(const double value)
	{
		this->allocate(0, value);
		return *(new Allocator(*this, 1));
	}
	void Matrixx::allocate(const int sequence, const double value)
	{
		if (sequence < height * width) {
			rows[sequence / width][sequence % width] = preventNegativeZero(value);
		}
	}

	Matrixx& Matrixx::operator=(const Matrixx& rightMatrix)
	{
		if (this == &rightMatrix) {
			return *this;
		}

		Matrixx copyMatrix(rightMatrix);
		linalg::swap(*this, copyMatrix);
		return *this;
	}
	void swap(Matrixx& leftMatrix, Matrixx& rightMatrix) noexcept
	{
		std::swap(leftMatrix.rows, rightMatrix.rows);
		std::swap(leftMatrix.height, rightMatrix.height);
		std::swap(leftMatrix.width, rightMatrix.width);
	}

	Matrixx& Matrixx::operator+=(const Matrixx& rightMatrix)
	{
		if (height != rightMatrix.height) {
			std::string errorString =
				"Heights do not match : row" + std::to_string(height) + " + row" + std::to_string(rightMatrix.height) + "\n";
			try {
				rows[0] + rightMatrix[0];
			}
			catch (const std::logic_error& e) {
				errorString += e.what();
			}
			throw std::logic_error(errorString);
		}
		for (int row = 0; row < height; row++) {
			rows[row] += rightMatrix[row];
		}
		return *this;
	}
	Matrixx& Matrixx::operator-=(const Matrixx& rightMatrix)
	{
		if (height != rightMatrix.height) {
			std::string errorString =
				"Heights do not match : row" + std::to_string(height) + " - row" + std::to_string(rightMatrix.height) + "\n";
			try {
				rows[0] - rightMatrix[0];
			}
			catch (const std::logic_error& e) {
				errorString += e.what();
			}
			throw std::logic_error(errorString);
		}
		for (int row = 0; row < height; row++) {
			rows[row] -= rightMatrix[row];
		}
		return *this;
	}
	Matrixx& Matrixx::operator*=(const double multiplier)
	{
		for (int row = 0; row < height; row++) {
			rows[row] *= multiplier;
		}
		return *this;
	}

	Matrixx& Matrixx::operator*=(const Matrixx& rightMatrix)
	{
		if (width != rightMatrix.height) {
			throw std::logic_error("Cannot multiply matrices : ("
				+ std::to_string(height) + " x " + std::to_string(width) + ") x ("
				+ std::to_string(rightMatrix.height) + " x " + std::to_string(rightMatrix.width) + ")");
		}
		Matrixx resultMatrix(height, rightMatrix.width);
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < rightMatrix.width; col++) {
				double product = 0.0;
				for (int join = 0; join < width; join++) {
					product += rows[row][join] * rightMatrix[join][col];
				}
				resultMatrix[row][col] = preventNegativeZero(product);
			}
		}
		linalg::swap(*this, resultMatrix);
		return *this;
	}

	Matrixx& Matrixx::operator&=(const Matrixx& rightMatrix)
	{
		if (height != rightMatrix.height) {
			throw std::logic_error(
				"Heights do not match : row" + std::to_string(height) + " & row" + std::to_string(rightMatrix.height) + "\n");
		}
		Matrixx appendedMatrix(height, width + rightMatrix.width);
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < width; col++) {
				appendedMatrix[row][col] = rows[row][col];
			}
		}
		for (int row = 0; row < height; row++) {
			for (int col = width; col < width + rightMatrix.width; col++) {
				appendedMatrix[row][col] = rightMatrix[row][col - width];
			}
		}
		linalg::swap(*this, appendedMatrix);
		return *this;
	}
	Matrixx& Matrixx::operator&=(const Vectorr& rightVector)
	{
		if (height != rightVector.height) {
			throw std::logic_error(
				"Heights do not match : row" + std::to_string(height) + " & row" + std::to_string(rightVector.height) + "\n");
		}
		Matrixx appendedMatrix(height, width + 1);
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < width; col++) {
				appendedMatrix[row][col] = rows[row][col];
			}
		}
		for (int row = 0; row < height; row++) {
			appendedMatrix[row][width] = rightVector[row];
		}
		linalg::swap(*this, appendedMatrix);
		return *this;
	}
	Matrixx& Matrixx::operator|=(const Matrixx& lowerMatrix)
	{
		if (width != lowerMatrix.width) {
			throw std::logic_error(
				"Widths do not match : col" + std::to_string(width) + " | col" + std::to_string(lowerMatrix.width) + "\n");
		}
		Matrixx appendedMatrix(height + lowerMatrix.height, width);
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < width; col++) {
				appendedMatrix[row][col] = rows[row][col];
			}
		}
		for (int row = height; row < height + lowerMatrix.height; row++) {
			for (int col = 0; col < width; col++) {
				appendedMatrix[row][col] = lowerMatrix[row - height][col];
			}
		}
		linalg::swap(*this, appendedMatrix);
		return *this;
	}
	Matrixx& Matrixx::operator|=(const Roww& lowerRow)
	{
		if (width != lowerRow.width) {
			throw std::logic_error(
				"Widths do not match : col" + std::to_string(width) + " | col" + std::to_string(lowerRow.width) + "\n");
		}
		Matrixx appendedMatrix(height + 1, width);
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < width; col++) {
				appendedMatrix[row][col] = rows[row][col];
			}
		}
		for (int col = 0; col < width; col++) {
			appendedMatrix[height][col] = lowerRow[col];
		}
		linalg::swap(*this, appendedMatrix);
		return *this;
	}

	Matrixx operator+(const Matrixx& leftMatrix, const Matrixx& rightMatrix)
	{
		Matrixx resultMatrix(leftMatrix);
		resultMatrix += rightMatrix;
		return resultMatrix;
	}
	Matrixx operator-(const Matrixx& leftMatrix, const Matrixx& rightMatrix)
	{
		Matrixx resultMatrix(leftMatrix);
		resultMatrix -= rightMatrix;
		return resultMatrix;
	}
	Matrixx operator*(const double multiplier, const Matrixx& rightMatrix)
	{
		Matrixx resultMatrix(rightMatrix);
		resultMatrix *= multiplier;
		return resultMatrix;
	}
	Matrixx operator*(const Matrixx& leftMatrix, const double multiplier)
	{
		return operator*(multiplier, leftMatrix);
	}
	Matrixx operator*(const Matrixx& leftMatrix, const Matrixx& rightMatrix)
	{
		Matrixx resultMatrix(leftMatrix);
		resultMatrix *= rightMatrix;
		return resultMatrix;
	}

	Matrixx operator&(const Matrixx& leftMatrix, const Matrixx& rightMatrix)
	{
		Matrixx appendedMatrix(leftMatrix);
		appendedMatrix &= rightMatrix;
		return appendedMatrix;
	}
	Matrixx operator&(const Matrixx& leftMatrix, const Vectorr& rightVector)
	{
		Matrixx appendedMatrix(leftMatrix);
		appendedMatrix &= rightVector;
		return appendedMatrix;
	}
	Matrixx operator&(const Vectorr& leftVector, const Matrixx& rightMatrix)
	{
		Matrixx appendedMatrix(leftVector);
		appendedMatrix &= rightMatrix;
		return appendedMatrix;
	}
	Matrixx operator&(const Vectorr& leftVector, const Vectorr& rightVector)
	{
		Matrixx appendedMatrix(leftVector);
		appendedMatrix &= rightVector;
		return appendedMatrix;
	}

	Matrixx operator|(const Matrixx& upperMatrix, const Matrixx& lowerMatrix)
	{
		Matrixx appendedMatrix(upperMatrix);
		appendedMatrix |= lowerMatrix;
		return appendedMatrix;
	}
	Matrixx operator|(const Matrixx& upperMatrix, const Roww& lowerRow)
	{
		Matrixx appendedMatrix(upperMatrix);
		appendedMatrix |= lowerRow;
		return appendedMatrix;
	}
	Matrixx operator|(const Roww& upperRow, const Matrixx& lowerMatrix)
	{
		Matrixx appendedMatrix(upperRow);
		appendedMatrix |= lowerMatrix;
		return appendedMatrix;
	}
	Matrixx operator|(const Roww& upperRow, const Matrixx& lowerRow)
	{
		Matrixx appendedMatrix(upperRow);
		appendedMatrix |= lowerRow;
		return appendedMatrix;
	}
	

	const int Matrixx::getHeight() const
	{
		return height;
	}
	const int Matrixx::getWidth() const
	{
		return width;
	}

	Roww Matrixx::getRow(int row) const
	{
		Roww copyRow(row);
		for (int col = 0; col < width; col++) {
			copyRow[col] = rows[row][col];
		}
		return copyRow;
	}
	Vectorr Matrixx::getColumn(int col) const
	{
		Vectorr copyVector(width);
		for (int row = 0; row < height; row++) {
			copyVector[row] = rows[row][col];
		}
		return copyVector;
	}

	const std::string Matrixx::str() const
	{
		std::string matrixString = "(" + std::to_string(height) + " x " + std::to_string(width) + " Matrix)\n";
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
		: Roww(copyRow.width)
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

	Roww Roww::operator+() const
	{
		return Roww(*this);
	}
	Roww Roww::operator-() const
	{
		Roww negativeRow(*this);
		for (int col = 0; col < width; col++) {
			negativeRow[col] = preventNegativeZero(-entries[col]);
		}
		return negativeRow;
	}

	Allocator& Roww::operator<<(const double value)
	{
		this->allocate(0, value);
		return *(new Allocator(*this, 1));
	}
	void Roww::allocate(const int sequence, const double value)
	{
		if (sequence < width) {
			entries[sequence] = preventNegativeZero(value);
		}
	}

	Roww& Roww::operator=(const Roww& rightRow)
	{
		if (this == &rightRow) {
			return *this;
		}

		Roww copyRow(rightRow);
		linalg::swap(*this, copyRow);
		return *this;
	}
	void swap(Roww& leftRow, Roww& rightRow) noexcept
	{
		std::swap(leftRow.entries, rightRow.entries);
		std::swap(leftRow.width, rightRow.width);
	}

	Roww& Roww::operator+=(const Roww& rightRow)
	{
		if (width != rightRow.width) {
			throw std::logic_error(
				"Widths do not match : col" + std::to_string(width) + " + col" + std::to_string(rightRow.width) + "\n");
		}
		for (int col = 0; col < width; col++) {
			entries[col] = preventNegativeZero(entries[col] + rightRow[col]);
		}
		return *this;
	}
	Roww& Roww::operator-=(const Roww& rightRow)
	{
		if (width != rightRow.width) {
			throw std::logic_error(
				"Widths do not match : col" + std::to_string(width) + " - col" + std::to_string(rightRow.width) + "\n");
		}
		for (int col = 0; col < width; col++) {
			entries[col] = preventNegativeZero(entries[col] - rightRow[col]);
		}
		return *this;
	}
	Roww& Roww::operator*=(const double multiplier)
	{
		for (int col = 0; col < width; col++) {
			entries[col] = preventNegativeZero(multiplier * entries[col]);
		}
		return *this;
	}

	Roww operator+(const Roww& leftRow, const Roww& rightRow)
	{
		Roww resultRow(leftRow);
		resultRow += rightRow;
		return resultRow;
	}
	Roww operator-(const Roww& leftRow, const Roww& rightRow)
	{
		Roww resultRow(leftRow);
		resultRow -= rightRow;
		return resultRow;
	}
	Roww operator*(const double multiplier, const Roww& rightRow)
	{
		Roww resultRow(rightRow);
		resultRow *= multiplier;
		return resultRow;
	}
	Roww operator*(const Roww& leftRow, const double multiplier)
	{
		return operator*(multiplier, leftRow);
	}

	const int Roww::getWidth() const
	{
		return width;
	}

	const std::string Roww::str() const
	{
		std::string rowString = "[\t";
		for (int col = 0; col < width; col++) {
			rowString += std::to_string(entries[col]) + "\t";
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
		: Vectorr(copyVector.height)
	{
		for (int row = 0; row < height; row++) {
			entries[row] = copyVector[row];
		}
	}
	linalg::Vectorr::~Vectorr()
	{
		delete[] entries;
	}
	void Vectorr::init(int height)
	{
		this->height = height;
		entries = new double[height];
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
		if (sequence < height) {
			entries[sequence] = preventNegativeZero(value);
		}
	}

	Vectorr Vectorr::operator+() const
	{
		return Vectorr(*this);
	}
	Vectorr Vectorr::operator-() const
	{
		Vectorr negativeVector(*this);
		for (int row = 0; row < height; row++) {
			negativeVector[row] = preventNegativeZero(-entries[row]);
		}
		return negativeVector;
	}

	Vectorr& Vectorr::operator=(const Vectorr& rightVector)
	{
		if (this == &rightVector) {
			return *this;
		}

		Vectorr copyRow(rightVector);
		linalg::swap(*this, copyRow);
		return *this;
	}
	void swap(Vectorr& leftVector, Vectorr& rightVector) noexcept
	{
		std::swap(leftVector.entries, rightVector.entries);
		std::swap(leftVector.height, rightVector.height);
	}

	Vectorr& Vectorr::operator+=(const Vectorr& rightVector)
	{
		if (height != rightVector.height) {
			throw std::logic_error(
				"Heights do not match : row" + std::to_string(height) + " + col" + std::to_string(rightVector.height) + "\n");
		}
		for (int row = 0; row < height; row++) {
			entries[row] = preventNegativeZero(entries[row] + rightVector[row]);
		}
		return *this;
	}
	Vectorr& Vectorr::operator-=(const Vectorr& rightVector)
	{
		if (height != rightVector.height) {
			throw std::logic_error(
				"Heights do not match : row" + std::to_string(height) + " - col" + std::to_string(rightVector.height) + "\n");
		}
		for (int row = 0; row < height; row++) {
			entries[row] = preventNegativeZero(entries[row] - rightVector[row]);
		}
		return *this;
	}
	Vectorr& Vectorr::operator*=(const double multiplier)
	{
		for (int row = 0; row < height; row++) {
			entries[row] = preventNegativeZero(multiplier * entries[row]);
		}
		return *this;
	}

	Vectorr operator+(const Vectorr& leftVector, const Vectorr& rightVector)
	{
		Vectorr resultVector(leftVector);
		resultVector += rightVector;
		return resultVector;
	}
	Vectorr operator-(const Vectorr& leftVector, const Vectorr& rightVector)
	{
		Vectorr resultVector(leftVector);
		resultVector -= rightVector;
		return resultVector;
	}
	Vectorr operator*(const double multiplier, const Vectorr& rightVector)
	{
		Vectorr resultVector(rightVector);
		resultVector *= multiplier;
		return resultVector;
	}
	Vectorr operator*(const Vectorr& leftVector, const double multiplier)
	{
		return operator*(multiplier, leftVector);
	}

	const int Vectorr::getHeight() const
	{
		return height;
	}

	const std::string Vectorr::str() const
	{
		std::string vectorString = "(" + std::to_string(height) + " row Vector)\n";
		for (int row = 0; row < height; row++) {
			vectorString += "[\t" + std::to_string(entries[row]) + "\t]\n";
		}
		return vectorString;
	}

	
	
	





	double preventNegativeZero(double value)
	{
		return (value == -0.0) ? 0.0 : value;
	}

	

	

	std::ostream& operator<<(std::ostream& outputStream, const Matrixx& outputMatrix)
	{
		outputStream << outputMatrix.str();
		return outputStream;
	}
	std::ostream& operator<<(std::ostream& outputStream, const Roww& outputRow)
	{
		outputStream << outputRow.str();
		return outputStream;
	}
	std::ostream& operator<<(std::ostream& outputStream, const Vectorr& outputVector)
	{
		outputStream << outputVector.str();
		return outputStream;
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