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




	const double Allocatable::convertNegativeZero(const double value) const
	{
		return (value == -0.0) ? 0.0 : value;
	}





	Matrixx::Matrixx(const int height, const int width)
	{
		int exceptNum = ExceptionHandler::checkValidHeight(height);
		exceptNum += ExceptionHandler::checkValidWidth(width);
		if (exceptNum > 0) {
			LengthArgument lengthArg(height, width);
			ExceptionHandler handler(ExceptionState::LengthError, exceptNum);
			handler.addArgument(lengthArg);
			handler.handleException();
		}

		this->height = height;
		this->width = width;
		rows = new Roww[height];
		for (int row = 0; row < height; row++) {
			rows[row].init(width);
		}
	}
	Matrixx::Matrixx(const Matrixx& copyMatrix)
		: Matrixx(copyMatrix.height, copyMatrix.width)
	{
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < width; col++) {
				rows[row][col] = copyMatrix[row][col];
			}
		}
	}
	Matrixx::Matrixx(const Roww& copyRow)
		: Matrixx(1, copyRow.width)
	{
		for (int col = 0; col < height; col++) {
			rows[0][col] = copyRow[col];
		}
	}
	Matrixx::Matrixx(const Vectorr& copyVector)
		: Matrixx(copyVector.height, 1)
	{
		for (int row = 0; row < height; row++) {
			rows[row][0] = copyVector[row];
		}
	}
	Matrixx::~Matrixx()
	{
		delete[] rows;
	}

	Matrixx Matrixx::block(const int beginRow, const int beginCol, const int blockHeight, const int blockWidth) const
	{
		int exceptNum = ExceptionHandler::checkRowIndex(beginRow, height);
		exceptNum += ExceptionHandler::checkColumnIndex(beginCol, width);
		if (exceptNum > 0) {
			RowIndexArgument rowIndexArg(beginRow, height);
			ColumnIndexArgument colIndexArg(beginCol, width);
			ExceptionHandler handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(rowIndexArg);
			handler.addArgument(colIndexArg);
			handler.handleException();
		}
		exceptNum = ExceptionHandler::checkRowIndex(beginRow + blockHeight - 1, height);
		exceptNum += ExceptionHandler::checkColumnIndex(beginCol + blockWidth - 1, width);
		if (exceptNum > 0) {
			RowIndexArgument rowIndexArg(beginRow + blockHeight - 1, height);
			ColumnIndexArgument colIndexArg(beginCol + blockWidth - 1, width);
			ExceptionHandler handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(rowIndexArg);
			handler.addArgument(colIndexArg);
			handler.handleException();
		}

		Matrixx blockMatrix(blockHeight, blockWidth);
		for (int row = 0; row < blockHeight; row++) {
			for (int col = 0; col < blockWidth; col++) {
				blockMatrix[row][col] = rows[beginRow + row][beginCol + col];
			}
		}
		return blockMatrix;
	}

	Matrixx Matrixx::identity(const int length)
	{
		Matrixx identityMatrix(length, length);
		for (int row = 0; row < length; row++) {
			for (int col = 0; col < length; col++) {
				identityMatrix[row][col] = (row == col) ? 1.0 : 0.0;
			}
		}
		return identityMatrix;
	}

	Matrixx Matrixx::zero(const int height, const int width)
	{
		Matrixx identityMatrix(height, width);
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < width; col++) {
				identityMatrix[row][col] = 0.0;
			}
		}
		return identityMatrix;
	}

	Roww& Matrixx::operator[](const int row)
	{
		int exceptNum = ExceptionHandler::checkRowIndex(row, height);
		if (exceptNum > 0) {
			RowIndexArgument rowIndexArg(row, height);
			ExceptionHandler handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(rowIndexArg);
			handler.handleException();
		}

		return rows[row];
	}
	const Roww& Matrixx::operator[](const int row) const
	{
		int exceptNum = ExceptionHandler::checkRowIndex(row, height);
		if (exceptNum > 0) {
			RowIndexArgument rowIndexArg(row, height);
			ExceptionHandler handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(rowIndexArg);
			handler.handleException();
		}

		return rows[row];
	}

	double& Matrixx::operator()(const int row, const int col) const
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
		if (sequence < height * width) {
			rows[sequence / width][sequence % width] = convertNegativeZero(value);
		}
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
				negativeMatrix[row][col] = convertNegativeZero(-rows[row][col]);
			}
		}
		return negativeMatrix;
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

	Matrixx& Matrixx::operator+=(const Matrixx& rightMatrix)
	{
		int exceptNum = ExceptionHandler::checkHeight(height, rightMatrix.height);
		exceptNum += ExceptionHandler::checkWidth(width, rightMatrix.width);
		if (exceptNum > 0) {
			LengthArgument leftLengthArg(height, width);
			LengthArgument rightLengthArg(rightMatrix.height, rightMatrix.width);
			OperationArgument operationArg('+', leftLengthArg, rightLengthArg);
			ExceptionHandler handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}
		
		for (int row = 0; row < height; row++) {
			rows[row] += rightMatrix[row];
		}
		return *this;
	}
	Matrixx& Matrixx::operator-=(const Matrixx& rightMatrix)
	{
		int exceptNum = ExceptionHandler::checkHeight(height, rightMatrix.height);
		exceptNum += ExceptionHandler::checkWidth(width, rightMatrix.width);
		if (exceptNum > 0) {
			LengthArgument leftLengthArg(height, width);
			LengthArgument rightLengthArg(rightMatrix.height, rightMatrix.width);
			OperationArgument operationArg('-', leftLengthArg, rightLengthArg);
			ExceptionHandler handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
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
		int exceptNum = ExceptionHandler::checkJoinLength(width, rightMatrix.height);
		if (exceptNum > 0) {
			LengthArgument leftLengthArg(height, width);
			LengthArgument rightLengthArg(rightMatrix.height, rightMatrix.width);
			OperationArgument operationArg('*', leftLengthArg, rightLengthArg);
			ExceptionHandler handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}
		
		Matrixx resultMatrix(height, rightMatrix.width);
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < rightMatrix.width; col++) {
				double product = 0.0;
				for (int join = 0; join < width; join++) {
					product += rows[row][join] * rightMatrix[join][col];
				}
				resultMatrix[row][col] = convertNegativeZero(product);
			}
		}
		linalg::swap(*this, resultMatrix);
		return *this;
	}

	Matrixx& Matrixx::operator&=(const Matrixx& rightMatrix)
	{
		int exceptNum = ExceptionHandler::checkHeight(height, rightMatrix.height);
		if (exceptNum > 0) {
			LengthArgument leftLengthArg(height, width);
			LengthArgument rightLengthArg(rightMatrix.height, rightMatrix.width);
			OperationArgument operationArg('&', leftLengthArg, rightLengthArg);
			ExceptionHandler handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
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
		int exceptNum = ExceptionHandler::checkHeight(height, rightVector.height);
		if (exceptNum > 0) {
			LengthArgument leftLengthArg(height, width);
			LengthArgument rightLengthArg(rightVector.height, 1);
			OperationArgument operationArg('&', leftLengthArg, rightLengthArg);
			ExceptionHandler handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
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
		int exceptNum = ExceptionHandler::checkWidth(width, lowerMatrix.width);
		if (exceptNum > 0) {
			LengthArgument upperLengthArg(height, width);
			LengthArgument lowerLengthArg(lowerMatrix.height, lowerMatrix.width);
			OperationArgument operationArg('|', upperLengthArg, lowerLengthArg);
			ExceptionHandler handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
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
		int exceptNum = ExceptionHandler::checkWidth(width, lowerRow.width);
		if (exceptNum > 0) {
			LengthArgument upperLengthArg(height, width);
			LengthArgument lowerLengthArg(1, lowerRow.width);
			OperationArgument operationArg('|', upperLengthArg, lowerLengthArg);
			ExceptionHandler handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
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
	Matrixx operator|(const Roww& upperRow, const Roww& lowerRow)
	{
		Matrixx appendedMatrix(upperRow);
		appendedMatrix |= lowerRow;
		return appendedMatrix;
	}
	std::ostream& operator<<(std::ostream& outputStream, const Matrixx& outputMatrix)
	{
		outputStream << outputMatrix.str();
		return outputStream;
	}
	
	const int Matrixx::getHeight() const
	{
		return height;
	}
	const int Matrixx::getWidth() const
	{
		return width;
	}

	Roww Matrixx::getRow(const int row) const
	{
		Roww copyRow(row);
		for (int col = 0; col < width; col++) {
			copyRow[col] = rows[row][col];
		}
		return copyRow;
	}
	Vectorr Matrixx::getColumn(const int col) const
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






	Roww::Roww(const int width)
	{
		init(width);
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
	void Roww::init(const int width)
	{
		int exceptNum = ExceptionHandler::checkValidWidth(width);
		if (exceptNum > 0) {
			LengthArgument lengthArg(1, width);
			ExceptionHandler handler(ExceptionState::LengthError, exceptNum);
			handler.addArgument(lengthArg);
			handler.handleException();
		}

		this->width = width;
		entries = new double[width];
	}

	double& Roww::operator[](const int col)
	{
		int exceptNum = ExceptionHandler::checkColumnIndex(col, width);
		if (exceptNum > 0) {
			ColumnIndexArgument colIndexArg(col, width);
			ExceptionHandler handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(colIndexArg);
			handler.handleException();
		}

		return entries[col];
	}
	const double& Roww::operator[](const int col) const
	{
		int exceptNum = ExceptionHandler::checkColumnIndex(col, width);
		if (exceptNum > 0) {
			ColumnIndexArgument colIndexArg(col, width);
			ExceptionHandler handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(colIndexArg);
			handler.handleException();
		}

		return entries[col];
	}

	Allocator& Roww::operator<<(const double value)
	{
		this->allocate(0, value);
		return *(new Allocator(*this, 1));
	}
	void Roww::allocate(const int sequence, const double value)
	{
		if (sequence < width) {
			entries[sequence] = convertNegativeZero(value);
		}
	}

	Roww Roww::operator+() const
	{
		return Roww(*this);
	}
	Roww Roww::operator-() const
	{
		Roww negativeRow(width);
		for (int col = 0; col < width; col++) {
			negativeRow[col] = convertNegativeZero(-entries[col]);
		}
		return negativeRow;
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

	Roww& Roww::operator+=(const Roww& rightRow)
	{
		int exceptNum = ExceptionHandler::checkWidth(width, rightRow.width);
		if (exceptNum > 0) {
			LengthArgument leftLengthArg(1, width);
			LengthArgument rightLengthArg(1, rightRow.width);
			OperationArgument operationArg('+', leftLengthArg, rightLengthArg);
			ExceptionHandler handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}
		
		for (int col = 0; col < width; col++) {
			entries[col] = convertNegativeZero(entries[col] + rightRow[col]);
		}
		return *this;
	}
	Roww& Roww::operator-=(const Roww& rightRow)
	{
		int exceptNum = ExceptionHandler::checkWidth(width, rightRow.width);
		if (exceptNum > 0) {
			LengthArgument leftLengthArg(1, width);
			LengthArgument rightLengthArg(1, rightRow.width);
			OperationArgument operationArg('-', leftLengthArg, rightLengthArg);
			ExceptionHandler handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}
		
		for (int col = 0; col < width; col++) {
			entries[col] = convertNegativeZero(entries[col] - rightRow[col]);
		}
		return *this;
	}
	Roww& Roww::operator*=(const double multiplier)
	{
		for (int col = 0; col < width; col++) {
			entries[col] = convertNegativeZero(multiplier * entries[col]);
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
	std::ostream& operator<<(std::ostream& outputStream, const Roww& outputRow)
	{
		outputStream << outputRow.str();
		return outputStream;
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

	
	


	Vectorr::Vectorr(const int height)
	{
		init(height);
	}
	Vectorr::Vectorr(const Vectorr& copyVector)
		: Vectorr(copyVector.height)
	{
		for (int row = 0; row < height; row++) {
			entries[row] = copyVector[row];
		}
	}
	Vectorr::~Vectorr()
	{
		delete[] entries;
	}
	void Vectorr::init(const int height)
	{
		int exceptNum = ExceptionHandler::checkValidHeight(height);
		if (exceptNum > 0) {
			LengthArgument lengthArg(height, 1);
			ExceptionHandler handler(ExceptionState::LengthError, exceptNum);
			handler.addArgument(lengthArg);
			handler.handleException();
		}

		this->height = height;
		entries = new double[height];
	}

	double& Vectorr::operator[](const int row)
	{
		int exceptNum = ExceptionHandler::checkRowIndex(row, height);
		if (exceptNum > 0) {
			RowIndexArgument rowIndexArg(row, height);
			ExceptionHandler handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(rowIndexArg);
			handler.handleException();
		}

		return entries[row];
	}
	const double& Vectorr::operator[](const int row) const
	{
		int exceptNum = ExceptionHandler::checkRowIndex(row, height);
		if (exceptNum > 0) {
			RowIndexArgument rowIndexArg(row, height);
			ExceptionHandler handler(ExceptionState::OutOfRange, exceptNum);
			handler.addArgument(rowIndexArg);
			handler.handleException();
		}

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
			entries[sequence] = convertNegativeZero(value);
		}
	}

	Vectorr Vectorr::operator+() const
	{
		return Vectorr(*this);
	}
	Vectorr Vectorr::operator-() const
	{
		Vectorr negativeVector(height);
		for (int row = 0; row < height; row++) {
			negativeVector[row] = convertNegativeZero(-entries[row]);
		}
		return negativeVector;
	}

	Vectorr& Vectorr::operator=(const Vectorr& rightVector)
	{
		if (this == &rightVector) {
			return *this;
		}

		Vectorr copyVector(rightVector);
		swap(*this, copyVector);
		return *this;
	}
	void swap(Vectorr& leftVector, Vectorr& rightVector) noexcept
	{
		std::swap(leftVector.entries, rightVector.entries);
		std::swap(leftVector.height, rightVector.height);
	}

	Vectorr& Vectorr::operator+=(const Vectorr& rightVector)
	{
		int exceptNum = ExceptionHandler::checkHeight(height, rightVector.height);
		if (exceptNum > 0) {
			LengthArgument leftLengthArg(height, 1);
			LengthArgument rightLengthArg(rightVector.height, 1);
			OperationArgument operationArg('+', leftLengthArg, rightLengthArg);
			ExceptionHandler handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}
		
		for (int row = 0; row < height; row++) {
			entries[row] = convertNegativeZero(entries[row] + rightVector[row]);
		}
		return *this;
	}
	Vectorr& Vectorr::operator-=(const Vectorr& rightVector)
	{
		int exceptNum = ExceptionHandler::checkHeight(height, rightVector.height);
		if (exceptNum > 0) {
			LengthArgument leftLengthArg(height, 1);
			LengthArgument rightLengthArg(rightVector.height, 1);
			OperationArgument operationArg('-', leftLengthArg, rightLengthArg);
			ExceptionHandler handler(ExceptionState::ArithmeticException, exceptNum);
			handler.addArgument(operationArg);
			handler.handleException();
		}
		
		for (int row = 0; row < height; row++) {
			entries[row] = convertNegativeZero(entries[row] - rightVector[row]);
		}
		return *this;
	}
	Vectorr& Vectorr::operator*=(const double multiplier)
	{
		for (int row = 0; row < height; row++) {
			entries[row] = convertNegativeZero(multiplier * entries[row]);
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
	std::ostream& operator<<(std::ostream& outputStream, const Vectorr& outputVector)
	{
		outputStream << outputVector.str();
		return outputStream;
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
}